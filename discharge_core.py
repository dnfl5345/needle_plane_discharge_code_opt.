# -*- coding: utf-8 -*-
"""
==============================================================================
 고전압공학 방전현상 시뮬레이터 v18 - 물리 코어 (discharge_core.py)
==============================================================================
 v17 물리 모델(전극 구성, 4성분 기체, 반응률 테이블, 라디칼 화학, 개시통계 등)은
 그대로 유지하고 '계산 속도'만 개선한 코어입니다. GUI(needle_plane_discharge_hmi.py)
 와 헤드리스 배치 실행기(run_batch.py)가 공통으로 import 합니다.

 [v18 변경사항 - 속도]
  1) 포아송 방정식: numpy np.roll 기반 적-흑 SOR(서브스텝당 60회, 미수렴)을
     희소행렬 직접해법(scipy SuperLU, 기하 고정 -> LU 분해 1회, 매 스텝 전진/후진
     대입만)으로 교체. v17 SOR 의 고정점(수렴해)과 '동일한 선형계'를 풀므로
     결과는 v17 을 무한히 반복시킨 것과 같고, 서브스텝당 65 ms -> ~1 ms.
     (scipy 가 없으면 v17 SOR 로 자동 폴백: poisson_used == "sor")
  2) 반응률 보간: 서브스텝마다 15회 반복되던 log-log np.interp(각 log10/clip/10**)
     을 균일 log(E/N) 격자의 조밀 LUT(8000점) 1회 인덱싱 + 종별 gather 로 교체.
     조성 고정이므로 nu_ion = sum_X k_X n_X 등 '합성 테이블'을 초기화 시 1회 생성.
  3) 유전완화 시간 제한: v17 은 tau_d = eps0/(q·max(ne)·max(mu)) 로 '격자 최대
     이동도(저전계 0.1 m2/Vs)'와 '최대 밀도'를 곱해 과도하게 작게 평가.
     v18 기본(dt_rule="v18")은 물리적으로 올바른 tau_d = eps0/max(sigma),
     sigma = q(mu_e n_e + mu_+ n_+ + mu_- n_-) 를 사용 (dt 약 2~3배 확대).
     dt_rule="v17" 로 두면 v17 과 동일한 제한식을 사용 (검증용).
  4) 초기화: 2500+2000+400 회 SOR(6.5 s) -> 직접해법 3회(<0.1 s).
  5) run_until(t_end, ...) 헤드리스 실행 API 및 벽시계 계측(t_wall_*) 추가.
  6) 갭 브리징 이후의 dt 붕괴 방지 (v17 에서 '무한히 느려지는' 직접 원인)
     - v17: 스트리머가 대향전극에 닿은 뒤 밀도가 상한(1e22)에 걸리면 종별 독립 clip
       으로 인위적 순전하 -> |E| 가 1e6 kV/cm 대로 발산 -> CFL dt 가 fs 대로 붕괴.
     - v18: 전리 생성량을 상한 내로 제한하고 상한 초과 시 세 종을 같은 비율로 축소해
       전하중성을 보존(_neutral_cap). 축상 ne 가 전 구간 문턱(1e18) 초과 시
       '브리징' 판정(bridged, t_bridge) -> 배치/GUI 가 자동으로 잔광 단계 전환 가능.
       dt 가 dt_floor(0.01 ps) 아래로 20 서브스텝 연속이면 stalled 플래그.
  7) 외부 회로 직렬저항 R_series (선택, 기본 0 = v17 과 동일한 이상 전압원)
     - Sato 공식 방전전류 I(t) = ∫ σ E·E_L1 dV (축대칭 환산) 를 매 스텝 계산하고
       V_needle = V_src(t) - R·I 로 전극전압을 갱신 -> 브리징 후 전압붕괴(스파크 전이)
       가 물리적으로 재현되어 임펄스 전 구간(300 ns) 실행이 가능.
     - circ_hist 에 (t, V_src, V_needle, I) 이력 기록. 변위전류(C dV/dt)는 미포함.
  8) 적응형 반암시적 포아송 (semi_implicit=True, 기본)
     - 브리징 후 도전 채널(ne ~ 1e21-1e22)에서는 유전완화 한계 dt < 0.5 ε0/σ 가
       0.03-0.1 ps 로 붕괴 -> 명시적 스킴으로는 300 ns 에 수 시간~수십 시간.
     - 유전완화 한계가 다른 한계(CFL/전리/파형)보다 'LU 재분해 비용비' 이상 작을 때만
       ∇·[(ε_r ε0 + dt σ)∇V] = -ρ 를 매 스텝 재조립·재분해(≈40 ms @161x125)하여 풀고
       dt 는 CFL/전리 한계로 결정 (Ventzek 1994 / Hagelaar-Kroesen 2000 형식).
       그 외 구간은 고정 LU 명시적 스킴(≈1 ms) 유지. dt_hist 제한요인 "semi:..." 표기.
  9) (버그 수정) BOLSIG+ CSV 로 일부 반응률만 교체할 때 CSV 의 E/N 격자 길이가
     기본 테이블(13점)과 다르면 v17 은 np.interp 길이 불일치로 예외 -> v18 은
     교체되지 않은 기본 테이블을 새 E/N 격자로 log-log 재보간.

 [실행]  pip install numpy scipy matplotlib
==============================================================================
"""

import time
import numpy as np

try:
    import scipy.sparse as _sp
    import scipy.sparse.linalg as _spla
    HAVE_SCIPY = True
except Exception:                                   # pragma: no cover
    HAVE_SCIPY = False

EPS0 = 8.8541878128e-12
QE   = 1.602176634e-19
N_LOSCHMIDT = 2.446e25
TD   = 1e-21                       # 1 Townsend = 1e-21 V·m^2

GAS_AB = {"AIR": (15.0, 365.0), "AR": (12.0, 180.0), "H2O": (13.0, 290.0)}

PLATE_PRESET = {
    "Metal (grounded)": ("GND_METAL",  1.0, "#707070"),
    "Metal (floating)": ("FLOAT_METAL", 1.0, "#8a8a8a"),
    "Quartz":           ("DIELECTRIC", 3.8, "#a6cee3"),
    "Ceramic":          ("DIELECTRIC", 9.5, "#e0c9a6"),
    "Custom (enter εr)": ("DIELECTRIC", None, "#c5b0d5"),
}

CONFIG_MAP = {"Needle - Plane":          "NEEDLE_PLANE",
              "Needle - Needle":         "NEEDLE_NEEDLE",
              "3 Needles - Plane":       "THREE_NEEDLE",
              "Needle - Rogowski plane": "NEEDLE_ROGOWSKI",
              "Rogowski - Rogowski":     "PLANE_PLANE",
              "Sphere - Sphere":         "SPHERE_SPHERE"}

# =============================================================================
#  단면적 기반 반응률 테이블  k(E/N) [m^3/s], 전자이동도 muN [1/(m·V·s)]
#  * BOLSIG+ (Phelps: N2, O2, Ar / Itikawa: H2O 단면적) Boltzmann 해석 결과의
#    대표값(교육/경향분석용). 정량 연구 시 자체 BOLSIG+ 출력 CSV로 교체할 것.
#  * E/N 범위 밖은 경계값으로 고정(외삽 없음). log-log 선형보간.
# =============================================================================
RATE_TABLES_DEFAULT = {
    "EN":         [1.0,    2.0,    5.0,    10.0,   20.0,   50.0,   100.0,  150.0,  200.0,  300.0,  500.0,  800.0,  1200.0],
    # 전리 (ionization)
    "k_ion_N2":   [1e-30,  1e-30,  1e-30,  1e-29,  1e-26,  3e-21,  1.5e-18, 3e-17, 1.8e-16, 9e-16, 3.5e-15, 8e-15, 1.4e-14],
    "k_ion_O2":   [1e-30,  1e-30,  1e-29,  1e-27,  1e-24,  3e-20,  8e-18,  9e-17,  4e-16,  1.5e-15, 4.5e-15, 9e-15, 1.5e-14],
    "k_ion_Ar":   [1e-30,  1e-30,  1e-30,  1e-28,  1e-25,  5e-20,  3e-17,  3e-16,  1e-15,  3e-15,  8e-15,  1.4e-14, 2e-14],
    "k_ion_H2O":  [1e-30,  1e-30,  1e-30,  1e-28,  1e-24,  1e-20,  3e-18,  4e-17,  1.5e-16, 7e-16, 2.5e-15, 6e-15, 1e-14],
    # 부착 (attachment: O2 해리성 + 3체 유효, H2O 해리성)
    "k_att_O2":   [1.5e-16, 1.3e-16, 1.0e-16, 8e-17, 6e-17, 5e-17, 6e-17,  8e-17,  9e-17,  1.0e-16, 1.1e-16, 1.2e-16, 1.2e-16],
    "k_att_H2O":  [2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17,  2e-17],
    # 해리 (dissociation -> 라디칼 생성)
    "k_diss_O2_6eV": [1e-30, 1e-30, 8e-23, 8e-20, 2.0e-17, 2.0e-16, 5.0e-16, 7.0e-16, 9.0e-16, 1.2e-15, 1.5e-15, 1.7e-15, 1.8e-15],
    "k_diss_O2_8eV": [1e-30, 1e-30, 2e-23, 2e-20, 1.0e-17, 2.0e-16, 7.0e-16, 1.3e-15, 2.1e-15, 3.3e-15, 4.5e-15, 5.3e-15, 6.2e-15],
    "k_diss_N2":  [1e-30,  1e-30,  1e-30,  1e-26,  1e-21,  1e-17,  2e-16,  6e-16,  1.2e-15, 2.5e-15, 5e-15, 8e-15,  1.1e-14],
    "k_diss_H2O": [1e-30,  1e-30,  1e-24,  1e-20,  1e-17,  2e-16,  8e-16,  1.4e-15, 2e-15, 3e-15,  4.5e-15, 6e-15, 7e-15],
    # O2- 충돌이탈 (O2- + M -> e + O2 + M), ~100Td에서 nu_det ~ 1e8 /s 오더
    "k_det_O2m":  [1e-24,  3e-24,  1e-23,  1e-22,  3e-21,  3e-19,  4e-18,  1.2e-17, 2.5e-17, 5e-17,  1e-16,  1.5e-16, 2e-16],
    # 여기 (N2 C3Πu, 11.0 eV -> SPS 발광원)
    "k_exc_N2":   [1e-30,  1e-30,  1e-28,  1e-24,  1e-20,  3e-17,  4e-16,  1.5e-15, 3e-15,  6e-15,  1.0e-14, 1.4e-14, 1.7e-14],
    # 전자 환산이동도 mu_e·N
    "muN":        [2.4e24, 2.2e24, 2.0e24, 1.8e24, 1.6e24, 1.4e24, 1.2e24, 1.15e24, 1.1e24, 1.0e24, 9e23,  8e23,   7e23],
}

# 활성종 화학상수 [m^3/s, 3체는 m^6/s] (300 K 문헌 대표값)
CHEM = {
    "kq_OD_N2": 3.1e-17,   # O(1D)+N2 -> O(3P)+N2   (소광)
    "kq_OD_O2": 4.0e-17,   # O(1D)+O2 -> O(3P)+O2   (소광)
    "k_OD_H2O": 2.2e-16,   # O(1D)+H2O -> 2OH
    "k_ND_O2":  6.0e-18,   # N(2D)+O2 -> NO+O(3P)
    "kq_ND_N2": 1.7e-20,   # N(2D)+N2 -> N(4S)+N2   (소광)
    "k_NS_O2":  9.0e-23,   # N(4S)+O2 -> NO+O(3P)   (상온: 매우 느림)
    "k_NS_NO":  3.0e-17,   # N(4S)+NO -> N2+O(3P)
    "k_NS_OH":  4.7e-17,   # N(4S)+OH -> NO+H
    "k_OP_OH":  3.3e-17,   # O(3P)+OH -> O2+H
    "k_OH_OH":  2.0e-18,   # OH+OH -> 소멸
    "k3_O_O2":  6.0e-46,   # O(3P)+O2+M -> O3+M
    "k3_O_O":   5.0e-46,   # O(3P)+O(3P)+M -> O2+M
    "k3_H_O2":  6.0e-44,   # H+O2+M -> HO2+M
}
RAD_SPECIES = ["O(P)", "O(D)", "N(S)", "N(D)", "OH", "NO", "O3", "H"]
RAD_COLOR = {"O(P)": "tab:blue", "O(D)": "tab:cyan",
             "N(S)": "tab:orange", "N(D)": "tab:brown",
             "OH": "tab:red", "NO": "tab:green", "O3": "m", "H": "k"}
# (색, 선종, 마커): 흑백에서도 종 구분 가능
RAD_STYLE = {"O(P)": ("tab:blue",  "-",  "o"), "O(D)": ("tab:cyan",  "--", "s"),
             "N(S)": ("tab:orange","-.", "^"), "N(D)": ("tab:brown", ":",  "v"),
             "OH":   ("tab:red",   (0,(5,1)),      "D"),
             "NO":   ("tab:green", (0,(3,1,1,1)),  "P"),
             "O3":   ("m",         (0,(1,1)),      "X"),
             "H":    ("k",         (0,(5,2,1,2)),  "*")}


def parse_rate_csv(path):
    """BOLSIG+ 반응률 CSV 파싱 -> (override_dict, 교체된 반응률 이름 목록)"""
    import csv as _csv
    with open(path, encoding="utf-8-sig") as f:
        rows = list(_csv.reader(f))
    head = [h.strip() for h in rows[0]]
    if not head or head[0] != "EN_Td":
        raise ValueError("First column header must be EN_Td.")
    data = np.array([[float(v) for v in r] for r in rows[1:] if r])
    ov = {"EN": data[:, 0]}
    loaded = []
    for c, name in enumerate(head[1:], start=1):
        if name in RATE_TABLES_DEFAULT:
            ov[name] = data[:, c]
            loaded.append(name)
        elif name == "k_diss_O2":              # 구버전 열: 채널 분배(하위호환)
            ov["k_diss_O2_6eV"] = 0.4 * data[:, c]
            ov["k_diss_O2_8eV"] = 0.6 * data[:, c]
            loaded += ["k_diss_O2_6eV", "k_diss_O2_8eV"]
    if not loaded:
        raise ValueError("No recognizable rate columns.\nAllowed: " +
                         ", ".join(k for k in RATE_TABLES_DEFAULT if k != "EN"))
    return ov, loaded


# =============================================================================
#  반응률 LUT: 균일 log10(E/N) 격자에 log-log 선형보간한 조밀 테이블
#  * 서브스텝마다 index() 1회(log10+clip) 후 종별 get()은 gather+FMA 만 수행
#  * v17 _k() (13점 log-log 보간)과의 차이: 격자 간격 ~0.0008 decade -> 상대오차 <1e-4
# =============================================================================
class RateLUT:
    def __init__(self, RT, npts=8000):
        EN = np.asarray(RT["EN"], float)
        lEN = np.log10(EN)
        self.lo, self.hi = float(EN[0]), float(EN[-1])
        self.l0 = float(lEN[0])
        self.n = int(npts)
        self.dl = (float(lEN[-1]) - self.l0) / (self.n - 1)
        lg = self.l0 + self.dl * np.arange(self.n)
        self.tab = {}
        for k, v in RT.items():
            if k == "EN":
                continue
            v = np.asarray(v, float)
            self.tab[k] = 10.0 ** np.interp(lg, lEN, np.log10(np.maximum(v, 1e-30)))

    def add_combo(self, name, parts):
        """합성 테이블: name = sum_i w_i * tab[p_i]  (parts: [(p_i, w_i), ...])"""
        self.tab[name] = sum(w * self.tab[p] for p, w in parts)

    def index(self, EN):
        x = (np.log10(np.clip(EN, self.lo, self.hi)) - self.l0) / self.dl
        i0 = x.astype(np.int64)                      # x >= 0 -> floor
        np.clip(i0, 0, self.n - 2, out=i0)
        return i0, x - i0

    def get(self, name, i0, f):
        t = self.tab[name]
        t0 = t[i0]
        return t0 + (t[i0 + 1] - t0) * f


# =============================================================================
#  물리 시뮬레이션 코어
# =============================================================================
class NeedlePlaneDischarge:
    def __init__(self, gap_mm=10.0, voltage_kV=30.0, polarity=+1,
                 pressure_atm=1.0,
                 r_main_mm=1.0, r_tip_mm=0.15, taper_mm=3.0,
                 nx=161, n_background=1e13, seed_density=5e18,
                 wave_mode="IMPULSE", t_front_us=0.01, t_tail_us=0.1,
                 n2_pct=78.0, o2_pct=21.0, ar_pct=1.0, h2o_pct=0.0,
                 gamma_see=0.05,
                 dt_user_ns=0.0, snap_ns=5.0,
                 plate_mode="GND_METAL", plate_thick_mm=0.0, plate_epsr=3.8,
                 rate_model="TABLE", rate_override=None,
                 cath_J0=1e-2, n_ion_bg=1e9, q_bg=1e7,
                 config="NEEDLE_PLANE", needle_sep_mm=4.0, sphere_R_mm=5.0,
                 rog_halfw_mm=6.0, rog_edge_mm=2.0,
                 poisson="direct", dt_rule="v18", lut_pts=8000,
                 R_series_ohm=0.0, n_cap=1e22, dt_floor_ps=0.01,
                 semi_implicit=True):

        # ---------- 전기/기체 ----------
        self.gap      = gap_mm * 1e-3
        self.Vpeak    = polarity * abs(voltage_kV) * 1e3
        self.polarity = polarity
        self.p_atm    = max(pressure_atm, 0.05)
        self.p_torr   = 760.0 * self.p_atm
        self.Ngas     = N_LOSCHMIDT * self.p_atm

        # 조성 입력: N2/O2/Ar/H2O [%] -> 합으로 자동 정규화 (몰분율)
        raw = np.array([max(n2_pct, 0.0), max(o2_pct, 0.0),
                        max(ar_pct, 0.0), max(h2o_pct, 0.0)], float)
        tot = raw.sum()
        if tot <= 0.0:
            raw, tot = np.array([78.0, 21.0, 1.0, 0.0]), 100.0
        self.x_n2, self.x_o2, self.x_ar, self.x_h2o = (raw / tot).tolist()
        self.nN2  = self.x_n2  * self.Ngas
        self.nO2  = self.x_o2  * self.Ngas
        self.nAr  = self.x_ar  * self.Ngas
        self.nH2O = self.x_h2o * self.Ngas
        # Townsend 근사 모드용 등가 분율
        self.f_air = self.x_n2 + self.x_o2
        self.f_ar  = self.x_ar
        self.f_h2o = self.x_h2o

        # 반응률 테이블 (기본값 + 사용자 CSV override 병합)
        self.rate_model = rate_model              # "TABLE" or "TOWNSEND"
        self.RT = {k: np.array(v, float) for k, v in RATE_TABLES_DEFAULT.items()}
        if rate_override:
            ov = {k: np.array(v, float) for k, v in rate_override.items()}
            if "EN" in ov and (ov["EN"].shape != self.RT["EN"].shape or
                               not np.allclose(ov["EN"], self.RT["EN"])):
                # (v18 수정) 교체되지 않은 기본 테이블을 새 E/N 격자로 log-log 재보간
                lE_old = np.log10(self.RT["EN"]); lE_new = np.log10(ov["EN"])
                for k in list(self.RT):
                    if k != "EN" and k not in ov:
                        self.RT[k] = 10.0 ** np.interp(
                            lE_new, lE_old, np.log10(np.maximum(self.RT[k], 1e-30)))
            for k, v in ov.items():
                self.RT[k] = v
        self._logEN = np.log10(self.RT["EN"])
        self.LUT = RateLUT(self.RT, npts=lut_pts)
        L = self.LUT
        L.add_combo("nu_i", [("k_ion_N2", self.nN2), ("k_ion_O2", self.nO2),
                             ("k_ion_Ar", self.nAr), ("k_ion_H2O", self.nH2O)])
        L.add_combo("nu_a", [("k_att_O2", self.nO2), ("k_att_H2O", self.nH2O)])
        L.add_combo("nu_det", [("k_det_O2m", self.Ngas)])
        L.add_combo("nu_d6", [("k_diss_O2_6eV", self.nO2)])
        L.add_combo("nu_d8", [("k_diss_O2_8eV", self.nO2)])
        L.add_combo("nu_dN2", [("k_diss_N2", self.nN2)])
        L.add_combo("nu_dH2O", [("k_diss_H2O", self.nH2O)])
        L.add_combo("nu_exc", [("k_exc_N2", self.nN2)])
        L.add_combo("mu", [("muN", 1.0 / self.Ngas)])

        self.gamma_see = max(gamma_see, 0.0)
        self.dt_user   = max(dt_user_ns, 0.0) * 1e-9
        self.snap_dt   = max(snap_ns, 0.1) * 1e-9
        self.next_snap = 0.0
        self.profiles  = []
        self.streak    = []      # [(t, 축상 ne 프로파일)]  z-t 스트릭용
        self.streak_em = []      # [(t, 축상 발광강도)]      발광 스트릭용
        self.rad_hist  = []      # [(t, {종: 총개수})]      활성종 이력용

        self.afterglow = False                    # 잔광 모드 플래그

        # ---------- 수치 옵션 / 계측 (v18) ----------
        self.dt_rule = dt_rule if dt_rule in ("v17", "v18") else "v18"
        self.poisson_req = poisson
        self.poisson_used = "sor"
        self._lu = None
        self.n_substeps = 0
        self.t_wall_field = 0.0                   # 포아송 해법 누적 벽시계 [s]
        self.t_wall_step  = 0.0                   # step() 누적 벽시계 [s]
        self.dt_hist = []                         # (t, dt, 제한요인) 표본
        self.n_cap = max(float(n_cap), 1e15)      # 밀도 상한 (v17: 1e22 고정)
        self.dt_floor = max(float(dt_floor_ps), 0.0) * 1e-12
        self.R_series = max(float(R_series_ohm), 0.0)   # 직렬 회로저항 [Ohm] (0: 이상 전압원=v17)
        self.I_dis = 0.0                          # 방전전류 (Sato 공식, 축대칭 환산) [A]
        self.V_src_now = 0.0                      # 회로 전원전압 V_src(t)
        self.circ_hist = []                       # [(t, V_src, V_needle, I)] 전류/전압 이력
        self.t_bridge = None                      # 갭 브리징(선단이 대향전극 도달) 시각
        self.bridged = False
        self.stalled = False                      # dt 가 dt_floor 아래로 붕괴
        self._n_stall = 0
        self.n_bridge_thr = 1e18                  # 브리징 판정 축상 ne 문턱 [m^-3]
        self.semi_implicit = bool(semi_implicit)  # 유전완화 한계가 지배할 때 반암시적 포아송 사용
        self.n_semi_steps = 0                     # 반암시적으로 푼 서브스텝 수
        self.t_wall_semi = 0.0                    # 반암시적 조립+분해 누적 [s]
        self._dt_prev_semi = None
        self.t_wall_factor = 0.0
        self.semi_reuse_tol = 0.10                # 반암시적 LU 재사용 허용 ε_eff 상대변화
        self._semi_lu = None                      # (lu, bN, bP, dt*sigma 배열)
        self.n_semi_factor = 0                    # 반암시적 재분해 횟수
        self.G_cond = 0.0                         # 채널 컨덕턴스 [S] (Sato 선형화)

        # ---------- 펄스 파형 ----------
        self.wave_mode = wave_mode
        self.T_front = max(t_front_us, 1e-4) * 1e-6
        self.T_tail  = max(t_tail_us, 2 * t_front_us) * 1e-6
        a = 1.0 / (1.4 * self.T_tail)
        b = 1.0 / (0.33 * self.T_front)
        if b <= 1.5 * a:
            b = 50.0 * a
        tp = np.log(b / a) / (b - a)
        self._wa, self._wb = a, b
        self._wk = 1.0 / (np.exp(-a * tp) - np.exp(-b * tp))
        self.t_peak = tp

        # ---------- 바늘/평판 기하 ----------
        self.r_main = max(r_main_mm * 1e-3, 5e-5)
        self.rtip   = min(max(r_tip_mm * 1e-3, 1e-5), self.r_main)
        self.taper  = max(taper_mm * 1e-3, 2 * self.rtip)
        self.needle_len = max(0.45 * self.gap, self.taper + 4 * self.rtip)

        self.plate_mode = plate_mode
        self.plate_epsr = max(plate_epsr, 1.0)
        t_req = max(plate_thick_mm, 0.0) * 1e-3

        # ---------- 전극 구성 ----------
        self.config = str(config)
        self.d_sep = max(needle_sep_mm, 0.5) * 1e-3     # 3침 간격
        self.R_sph = max(sphere_R_mm, 0.5) * 1e-3       # 구 반경
        self.w_rog = max(rog_halfw_mm, 1.0) * 1e-3      # Rogowski 평탄 반폭
        self.T_rog = max(rog_edge_mm, 0.5) * 1e-3       # Rogowski 연부 반경(두께)
        if self.config in ("NEEDLE_NEEDLE", "SPHERE_SPHERE",
                           "NEEDLE_ROGOWSKI", "PLANE_PLANE"):
            self.plate_mode = "GND_METAL"               # 하부 형상 금속 강제
            t_req = 0.0

        self.width = 2.0 * self.gap
        if self.config == "THREE_NEEDLE":
            self.width = max(self.width, 2.0*(self.d_sep + self.r_main) + self.gap)
        elif self.config == "SPHERE_SPHERE":
            self.width = max(self.width, 4.0*self.R_sph + self.gap)
        elif self.config in ("NEEDLE_ROGOWSKI", "PLANE_PLANE"):
            self.width = max(self.width,
                             2.0*(self.w_rog + 4.0*self.T_rog) + self.gap)
        self.nx = int(nx)
        self.h  = self.width / (self.nx - 1)
        self.n_solid = int(round(t_req / self.h))
        if self.plate_mode in ("FLOAT_METAL", "DIELECTRIC"):
            self.n_solid = max(self.n_solid, 1)
        self.t_solid = self.n_solid * self.h
        self.j_s = max(self.n_solid, 1)

        # 하부 기준면(y_bot: 축상 하부전극 표면)과 상부 구조 높이
        if self.config == "NEEDLE_NEEDLE":
            self.y_bot = self.taper + 3.0e-3
        elif self.config == "SPHERE_SPHERE":
            self.y_bot = 2.0 * self.R_sph
        elif self.config in ("NEEDLE_ROGOWSKI", "PLANE_PLANE"):
            self.y_bot = self.T_rog
        else:
            self.y_bot = self.t_solid
        if self.config == "PLANE_PLANE":
            top_len = self.T_rog + 2.0e-3
        elif self.config == "SPHERE_SPHERE":
            top_len = 2.0 * self.R_sph + 2.0e-3
        else:
            top_len = self.needle_len
        self.height = self.y_bot + self.gap + top_len
        self.ny = int(round(self.height / self.h)) + 1
        self.height = (self.ny - 1) * self.h

        self.x = np.linspace(-self.width / 2, self.width / 2, self.nx)
        self.y = np.linspace(0.0, self.height, self.ny)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.y_tip = self.y_bot + self.gap

        self.solid = np.zeros((self.ny, self.nx), bool)
        if self.n_solid >= 1:
            self.solid[:self.n_solid, :] = True

        def _top_needle(xc):
            ycn = self.y_tip + self.rtip
            frac = np.clip((self.Y - ycn) / self.taper, 0.0, 1.0)
            r_prof = self.rtip + (self.r_main - self.rtip) * frac
            body = (self.Y >= ycn) & (np.abs(self.X - xc) <= r_prof)
            capm = ((self.X - xc)**2 + (self.Y - ycn)**2) <= self.rtip**2
            axis = (np.abs(self.X - xc) <= 0.51*self.h) & (self.Y >= self.y_tip)
            return body | capm | axis

        def _bot_needle(xc):
            ycn = self.y_bot - self.rtip
            frac = np.clip((ycn - self.Y) / self.taper, 0.0, 1.0)
            r_prof = self.rtip + (self.r_main - self.rtip) * frac
            body = (self.Y <= ycn) & (np.abs(self.X - xc) <= r_prof)
            capm = ((self.X - xc)**2 + (self.Y - ycn)**2) <= self.rtip**2
            axis = (np.abs(self.X - xc) <= 0.51*self.h) & (self.Y <= self.y_bot)
            return body | capm | axis

        self.gnd = np.zeros((self.ny, self.nx), bool)   # 하부 '형상' 접지전극
        cfg = self.config
        if cfg == "THREE_NEEDLE":
            self.needle = (_top_needle(-self.d_sep) | _top_needle(0.0) |
                           _top_needle(+self.d_sep))
        elif cfg == "NEEDLE_NEEDLE":
            self.needle = _top_needle(0.0)
            self.gnd = _bot_needle(0.0)
        elif cfg == "NEEDLE_ROGOWSKI":
            self.needle = _top_needle(0.0)
            sdist = np.abs(self.X) - self.w_rog
            Lr = 2.0 * self.T_rog                       # 램프 특성폭
            prof = np.where(sdist <= 0.0, self.T_rog,
                            np.where(sdist <= 2.0*Lr,
                                     self.T_rog*np.exp(-(sdist/Lr)**2), 0.0))
            self.gnd = self.Y <= prof + 1e-12
        elif cfg == "PLANE_PLANE":                      # Rogowski - Rogowski 쌍
            sdist = np.abs(self.X) - self.w_rog
            Lr = 2.0 * self.T_rog                       # 램프 특성폭
            prof = np.where(sdist <= 0.0, self.T_rog,
                            np.where(sdist <= 2.0*Lr,
                                     self.T_rog*np.exp(-(sdist/Lr)**2), 0.0))
            self.gnd = self.Y <= prof + 1e-12           # 하부 접지 Rogowski
            Yp = (self.y_tip + self.T_rog) - self.Y     # 상부 미러 좌표
            self.needle = ((Yp >= -1e-12) & (Yp <= prof + 1e-12)) | \
                          (self.Y >= self.y_tip + self.T_rog)
        elif cfg == "SPHERE_SPHERE":
            yct = self.y_tip + self.R_sph
            ycb = self.y_bot - self.R_sph
            stem = max(self.rtip, 1.01*self.h)
            self.needle = ((self.X**2 + (self.Y - yct)**2) <= self.R_sph**2) | \
                          ((np.abs(self.X) <= stem) & (self.Y >= yct))
            self.gnd    = ((self.X**2 + (self.Y - ycb)**2) <= self.R_sph**2) | \
                          ((np.abs(self.X) <= stem) & (self.Y <= ycb))
        else:                                            # NEEDLE_PLANE
            self.needle = _top_needle(0.0)

        def _adj(mask):
            nb = np.zeros_like(mask)
            nb[1:, :]  |= mask[:-1, :]
            nb[:-1, :] |= mask[1:, :]
            nb[:, 1:]  |= mask[:, :-1]
            nb[:, :-1] |= mask[:, 1:]
            return nb & (~mask)

        self.needle_adj = _adj(self.needle) & (~self.gnd) & (~self.solid)
        self.gnd_adj    = _adj(self.gnd) & (~self.needle) & (~self.solid)
        if self.gnd.any():
            self.bot_adj = self.gnd_adj
        else:
            self.bot_adj = np.zeros_like(self.gnd)
            self.bot_adj[self.j_s, :] = ~self.needle[self.j_s, :]
        icc = self.nx // 2
        topj = np.where(self.needle[:, icc])[0]
        self.j_ax_hi = int(topj.min() - 1) if topj.size else self.ny - 2
        botj = np.where(self.gnd[:, icc] | self.solid[:, icc])[0]
        self.j_ax_lo = int(botj.max() + 1) if botj.size else self.j_s

        # ---------- 유전율/SOR 계수 ----------
        self.epsr = np.ones((self.ny, self.nx))
        if self.plate_mode == "DIELECTRIC":
            self.epsr[self.solid] = self.plate_epsr
        e = self.epsr
        def harm(a_, b_): return 2.0 * a_ * b_ / (a_ + b_)
        self.cN = np.ones_like(e); self.cS = np.ones_like(e)
        self.cE = np.ones_like(e); self.cW = np.ones_like(e)
        self.cN[:-1, :] = harm(e[:-1, :], e[1:, :]);  self.cN[-1, :] = e[-1, :]
        self.cS[1:, :]  = harm(e[1:, :],  e[:-1, :]); self.cS[0, :]  = e[0, :]
        self.cE[:, :-1] = harm(e[:, :-1], e[:, 1:]);  self.cE[:, -1] = e[:, -1]
        self.cW[:, 1:]  = harm(e[:, 1:],  e[:, :-1]); self.cW[:, 0]  = e[:, 0]
        self.csum = self.cN + self.cS + self.cE + self.cW

        # ---------- 수송/반응 계수 ----------
        self.mu_e  = 0.048 / self.p_atm          # Townsend 모드용 / 안정한계 참조
        self.De    = 0.18  / self.p_atm
        self.mu_ip = 2.0e-4 / self.p_atm
        self.mu_in = 2.2e-4 / self.p_atm
        self.beta  = 2.0e-13
        self.D_rad = 3.0e-5 / self.p_atm         # 라디칼 확산계수 (O in air 오더)
        self.k_O3_2b = CHEM["k3_O_O2"] * self.Ngas   # 3체 -> 유효 2체
        self.k_OO_2b = CHEM["k3_O_O"]  * self.Ngas

        # ---------- 상태 변수 ----------
        self.V     = np.zeros((self.ny, self.nx))
        self.sigma = np.zeros(self.nx)
        self.Q_acc = 0.0
        self.V_plate = 0.0

        self.J0 = max(cath_J0, 0.0)          # 음극 배경 전자방출 [A/m^2]
        self.n_ion_bg = max(n_ion_bg, 0.0)   # 배경 소이온 밀도 (O2-, +이온쌍) [m^-3]
        self.q_bg = max(q_bg, 0.0)           # 배경 이온쌍 생성률 [m^-3 s^-1]
        self.init_stats = []                 # [(t, Lambda, H, F)] 개시통계 이력
        self._H_hazard = 0.0
        self._t_prev_frame = 0.0
        self.M_townsend = 1.0                # e^{int alpha dz} (축상)
        self.gamma_MK = 0.0                  # gamma*(M-1): >=1 자립
        self.n_bg = max(n_background, 0.0)   # 0 허용 (Townsend 개시 모드)
        self.n_e  = np.full_like(self.V, self.n_bg)
        self.n_ip = np.full_like(self.V, max(self.n_bg, self.n_ion_bg))
        self.n_in = np.full_like(self.V, self.n_ion_bg)   # 배경 O2- 저장고
        self.rad = {sp: np.zeros_like(self.V) for sp in RAD_SPECIES}

        sig = max(2.5 * self.h, 0.12e-3)
        g = seed_density * np.exp(-((self.X**2 + (self.Y - self.y_tip)**2) / (2 * sig**2)))
        self.n_e  += g
        self.n_ip += g
        self._zero_in_solids()

        jj, ii = np.meshgrid(np.arange(self.ny), np.arange(self.nx), indexing="ij")
        interior = np.ones_like(self.V, bool)
        interior[0, :] = interior[-1, :] = False
        interior[:, 0] = interior[:, -1] = False
        interior &= ~self.needle
        if self.plate_mode in ("GND_METAL", "FLOAT_METAL"):
            interior &= ~self.solid
        self.red   = interior & (((jj + ii) % 2) == 0)
        self.black = interior & (((jj + ii) % 2) == 1)
        self._gas_mask = ~self.solid
        self._vol_w = np.pi * np.abs(self.X) * self.h * self.h   # 축대칭 부피가중

        # ---------- (v18) 포아송 직접해법 준비 ----------
        if self.poisson_req == "direct" and HAVE_SCIPY:
            self._build_poisson()

        self.t = 0.0
        self.dt_last = 0.0
        self.dt_limited = False
        self.dt_limiter = "-"
        self.V_src_now = self.applied_voltage(0.0)
        self.Va_now = self.V_src_now
        self._compute_laplace_field()
        self.V[:] = self.Va_now * np.clip((self.Y - self.y_bot) /
                                          max(self.height - self.y_bot, 1e-9), 0, 1)
        self.V[self.gnd] = 0.0
        self._apply_bc()
        self.solve_field(iters=2500, w=1.9)

        Vsave, va = self.V.copy(), self.Va_now
        self.Va_now = self.Vpeak
        self._apply_bc(); self.solve_field(iters=2000, w=1.9)
        self.E0max_peak = np.nanmax(self.Emag_masked())
        self.Va_now = va; self.V[:] = Vsave
        self._apply_bc(); self.solve_field(iters=400)
        self._record_profile()
        self._record_hist()

    # ------------------------------------------------------------------
    def _build_poisson(self):
        """v17 적-흑 SOR 의 고정점과 동일한 선형계  csum·V - Σc·V_nb = h²ρ/ε0  를
        미지수(내부 기체/유전체 셀)에 대해 희소행렬로 조립하고 LU 분해(기하 고정 -> 1회).
        Dirichlet 이웃(침=Va, 접지형상=0, 금속평판=0/V_plate, 바닥행)은 우변 계수
        벡터(bN: Va 배수, bP: V_plate 배수)로, Neumann 경계(측면/상단 미러)는 대각항
        감소로 처리한다. 이웃 분류(기하)는 1회만 수행하고 _assemble_poisson 이 재사용
        (반암시적 모드에서는 계수만 바꿔 매 스텝 재조립)."""
        ny, nx = self.ny, self.nx
        unk = (self.red | self.black) & (~self.gnd)
        nun = int(unk.sum())
        idx = np.full((ny, nx), -1, np.int64)
        idx[unk] = np.arange(nun)
        jj, ii = np.nonzero(unk)
        p = idx[jj, ii]
        self._pj, self._pi, self._pidx = jj, ii, p
        self.n_unknowns = nun
        self._nb = []
        for dj, di, cname in ((1, 0, "cN"), (-1, 0, "cS"), (0, 1, "cE"), (0, -1, "cW")):
            jn, in_ = jj + dj, ii + di                 # 내부 미지수의 이웃 -> 항상 격자 안
            q = idx[jn, in_]
            is_gnd = self.gnd[jn, in_]
            is_ndl = self.needle[jn, in_] & ~is_gnd
            is_mir = ((in_ == 0) | (in_ == nx - 1) | (jn == ny - 1)) & ~is_gnd & ~is_ndl
            is_unk = (q >= 0) & ~is_gnd & ~is_ndl & ~is_mir
            rest = ~(is_gnd | is_ndl | is_mir | is_unk)   # 금속평판 셀 / 바닥행: 0 또는 V_plate
            self._nb.append(dict(cname=cname, m_unk=is_unk, q_unk=q[is_unk], p_unk=p[is_unk],
                                 m_ndl=is_ndl, m_rest=rest, m_mir=is_mir))
        A, self._bN, self._bP = self._assemble_poisson(self.cN, self.cS, self.cE, self.cW, self.csum)
        t0 = time.perf_counter()
        self._lu = _spla.splu(A, permc_spec="MMD_AT_PLUS_A")
        self.t_wall_factor = time.perf_counter() - t0
        self.poisson_used = "direct"

    def _assemble_poisson(self, cN, cS, cE, cW, csum):
        jj, ii, p = self._pj, self._pi, self._pidx
        diag = csum[jj, ii].copy()
        bN = np.zeros(p.size); bP = np.zeros(p.size)
        rows = [p]; cols = [p]; vals = [None]
        floating = (self.plate_mode == "FLOAT_METAL")
        coef = {"cN": cN, "cS": cS, "cE": cE, "cW": cW}
        for d in self._nb:
            c = coef[d["cname"]][jj, ii]
            bN[d["m_ndl"]] += c[d["m_ndl"]]
            if floating:
                bP[d["m_rest"]] += c[d["m_rest"]]
            diag[d["m_mir"]] -= c[d["m_mir"]]                   # 미러 이웃 = 자기 자신
            rows.append(d["p_unk"]); cols.append(d["q_unk"]); vals.append(-c[d["m_unk"]])
        vals[0] = diag
        A = _sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                           shape=(p.size, p.size)).tocsc()
        return A, bN, bP

    def _semi_coefficients(self, dt, sigma):
        """반암시적 포아송: ∇·[(ε_r ε0 + dt σ)∇V] = -ρ  ->  유효 유전율 ε_eff = ε_r + dt σ/ε0
        면 계수는 v17 과 같은 조화평균."""
        e = self.epsr + (dt / EPS0) * sigma
        e[self.needle] = self.epsr[self.needle]
        e[self.gnd] = self.epsr[self.gnd]
        def harm(a_, b_): return 2.0 * a_ * b_ / (a_ + b_)
        cN = np.empty_like(e); cS = np.empty_like(e); cE = np.empty_like(e); cW = np.empty_like(e)
        cN[:-1, :] = harm(e[:-1, :], e[1:, :]);  cN[-1, :] = e[-1, :]
        cS[1:, :]  = harm(e[1:, :],  e[:-1, :]); cS[0, :]  = e[0, :]
        cE[:, :-1] = harm(e[:, :-1], e[:, 1:]);  cE[:, -1] = e[:, -1]
        cW[:, 1:]  = harm(e[:, 1:],  e[:, :-1]); cW[:, 0]  = e[:, 0]
        return cN, cS, cE, cW, cN + cS + cE + cW

    # ------------------------------------------------------------------
    def _compute_laplace_field(self):
        """단위 인가전압(Va=1 V, 공간전하 0)의 라플라스 전계 E_L1 (Sato 전류공식용).
        Sato:  I = ∫ q(n_+ mu_+ + n_e mu_e + n_- mu_-) E·E_L1 dV  [A]  (축대칭 환산)"""
        Vs, va, vp = self.V.copy(), self.Va_now, self.V_plate
        ne, nip, nin = self.n_e, self.n_ip, self.n_in
        self.n_e = np.zeros_like(Vs); self.n_ip = np.zeros_like(Vs); self.n_in = np.zeros_like(Vs)
        sig_s = self.sigma.copy(); self.sigma = np.zeros_like(self.sigma)
        self.Va_now, self.V_plate = 1.0, 0.0
        self.V[:] = 0.0
        self._apply_bc()
        self.solve_field(iters=3000, w=1.9)
        self.ExL1, self.EyL1 = self.Ex.copy(), self.Ey.copy()
        self.n_e, self.n_ip, self.n_in, self.sigma = ne, nip, nin, sig_s
        self.V[:] = Vs; self.Va_now, self.V_plate = va, vp
        self._apply_bc()

    def _sigma(self, mu_e_loc):
        """전도도 σ = q(μe ne + μ+ n+ + μ- n-) [S/m]"""
        return QE * (mu_e_loc * self.n_e + self.mu_ip * self.n_ip + self.mu_in * self.n_in)

    def discharge_current(self, mu_e_loc, sigma=None):
        """Sato 공식 방전(전도)전류 [A]: I = ∫ σ E·E_L1 dV, 2D 데카르트 -> 축대칭(|x|=r) 환산.
        (전극 정전용량의 변위전류 C·dV/dt 는 포함하지 않음)"""
        sig = self._sigma(mu_e_loc) if sigma is None else sigma
        EdotEL = self.Ex * self.ExL1 + self.Ey * self.EyL1
        return float((sig * EdotEL * self._vol_w).sum())

    def _update_applied_voltage(self, sigma=None):
        """회로: V_needle = V_src(t) - R·I.  Sato 전류를 E = Va·E_L1 + E_ρ 로 선형화하면
        I = G·Va + I_ρ (G = ∫σ E_L1·E_L1 dV 채널 컨덕턴스) 이므로 암시적으로
        Va = (V_src - R I_ρ)/(1 + R G) 를 사용 (진동·클램프 없이 안정)."""
        self.V_src_now = self.applied_voltage(self.t)
        if self.R_series > 0.0 and sigma is not None:
            w = self._vol_w
            G = float((sigma * (self.ExL1**2 + self.EyL1**2) * w).sum())
            Erx = self.Ex - self.Va_now * self.ExL1          # 공간전하 전계 (Va 무관)
            Ery = self.Ey - self.Va_now * self.EyL1
            I_rho = float((sigma * (Erx * self.ExL1 + Ery * self.EyL1) * w).sum())
            va = (self.V_src_now - self.R_series * I_rho) / (1.0 + self.R_series * G)
            if self.V_src_now >= 0:                          # 극성 반전 방지 (안전장치)
                va = min(max(va, 0.0), self.V_src_now)
            else:
                va = max(min(va, 0.0), self.V_src_now)
            self.Va_now = va
            self.G_cond = G
            self.I_dis = G * va + I_rho
        else:
            self.Va_now = self.V_src_now

    def head_z(self, thr=None):
        """축상 ne > thr 인 가장 낮은 z [m] (스트리머 선단 위치; 없으면 nan)"""
        thr = self.n_bridge_thr if thr is None else thr
        ic = self.nx // 2
        ne = self.n_e[self.j_ax_lo:self.j_ax_hi + 1, ic]
        idx = np.nonzero(ne > thr)[0]
        if idx.size == 0:
            return float("nan")
        return float(self.y[self.j_ax_lo + idx.min()] - self.y_bot)

    def _check_bridge(self):
        """축상 기체구간 전체가 ne > 문턱이면 갭 브리징(스트리머-스파크 전이 시작)으로 판정"""
        if self.bridged:
            return
        ic = self.nx // 2
        ne = self.n_e[self.j_ax_lo:self.j_ax_hi + 1, ic]
        if ne.size and ne.min() > self.n_bridge_thr:
            self.bridged = True
            self.t_bridge = self.t

    def _neutral_cap(self):
        """밀도 상한 적용 시 세 하전종을 같은 비율로 축소해 전하중성(ρ)을 보존.
        (v17: 종별 독립 clip -> 상한 도달 셀에서 인위적 순전하 -> 전계 발산·dt 붕괴)"""
        top = np.maximum(np.maximum(self.n_e, self.n_ip), self.n_in)
        over = top > self.n_cap
        if over.any():
            f = self.n_cap / top[over]
            self.n_e[over] *= f; self.n_ip[over] *= f; self.n_in[over] *= f

    # ------------------------------------------------------------------
    def _k(self, name, EN_Td):
        """단면적 기반 반응률 (v17 API 호환; 내부는 LUT 보간)"""
        i0, f = self.LUT.index(np.asarray(EN_Td, float))
        return self.LUT.get(name, i0, f)

    # ------------------------------------------------------------------
    def applied_voltage(self, t):
        if self.wave_mode == "DC":
            return self.Vpeak
        if t <= 0:
            return 0.0
        return self.Vpeak * self._wk * (np.exp(-self._wa * t) - np.exp(-self._wb * t))

    def waveform_curve(self, npts=400):
        if self.wave_mode == "DC":
            tmax = max(self.t * 1.5, 100e-9)
            tt = np.linspace(0, tmax, npts)
            return tt, np.full(npts, self.Vpeak)
        tmax = 3.0 * self.T_tail
        tt = np.linspace(0, tmax, npts)
        vv = self.Vpeak * self._wk * (np.exp(-self._wa * tt) - np.exp(-self._wb * tt))
        return tt, vv

    # ------------------------------------------------------------------
    def _zero_in_solids(self):
        for n in (self.n_e, self.n_ip, self.n_in, *self.rad.values()):
            n[self.needle] = 0.0
            n[self.gnd] = 0.0
            n[self.solid] = 0.0

    def _apply_bc(self):
        V = self.V
        if self.plate_mode == "FLOAT_METAL":
            V[0, :] = self.V_plate
            V[self.solid] = self.V_plate
        else:
            V[0, :] = 0.0
            if self.plate_mode == "GND_METAL" and self.n_solid >= 1:
                V[self.solid] = 0.0
        V[:, 0]  = V[:, 1]
        V[:, -1] = V[:, -2]
        V[-1, :] = V[-2, :]
        V[self.needle] = self.Va_now
        V[self.gnd] = 0.0

    def _rho_source(self):
        rho = QE * (self.n_ip - self.n_e - self.n_in)
        if self.plate_mode == "DIELECTRIC":
            rho[self.j_s, :] += self.sigma / self.h
        return rho

    def solve_field(self, iters=100, w=1.85, dt_semi=None, sigma=None):
        """포아송 해법. direct: 희소 LU 직접해(iters 무시), sor: v17 적-흑 SOR.
        dt_semi/sigma 가 주어지면 반암시적(Ventzek-Hagelaar) 계수로 매 스텝 재조립·재분해:
        다음 스텝(dt_semi) 동안의 전도전류에 의한 전하 이완을 전계에 미리 반영 -> 유전완화
        시간 제한(dt < ε0/σ) 없이 안정."""
        t0 = time.perf_counter()
        if self._lu is not None:
            rho = self._rho_source()
            b = (self.h**2 / EPS0) * rho[self._pj, self._pi]
            if dt_semi is not None and sigma is not None:
                dts = dt_semi * sigma
                reuse = False
                if self._semi_lu is not None:
                    lu, bN, bP, dts_old = self._semi_lu
                    # ε_eff = ε_r + dt σ/ε0 의 상대변화가 허용치 이내면 기존 분해 재사용
                    rel = np.abs(dts - dts_old) / (EPS0 * self.epsr + dts_old)
                    reuse = float(rel.max()) < self.semi_reuse_tol
                if not reuse:
                    cN, cS, cE, cW, csum = self._semi_coefficients(dt_semi, sigma)
                    A, bN, bP = self._assemble_poisson(cN, cS, cE, cW, csum)
                    lu = _spla.splu(A, permc_spec="MMD_AT_PLUS_A")
                    self._semi_lu = (lu, bN, bP, dts)
                    self.n_semi_factor += 1
                b += bN * self.Va_now
                if self.plate_mode == "FLOAT_METAL":
                    b += bP * self.V_plate
                self.V[self._pj, self._pi] = lu.solve(b)
                self.n_semi_steps += 1
                self.t_wall_semi += time.perf_counter() - t0
            else:
                b += self._bN * self.Va_now
                if self.plate_mode == "FLOAT_METAL":
                    b += self._bP * self.V_plate
                self.V[self._pj, self._pi] = self._lu.solve(b)
            self._apply_bc()
        else:
            self._solve_field_sor(iters, w)
        gy, gx = np.gradient(self.V, self.h)
        self.Ex, self.Ey = -gx, -gy
        self.Emag = np.hypot(self.Ex, self.Ey)
        self.t_wall_field += time.perf_counter() - t0

    def _solve_field_sor(self, iters=100, w=1.85):
        """v17 원본 SOR (scipy 미설치 시 폴백)"""
        src = (self.h**2) * self._rho_source() / EPS0
        V = self.V
        for _ in range(int(iters)):
            for m in (self.red, self.black):
                nbs = (self.cN * np.roll(V, -1, 0) + self.cS * np.roll(V, 1, 0) +
                       self.cE * np.roll(V, -1, 1) + self.cW * np.roll(V, 1, 1))
                V[m] = (1 - w) * V[m] + w * (nbs[m] + src[m]) / self.csum[m]
            self._apply_bc()

    def poisson_residual(self):
        """현재 V 의 이산 포아송 잔차 max|csum·V - Σc·V_nb - src| / max|src| (검증용)"""
        src = (self.h**2) * self._rho_source() / EPS0
        V = self.V
        nbs = (self.cN * np.roll(V, -1, 0) + self.cS * np.roll(V, 1, 0) +
               self.cE * np.roll(V, -1, 1) + self.cW * np.roll(V, 1, 1))
        m = (self.red | self.black) & (~self.gnd)
        res = self.csum * V - nbs - src
        scale = max(np.abs(self.csum * V)[m].max(), 1e-300)
        return float(np.abs(res[m]).max() / scale)

    def Emag_masked(self):
        E = self.Emag.copy()
        E[self.needle] = np.nan
        E[self.gnd] = np.nan
        if self.plate_mode in ("GND_METAL", "FLOAT_METAL"):
            E[self.solid] = np.nan
        return E

    # ------------------------------------------------------------------
    def _alpha_eta_townsend(self, Emag):
        Ecm = np.maximum(Emag, 1.0) / 100.0
        inv = self.p_torr / np.maximum(Ecm, 1e-2)
        alpha_cm = 0.0
        for name, f in (("AIR", self.f_air), ("AR", self.f_ar), ("H2O", self.f_h2o)):
            if f > 0:
                A, B = GAS_AB[name]
                alpha_cm = alpha_cm + f * A * self.p_torr * np.exp(-B * inv)
        eta_val = (1.0 * (self.x_o2 / 0.21) + 20.0 * self.f_h2o) * (self.p_torr / 760.0)
        return alpha_cm * 100.0, np.full_like(Ecm, eta_val * 100.0)

    def _transport(self, n, vx, vy, D, dt):
        h = self.h
        vxa = 0.5 * (vx[:, :-1] + vx[:, 1:])
        Fx  = np.where(vxa > 0, n[:, :-1], n[:, 1:]) * vxa
        vya = 0.5 * (vy[:-1, :] + vy[1:, :])
        Fy  = np.where(vya > 0, n[:-1, :], n[1:, :]) * vya
        div = np.zeros_like(n)
        div[:, 1:-1] += (Fx[:, 1:] - Fx[:, :-1]) / h
        div[1:-1, :] += (Fy[1:, :] - Fy[:-1, :]) / h
        lap = np.zeros_like(n)
        lap[1:-1, 1:-1] = (n[1:-1, 2:] + n[1:-1, :-2] +
                           n[2:, 1:-1] + n[:-2, 1:-1] - 4 * n[1:-1, 1:-1]) / h**2
        n2 = n + dt * (-div + D * lap)
        n2[self.needle] = 0.0
        n2[self.gnd] = 0.0
        return np.clip(n2, 0.0, self.n_cap)

    def _lap(self, n):
        out = np.zeros_like(n)
        out[1:-1, 1:-1] = (n[1:-1, 2:] + n[1:-1, :-2] +
                           n[2:, 1:-1] + n[:-2, 1:-1] - 4 * n[1:-1, 1:-1]) / self.h**2
        return out

    # ------------------------------------------------------------------
    def _radical_chem(self, dt):
        """상태별 활성종 화학. 여기종(O(D),N(D))·H 1차소광은 지수형 해석 갱신"""
        OP = self.rad["O(P)"]; OD = self.rad["O(D)"]
        NS = self.rad["N(S)"]; ND = self.rad["N(D)"]
        OH = self.rad["OH"];   NO = self.rad["NO"]
        O3 = self.rad["O3"]

        # ---- (1) O(1D): 배경기체 소광 -> 지수형(무조건 안정) ----
        nu1 = CHEM["kq_OD_N2"] * self.nN2      # -> O(P)
        nu2 = CHEM["kq_OD_O2"] * self.nO2      # -> O(P)
        nu3 = CHEM["k_OD_H2O"] * self.nH2O     # -> 2 OH
        nuOD = nu1 + nu2 + nu3
        if nuOD > 0:
            dOD = OD * (-np.expm1(-nuOD * dt))
            self.rad["O(D)"] = OD - dOD
            OP = OP + dOD * ((nu1 + nu2) / nuOD)
            OH = OH + dOD * (2.0 * nu3 / nuOD)

        # ---- (2) N(2D): O2 반응 + N2 소광 -> 지수형 ----
        nu4 = CHEM["k_ND_O2"]  * self.nO2      # -> NO + O(P)
        nu5 = CHEM["kq_ND_N2"] * self.nN2      # -> N(S)
        nuND = nu4 + nu5
        if nuND > 0:
            dND = ND * (-np.expm1(-nuND * dt))
            self.rad["N(D)"] = ND - dND
            NO = NO + dND * (nu4 / nuND)
            OP = OP + dND * (nu4 / nuND)
            NS = NS + dND * (nu5 / nuND)

        # ---- (3) 명시적 2차 반응 (느린 반응들) ----
        r6  = CHEM["k_NS_O2"] * NS * self.nO2 * dt        # NS+O2->NO+OP
        r7  = CHEM["k_NS_NO"] * NS * NO * dt              # NS+NO->N2+OP
        r8  = CHEM["k_NS_OH"] * NS * OH * dt              # NS+OH->NO+H
        r9  = CHEM["k3_O_O2"] * self.Ngas * OP * self.nO2 * dt   # OP+O2+M->O3
        r10 = CHEM["k_OP_OH"] * OP * OH * dt              # OP+OH->O2+H
        r11 = CHEM["k_OH_OH"] * OH * OH * dt              # OH+OH 소멸
        r12 = CHEM["k3_O_O"]  * self.Ngas * OP * OP * dt  # OP+OP+M->O2

        fNS = np.minimum(NS / np.maximum(r6 + r7 + r8, 1e-300), 1.0)
        fOH = np.minimum(OH / np.maximum(r8 + r10 + 2*r11, 1e-300), 1.0)
        fOP = np.minimum(OP / np.maximum(r9 + r10 + 2*r12, 1e-300), 1.0)
        fNO = np.minimum(NO / np.maximum(r7, 1e-300), 1.0)
        r6 *= fNS; r7 *= fNS * fNO; r8 *= fNS * fOH
        r9 *= fOP; r10 *= fOP * fOH; r11 *= fOH; r12 *= fOP

        self.rad["N(S)"] = NS - r6 - r7 - r8
        self.rad["O(P)"] = OP + r6 + r7 - r9 - r10 - 2*r12
        self.rad["OH"]   = OH - r8 - r10 - 2*r11
        self.rad["NO"]   = NO + r6 + r8 - r7
        self.rad["O3"]   = O3 + r9
        self.rad["H"]    = self.rad["H"] + r8 + r10

        # ---- (4) H: H+O2+M 소멸 -> 지수형 ----
        kH = CHEM["k3_H_O2"] * self.nO2 * self.Ngas
        self.rad["H"] = self.rad["H"] * np.exp(-kH * dt)

        for sp in RAD_SPECIES:
            np.clip(self.rad[sp], 0.0, 1e24, out=self.rad[sp])

    # ------------------------------------------------------------------
    def _record_profile(self):
        ic = self.nx // 2
        E_ax  = np.where(self.needle[:, ic] | self.gnd[:, ic],
                         np.nan, self.Emag[:, ic])
        rho_ax = QE * (self.n_ip[:, ic] - self.n_e[:, ic] - self.n_in[:, ic])
        self.profiles.append((self.t, E_ax.copy(), self.n_e[:, ic].copy(),
                              self.n_ip[:, ic].copy(), rho_ax.copy()))
        self.next_snap = self.t + self.snap_dt

    def _record_hist(self):
        """스트릭(축상 ne) 및 활성종 총개수(축대칭 환산) 이력 기록"""
        ic = self.nx // 2
        self.streak.append((self.t, self.n_e[:, ic].copy()))
        # 발광강도 ∝ ne · k_exc(E/N) · [N2]  (N2 SPS, 순시 여기율 근사)
        EN_ax = self.Emag[:, ic] / (self.Ngas * TD)
        em = self.n_e[:, ic] * self._k("nu_exc", EN_ax)
        self.streak_em.append((self.t, em))
        if len(self.streak) > 2400:                # 메모리 상한: 솎아내기
            self.streak = self.streak[::2]
            self.streak_em = self.streak_em[::2]
        # 총개수 N = ∮ n·2πr dr dz 근사 (2D 데카르트 -> 축대칭 환산, |x|=r)
        w = self._vol_w
        tot = {sp: float((self.rad[sp] * w).sum()) for sp in RAD_SPECIES}
        self.rad_hist.append((self.t, tot))
        if len(self.rad_hist) > 4000:
            self.rad_hist = self.rad_hist[::2]

    # ------------------------------------------------------------------
    def step(self, n_sub=6, field_iters=60):
        """방전(스트리머) 모드: 전계+하전입자+라디칼 생성"""
        t_w0 = time.perf_counter()
        js = self.j_s
        L = self.LUT
        for _ in range(int(n_sub)):
            Ex, Ey, Em = self.Ex, self.Ey, self.Emag
            EN = Em / (self.Ngas * TD)                     # E/N [Td]
            i0, fr = L.index(EN)                           # LUT 인덱스 1회

            # ---- 반응률/이동도: 단면적 테이블 vs Townsend 근사 ----
            if self.rate_model == "TABLE":
                mu_loc = L.get("mu", i0, fr)               # 전계의존 이동도
                vmag = mu_loc * Em
                vex, vey = -mu_loc * Ex, -mu_loc * Ey
                nu_i = L.get("nu_i", i0, fr)               # [1/s]
                nu_a = L.get("nu_a", i0, fr)
                mu_ref = max(mu_loc.max(), 1e-6)
            else:
                alpha, eta = self._alpha_eta_townsend(Em)
                mu_loc = self.mu_e
                vmag = self.mu_e * Em
                vex, vey = -self.mu_e * Ex, -self.mu_e * Ey
                nu_i = alpha * vmag
                nu_a = eta * vmag
                mu_ref = self.mu_e

            vmax = max(vmag.max(), 1e2)
            rate_max = max(nu_i.max(), 1.0)

            lims = {"CFL": 0.35 * self.h / vmax,
                    "diffusion": 0.2 * self.h**2 / self.De,
                    "ionization": 0.3 / rate_max}
            if self.dt_rule == "v17":
                tau_d = EPS0 / (QE * max(self.n_e.max(), 1e10) * mu_ref)
            else:
                # (v18) 물리적 유전완화: tau_d = eps0 / max(sigma), sigma = q Σ mu_s n_s
                sig = QE * (mu_loc * self.n_e + self.mu_ip * self.n_ip + self.mu_in * self.n_in)
                tau_d = EPS0 / max(sig.max(), QE * 1e10 * mu_ref)
            if self.wave_mode != "DC":
                lims["waveform"] = self.T_front / 40.0
            if self.dt_user > 0:
                lims["manual"] = self.dt_user
            lim_other = min(lims, key=lims.get)
            dt_other = lims[lim_other]
            dt_diel = 0.5 * tau_d
            use_semi = False
            if self.semi_implicit and self._lu is not None:
                # 반암시적 분해 비용(≈ LU 분해시간 / 명시적 스텝시간) 이상으로 dt 이득이 있을 때만 사용
                t_expl = (self.t_wall_step - self.t_wall_semi) / max(self.n_substeps - self.n_semi_steps, 1)
                ratio = 1.0 + self.t_wall_factor / max(t_expl, 1e-4) if self.n_substeps > 10 else 8.0
                use_semi = dt_diel * max(ratio, 2.0) < dt_other
            if use_semi:
                dt = dt_other
                if self._dt_prev_semi is not None:
                    dt = min(dt, 1.25 * self._dt_prev_semi)   # 반암시적 모드: dt 급증 제한
                lim_name = "semi:" + lim_other
            else:
                if dt_diel < dt_other:
                    dt, lim_name = dt_diel, "dielectric"
                else:
                    dt, lim_name = dt_other, lim_other
            self._dt_prev_semi = dt if use_semi else None
            self.dt_limited = (self.dt_user > 0) and (dt < self.dt_user * 0.999)
            self.dt_last = dt
            self.dt_limiter = lim_name

            # ---- 이류/확산 ----
            self.n_e  = self._transport(self.n_e, vex, vey, self.De, dt)
            self.n_ip = self._transport(self.n_ip,  self.mu_ip * Ex,
                                        self.mu_ip * Ey, 1e-4, dt)
            self.n_in = self._transport(self.n_in, -self.mu_in * Ex,
                                        -self.mu_in * Ey, 1e-4, dt)

            # ---- 평판 전하 수집 ----
            dep = QE * self.h * (self.n_ip[:js, :].sum(0)
                                 - self.n_e[:js, :].sum(0)
                                 - self.n_in[:js, :].sum(0))
            if self.plate_mode == "DIELECTRIC":
                self.sigma += dep
            elif self.plate_mode == "FLOAT_METAL":
                self.Q_acc += dep.mean()
                self.V_plate = float(np.clip(self.Q_acc / (EPS0 / self.gap),
                                             -1.2 * abs(self.Vpeak),
                                              1.2 * abs(self.Vpeak)))
            for n in (self.n_e, self.n_ip, self.n_in):
                n[:js, :] = 0.0

            # ---- O2- 전계의존 이탈 (지수형 무조건 안정) 및 배경 이온쌍 생성 ----
            nu_det = L.get("nu_det", i0, fr)
            d_det = self.n_in * (-np.expm1(-np.clip(nu_det * dt, 0.0, 50.0)))
            self.n_in -= d_det
            self.n_e  += d_det
            if self.q_bg > 0:
                self.n_in += self.q_bg * dt               # e는 즉시 부착 가정
                self.n_ip += self.q_bg * dt

            # ---- 전리/부착 ----
            arg_i = np.clip(nu_i * dt, 0.0, 18.0)
            arg_a = np.clip(nu_a * dt, 0.0, 18.0)
            d_ion = self.n_e * np.expm1(arg_i)
            # 밀도 상한: 전자·양이온 쌍 생성이 상한을 넘지 않도록 제한 (중성 보존)
            np.minimum(d_ion, np.maximum(self.n_cap - np.maximum(self.n_e, self.n_ip), 0.0), out=d_ion)
            d_att = self.n_e * (-np.expm1(-arg_a))
            # ---- 라디칼 생성 (전자충돌 해리, 상태별 채널) ----
            s6 = self.n_e * np.clip(L.get("nu_d6", i0, fr) * dt, 0, 5.0)     # -> 2 O(P)
            s8 = self.n_e * np.clip(L.get("nu_d8", i0, fr) * dt, 0, 5.0)     # -> O(P)+O(D)
            sN = self.n_e * np.clip(L.get("nu_dN2", i0, fr) * dt, 0, 5.0)    # -> N(S)+N(D)
            sW = self.n_e * np.clip(L.get("nu_dH2O", i0, fr) * dt, 0, 5.0)   # -> OH+H
            self.rad["O(P)"] += 2.0 * s6 + s8
            self.rad["O(D)"] += s8
            self.rad["N(S)"] += sN
            self.rad["N(D)"] += sN
            self.rad["OH"]   += sW
            self.rad["H"]    += sW

            self.n_e  += d_ion - d_att
            self.n_ip += d_ion
            self.n_in += d_att

            # ---- 음극 배경 전자방출 J0 (Townsend 개시원) ----
            if self.J0 > 0:
                s_emit = self.J0 / (QE * self.h) * dt     # [m^-3]
                if self.Va_now < 0:                       # 침이 음극
                    self.n_e[self.needle_adj] += s_emit
                elif self.Va_now > 0:                     # 하부 전극이 음극
                    self.n_e[self.bot_adj] += s_emit

            # ---- 2차전자방출 (gamma) ----
            if self.gamma_see > 0:
                if self.Va_now < 0:
                    m = self.needle_adj
                    see = (self.gamma_see * self.n_ip[m] * self.mu_ip *
                           Em[m] / self.h * dt)
                    self.n_e[m] += see
                elif self.Va_now > 0 and self.gnd.any():
                    m = self.bot_adj
                    see = (self.gamma_see * self.n_ip[m] * self.mu_ip *
                           Em[m] / self.h * dt)
                    self.n_e[m] += see
                elif self.Va_now > 0:
                    inc = self.Ey[js, :] < 0
                    see = (self.gamma_see * self.n_ip[js, inc] * self.mu_ip *
                           np.abs(self.Ey[js, inc]) / self.h * dt)
                    self.n_e[js, inc] += see
                    q_out = np.zeros(self.nx)
                    q_out[inc] = QE * see * self.h
                    if self.plate_mode == "DIELECTRIC":
                        self.sigma += q_out
                    elif self.plate_mode == "FLOAT_METAL":
                        self.Q_acc += q_out.mean()

            # ---- 재결합 / 라디칼 화학 ----
            r1 = np.minimum(self.beta * self.n_e * self.n_ip * dt,
                            np.minimum(self.n_e, self.n_ip))
            r2 = np.minimum(self.beta * self.n_in * self.n_ip * dt,
                            np.minimum(self.n_in, self.n_ip))
            self.n_e  -= r1
            self.n_ip -= (r1 + r2)
            self.n_in -= r2
            self._radical_chem(dt)

            gas = self._gas_mask
            if self.n_bg > 0:
                self.n_e[gas]  = np.maximum(self.n_e[gas],  self.n_bg)
                self.n_ip[gas] = np.maximum(self.n_ip[gas], self.n_bg)
            self._zero_in_solids()
            for n in (self.n_e, self.n_ip, self.n_in):
                np.clip(n, 0, None, out=n)
            self._neutral_cap()

            self.t += dt
            sig_new = self._sigma(mu_loc) if (use_semi or self.R_series > 0.0) else None
            self._update_applied_voltage(sig_new)
            if use_semi:
                self.solve_field(dt_semi=dt, sigma=sig_new)
            else:
                self.solve_field(iters=field_iters)
            self.n_substeps += 1
            if self.dt_floor > 0.0:
                self._n_stall = self._n_stall + 1 if dt < self.dt_floor else 0
                if self._n_stall >= 20:
                    self.stalled = True

            if self.t >= self.next_snap:
                self._record_profile()
        self._townsend_metrics()
        self._initiation_stats()
        self._record_hist()
        self._check_bridge()
        if self.R_series <= 0.0:                     # 진단용 전류 (프레임마다)
            mu_d = (self.LUT.get("mu", *self.LUT.index(self.Emag / (self.Ngas * TD)))
                    if self.rate_model == "TABLE" else self.mu_e)
            self.I_dis = self.discharge_current(mu_d)
        self.circ_hist.append((self.t, self.V_src_now, self.Va_now, self.I_dis))
        if len(self.circ_hist) > 4000:
            self.circ_hist = self.circ_hist[::2]
        self.dt_hist.append((self.t, self.dt_last, self.dt_limiter))
        if len(self.dt_hist) > 20000:
            self.dt_hist = self.dt_hist[::2]
        self.t_wall_step += time.perf_counter() - t_w0

    # ------------------------------------------------------------------
    def run_until(self, t_end, n_sub=6, on_frame=None, wall_limit=None,
                  stop_on_bridge=False, post_bridge=0.0, stop_on_stall=True):
        """헤드리스 실행: t >= t_end 까지 step(n_sub) 반복.
        on_frame(sim) 을 프레임마다 호출(로그/스냅샷용), wall_limit[s] 초과 시 중단.
        stop_on_bridge: 갭 브리징 후 post_bridge[s] 만큼 더 진행하고 중단
        stop_on_stall : dt 가 dt_floor 아래로 붕괴(20 서브스텝 연속)하면 중단
        반환: dict(frames, substeps, wall_s, t_end_reached, stop_reason)"""
        t0 = time.perf_counter()
        frames = 0
        reason = "t_end"
        while self.t < t_end:
            self.step(n_sub=n_sub)
            frames += 1
            if on_frame is not None:
                on_frame(self)
            if wall_limit is not None and (time.perf_counter() - t0) > wall_limit:
                reason = "wall_limit"; break
            if stop_on_bridge and self.bridged and self.t >= self.t_bridge + post_bridge:
                reason = "bridged"; break
            if stop_on_stall and self.stalled:
                reason = "dt_stalled"; break
        return dict(frames=frames, substeps=self.n_substeps,
                    wall_s=time.perf_counter() - t0, t_end_reached=self.t >= t_end,
                    stop_reason=reason)

    # ------------------------------------------------------------------
    def _townsend_metrics(self):
        """축상 alpha_eff 적분으로 증배계수 M=e^{int alpha dz}, gamma(M-1) 산출"""
        ic = self.nx // 2
        Eax = self.Emag[:, ic]
        ENax = Eax / (self.Ngas * TD)
        if self.rate_model == "TABLE":
            i0, fr = self.LUT.index(ENax)
            v = np.maximum(self.LUT.get("mu", i0, fr) * Eax, 1e-3)
            nu = self.LUT.get("nu_i", i0, fr) - self.LUT.get("nu_a", i0, fr)
            a_eff = nu / v
        else:
            al, et = self._alpha_eta_townsend(Eax)
            a_eff = al - et
        a_eff = np.where(self.needle[:, ic] | self.gnd[:, ic],
                         0.0, np.maximum(a_eff, 0.0))
        j_lo = self.j_ax_lo
        j_hi = max(self.j_ax_hi, j_lo)
        integ = float(a_eff[j_lo:j_hi + 1].sum() * self.h)
        self.M_townsend = float(np.exp(min(integ, 200.0)))
        self.gamma_MK = self.gamma_see * (self.M_townsend - 1.0)

    # ------------------------------------------------------------------
    def _initiation_stats(self):
        """유효 개시율 Lambda(t)=∫_{a_eff>0} nu_det·n(O2-) dV, H, F 기록"""
        EN = self.Emag / (self.Ngas * TD)
        i0, fr = self.LUT.index(EN)
        if self.rate_model == "TABLE":
            v = np.maximum(self.LUT.get("mu", i0, fr) * self.Emag, 1e-3)
            nu_net = self.LUT.get("nu_i", i0, fr) - self.LUT.get("nu_a", i0, fr)
            a_eff = nu_net / v
        else:
            al, et = self._alpha_eta_townsend(self.Emag)
            a_eff = al - et
        trig = (a_eff > 0) & (~self.needle) & (~self.gnd) & (~self.solid)
        nu_det = self.LUT.get("nu_det", i0, fr)
        Lam = float((nu_det * self.n_in * self._vol_w)[trig].sum())     # [1/s]
        dtf = max(self.t - self._t_prev_frame, 0.0)
        self._t_prev_frame = self.t
        self._H_hazard += Lam * dtf
        F = 1.0 - np.exp(-min(self._H_hazard, 500.0))
        self.init_stats.append((self.t, Lam, self._H_hazard, F))
        if len(self.init_stats) > 4000:
            self.init_stats = self.init_stats[::2]

    def sample_delays(self, n=2000, rng=None):
        """누적위험 H(t) 역변환으로 통계적 지연시간 표본 추출.
        반환: (유효표본[s], 검열비율: t_end까지 미개시 확률)"""
        if len(self.init_stats) < 3:
            return np.array([]), 1.0
        tt = np.array([p[0] for p in self.init_stats])
        HH = np.array([p[2] for p in self.init_stats])
        if HH[-1] <= 0:
            return np.array([]), 1.0
        rng = rng or np.random.default_rng()
        target = -np.log(1.0 - rng.random(int(n)))
        ok = target <= HH[-1]
        samples = np.interp(target[ok], HH, tt)
        return samples, float(1.0 - ok.mean())

    # ------------------------------------------------------------------
    def step_afterglow(self, n_sub=4):
        """잔광 모드: 전계/하전수송 정지, 라디칼 확산+화학 (us 오더 dt)"""
        for _ in range(int(n_sub)):
            # 화학 안정 dt: 명시적(느린) 반응만 기준.
            # 여기종 O(D)/N(D)·H 소광은 지수형 해석 갱신이므로 dt를 제한하지 않음
            rmax = CHEM["k3_O_O2"] * self.nO2 * self.Ngas
            rmax = max(rmax, CHEM["k_NS_NO"] * self.rad["NO"].max())
            rmax = max(rmax, CHEM["k_NS_OH"] * self.rad["OH"].max())
            rmax = max(rmax, CHEM["k_OP_OH"] *
                       max(self.rad["O(P)"].max(), self.rad["OH"].max()))
            rmax = max(rmax, CHEM["k_OH_OH"] * self.rad["OH"].max())
            rmax = max(rmax, CHEM["k3_O_O"] * self.Ngas * self.rad["O(P)"].max())
            dt = min(0.2 / max(rmax, 1.0), 0.25 * self.h**2 / self.D_rad, 2e-6)
            self.dt_last = dt

            # 하전입자 감쇠 (무조건 안정: 지수/해석형)
            nu_a0 = (self._k("k_att_O2", np.array([1.0]))[0] * self.nO2 +
                     self._k("k_att_H2O", np.array([1.0]))[0] * self.nH2O)
            fac = np.exp(-nu_a0 * dt)
            d = self.n_e * (1.0 - fac)
            self.n_e *= fac
            self.n_in += d
            ntot = self.n_e + self.n_in
            rem = (self.beta * self.n_ip * ntot * dt /
                   (1.0 + self.beta * dt * np.maximum(self.n_ip, ntot)))
            rem = np.minimum(rem, np.minimum(self.n_ip, ntot))
            fr = rem / np.maximum(ntot, 1e-30)
            self.n_ip -= rem
            self.n_e  -= self.n_e * fr
            self.n_in -= self.n_in * fr

            # 라디칼 확산 + 화학
            for sp in RAD_SPECIES:
                self.rad[sp] = np.clip(
                    self.rad[sp] + dt * self.D_rad * self._lap(self.rad[sp]),
                    0.0, 1e24)
            self._radical_chem(dt)
            self._zero_in_solids()
            self.t += dt
        self._record_hist()
