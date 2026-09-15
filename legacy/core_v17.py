# -*- coding: utf-8 -*-
"""legacy/core_v17.py
v17 원본(needle_plane_discharge_hmi_v17_original.py)의 물리 코어를 tkinter/matplotlib
의존 없이 import 할 수 있도록 EPS0 ~ NeedlePlaneDischarge 클래스 끝까지 '원문 그대로'
잘라낸 파일입니다. 벤치마크 기준(reference)용이며 수정하지 마십시오.
"""
import os
import glob
import time
import numpy as np

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


def get_hv_dir():
    """바탕화면\\HighVoltage 폴더 탐색 (OneDrive/한글 바탕화면 포함), 없으면 생성"""
    home = os.path.expanduser("~")
    cands = [os.path.join(home, "Desktop", "HighVoltage"),
             os.path.join(home, "OneDrive", "Desktop", "HighVoltage"),
             os.path.join(home, "OneDrive", "바탕 화면", "HighVoltage"),
             os.path.join(home, "바탕화면", "HighVoltage"),
             os.path.join(home, "바탕 화면", "HighVoltage")]
    for c in cands:
        if os.path.isdir(c):
            return c
    for c in cands[:2]:                       # Desktop 후보에 생성 시도
        if os.path.isdir(os.path.dirname(c)):
            try:
                os.makedirs(c, exist_ok=True)
                return c
            except Exception:
                pass
    return home                               # 최후 fallback: 홈 디렉토리


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
                 rog_halfw_mm=6.0, rog_edge_mm=2.0):

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
            for k, v in rate_override.items():
                self.RT[k] = np.array(v, float)
        self._logEN = np.log10(self.RT["EN"])

        self.gamma_see = max(gamma_see, 0.0)
        self.dt_user   = max(dt_user_ns, 0.0) * 1e-9
        self.snap_dt   = max(snap_ns, 0.1) * 1e-9
        self.next_snap = 0.0
        self.profiles  = []
        self.streak    = []      # [(t, 축상 ne 프로파일)]  z-t 스트릭용
        self.streak_em = []      # [(t, 축상 발광강도)]      발광 스트릭용
        self.rad_hist  = []      # [(t, {종: 총개수})]      활성종 이력용

        self.afterglow = False                    # 잔광 모드 플래그

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

        self.t = 0.0
        self.dt_last = 0.0
        self.dt_limited = False
        self.Va_now = self.applied_voltage(0.0)
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
    def _k(self, name, EN_Td):
        """단면적 기반 반응률 log-log 보간 (범위 밖은 경계값 고정)"""
        tab = self.RT[name]
        lE = np.log10(np.clip(EN_Td, self.RT["EN"][0], self.RT["EN"][-1]))
        return 10.0 ** np.interp(lE, self._logEN, np.log10(np.maximum(tab, 1e-30)))

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

    def solve_field(self, iters=100, w=1.85):
        rho = QE * (self.n_ip - self.n_e - self.n_in)
        if self.plate_mode == "DIELECTRIC":
            rho = rho.copy()
            rho[self.j_s, :] += self.sigma / self.h
        src = (self.h**2) * rho / EPS0
        V = self.V
        for _ in range(int(iters)):
            for m in (self.red, self.black):
                nbs = (self.cN * np.roll(V, -1, 0) + self.cS * np.roll(V, 1, 0) +
                       self.cE * np.roll(V, -1, 1) + self.cW * np.roll(V, 1, 1))
                V[m] = (1 - w) * V[m] + w * (nbs[m] + src[m]) / self.csum[m]
            self._apply_bc()
        gy, gx = np.gradient(V, self.h)
        self.Ex, self.Ey = -gx, -gy
        self.Emag = np.hypot(self.Ex, self.Ey)

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
        return np.clip(n2, 0.0, 1e22)

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
        em = self.n_e[:, ic] * self._k("k_exc_N2", EN_ax) * self.nN2
        self.streak_em.append((self.t, em))
        if len(self.streak) > 2400:                # 메모리 상한: 솎아내기
            self.streak = self.streak[::2]
            self.streak_em = self.streak_em[::2]
        # 총개수 N = ∮ n·2πr dr dz 근사 (2D 데카르트 -> 축대칭 환산, |x|=r)
        w = np.pi * np.abs(self.X) * self.h * self.h
        tot = {sp: float((self.rad[sp] * w).sum()) for sp in RAD_SPECIES}
        self.rad_hist.append((self.t, tot))
        if len(self.rad_hist) > 4000:
            self.rad_hist = self.rad_hist[::2]

    # ------------------------------------------------------------------
    def step(self, n_sub=6, field_iters=60):
        """방전(스트리머) 모드: 전계+하전입자+라디칼 생성"""
        js = self.j_s
        for _ in range(int(n_sub)):
            Ex, Ey, Em = self.Ex, self.Ey, self.Emag
            EN = Em / (self.Ngas * TD)                     # E/N [Td]

            # ---- 반응률/이동도: 단면적 테이블 vs Townsend 근사 ----
            if self.rate_model == "TABLE":
                mu_loc = self._k("muN", EN) / self.Ngas    # 전계의존 이동도
                vmag = mu_loc * Em
                vex, vey = -mu_loc * Ex, -mu_loc * Ey
                nu_i = (self._k("k_ion_N2", EN) * self.nN2 +
                        self._k("k_ion_O2", EN) * self.nO2 +
                        self._k("k_ion_Ar", EN) * self.nAr +
                        self._k("k_ion_H2O", EN) * self.nH2O)   # [1/s]
                nu_a = (self._k("k_att_O2", EN) * self.nO2 +
                        self._k("k_att_H2O", EN) * self.nH2O)
                mu_ref = max(mu_loc.max(), 1e-6)
            else:
                alpha, eta = self._alpha_eta_townsend(Em)
                vmag = self.mu_e * Em
                vex, vey = -self.mu_e * Ex, -self.mu_e * Ey
                nu_i = alpha * vmag
                nu_a = eta * vmag
                mu_ref = self.mu_e

            vmax = max(vmag.max(), 1e2)
            rate_max = max(nu_i.max(), 1.0)

            dt_stab = 0.35 * self.h / vmax
            dt_stab = min(dt_stab, 0.2 * self.h**2 / self.De)
            dt_stab = min(dt_stab, 0.3 / rate_max)
            tau_d = EPS0 / (QE * max(self.n_e.max(), 1e10) * mu_ref)
            dt_stab = min(dt_stab, 0.5 * tau_d)
            if self.wave_mode != "DC":
                dt_stab = min(dt_stab, self.T_front / 40.0)
            if self.dt_user > 0:
                dt = min(self.dt_user, dt_stab)
                self.dt_limited = dt < self.dt_user * 0.999
            else:
                dt = dt_stab
                self.dt_limited = False
            self.dt_last = dt

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
            nu_det = self._k("k_det_O2m", EN) * self.Ngas
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
            d_att = self.n_e * (-np.expm1(-arg_a))
            # ---- 라디칼 생성 (전자충돌 해리, 상태별 채널) ----
            nu_a6  = self._k("k_diss_O2_6eV", EN) * self.nO2   # -> 2 O(P)
            nu_b8  = self._k("k_diss_O2_8eV", EN) * self.nO2   # -> O(P)+O(D)
            nu_dN2 = self._k("k_diss_N2", EN)  * self.nN2      # -> N(S)+N(D)
            nu_dH2O = self._k("k_diss_H2O", EN) * self.nH2O    # -> OH+H
            s6 = self.n_e * np.clip(nu_a6 * dt, 0, 5.0)
            s8 = self.n_e * np.clip(nu_b8 * dt, 0, 5.0)
            sN = self.n_e * np.clip(nu_dN2 * dt, 0, 5.0)
            sW = self.n_e * np.clip(nu_dH2O * dt, 0, 5.0)
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

            gas = ~self.solid
            self.n_e[gas]  = np.maximum(self.n_e[gas],  self.n_bg)
            self.n_ip[gas] = np.maximum(self.n_ip[gas], self.n_bg)
            self._zero_in_solids()
            for n in (self.n_e, self.n_ip, self.n_in):
                np.clip(n, 0, 1e22, out=n)

            self.t += dt
            self.Va_now = self.applied_voltage(self.t)
            self.solve_field(iters=field_iters)

            if self.t >= self.next_snap:
                self._record_profile()
        self._townsend_metrics()
        self._initiation_stats()
        self._record_hist()

    # ------------------------------------------------------------------
    def _townsend_metrics(self):
        """축상 alpha_eff 적분으로 증배계수 M=e^{int alpha dz}, gamma(M-1) 산출"""
        ic = self.nx // 2
        Eax = self.Emag[:, ic]
        ENax = Eax / (self.Ngas * TD)
        if self.rate_model == "TABLE":
            v = np.maximum(self._k("muN", ENax) / self.Ngas * Eax, 1e-3)
            nu = (self._k("k_ion_N2", ENax) * self.nN2 +
                  self._k("k_ion_O2", ENax) * self.nO2 +
                  self._k("k_ion_Ar", ENax) * self.nAr +
                  self._k("k_ion_H2O", ENax) * self.nH2O -
                  self._k("k_att_O2", ENax) * self.nO2 -
                  self._k("k_att_H2O", ENax) * self.nH2O)
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
        if self.rate_model == "TABLE":
            v = np.maximum(self._k("muN", EN) / self.Ngas * self.Emag, 1e-3)
            nu_net = (self._k("k_ion_N2", EN) * self.nN2 +
                      self._k("k_ion_O2", EN) * self.nO2 +
                      self._k("k_ion_Ar", EN) * self.nAr +
                      self._k("k_ion_H2O", EN) * self.nH2O -
                      self._k("k_att_O2", EN) * self.nO2 -
                      self._k("k_att_H2O", EN) * self.nH2O)
            a_eff = nu_net / v
        else:
            al, et = self._alpha_eta_townsend(self.Emag)
            a_eff = al - et
        trig = (a_eff > 0) & (~self.needle) & (~self.gnd) & (~self.solid)
        nu_det = self._k("k_det_O2m", EN) * self.Ngas
        w = np.pi * np.abs(self.X) * self.h * self.h          # 축대칭 부피가중
        Lam = float((nu_det * self.n_in * w)[trig].sum())     # [1/s]
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
