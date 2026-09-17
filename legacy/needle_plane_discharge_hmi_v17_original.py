# -*- coding: utf-8 -*-
"""
==============================================================================
 고전압공학 방전현상 시뮬레이터 v17 (다중 전극 구성)
==============================================================================
 [v17 변경사항]
  1) Plane-Plane 구성을 Rogowski-Rogowski 쌍으로 대체
     - 상·하 전극 모두 평탄부 + 가우시안 램프 연부(특성폭 2T, 연장 4T)의
       로고우스키형 단면 (Needle-Rogowski의 하부 전극도 동일 프로파일)
     - 평탄부 사이 균일전계, 연부 매끈전계 기준 과전계 ~10% 수준으로 억제
       (90° 원호 종단은 접지면 위 돌출 능선 효과로 ~30-40% 과전계 -> 개선)
  2) 기체 조성에 N2 분율 직접 입력 추가 (N2/O2/Ar/H2O 4성분)
     - 네 입력의 합으로 자동 정규화하여 몰분율 산출 (합이 100이 아니어도 됨)
  3) 압력 입력(Gas pressure [atm])은 기존과 동일하게 반영되며(기체밀도 N,
     환산이동도, Townsend 근사 p 스케일링), 상태창에 Pressure 표기 추가
 [v16 변경사항]
  0) 전극 구성(Configuration) 선택 기능
     - Needle-Plane(기존) / Needle-Needle / 3 Needles-Plane /
       Needle-Rogowski plane / Plane-Plane / Sphere-Sphere
     - 하부 형상 접지전극 마스크(self.gnd) 신설: Dirichlet V=0, 하전입자 흡수,
       음극 인접셀(bot_adj)에서 J0 방출·gamma 2차방출 동작
     - 대칭축 기체구간(j_ax_lo~j_ax_hi)과 z=0 기준면(y_bot)을 구성 독립적으로
       일반화 -> ⑦탭/Townsend M/개시통계가 모든 구성에서 동작
     - 침-침/구-구/로고우스키 구성은 하부가 형상 금속이므로 평판 재료 설정은
       자동으로 금속(접지)로 강제됨
 [v15 변경사항]
  1) GUI 표기를 전부 영문으로 변경 (위젯 글꼴: Times New Roman 12 pt)
     - 그래프 텍스트는 v14와 동일하게 Times New Roman 16 pt 유지
  2) 상태창을 영문으로 전면 재작성하면서, 이전 편집에서 누락되었던
     진단 항목(SEE γ, J0, M=e^∫αdz, γ(M-1), Λ(t), F(t), 스냅샷 수) 복원
 [v14 변경사항]
  1) 기체 조성: O2/Ar/H2O 분율 입력 -> N2 = 100 - (O2+Ar+H2O) 자동 산출·표시
     (기본값: 대기 조성 O2 21%, Ar 1%, H2O 0% -> N2 78%)
  2) 기타 재료 εr 기본값 0 (선택 시 실제 재료값 입력 요망; 내부에서 >=1로 보정)
  3) 순전하밀도 축라벨을 그리스 문자 ρ(mathtext)로 표기 (눈금은 ASCII 유지)
  4) 흑백 인쇄 대응: 다곡선 그래프에 색상 + 선종(실선/파선/일점쇄선/점선 등)
     + 마커(o, s, ^, v, D, ...)를 함께 부여하여 무채색에서도 곡선 구분 가능
  5) 전체 글꼴: Times New Roman, 최소 16 pt (GUI 위젯 및 모든 그래프 텍스트)
 [v13 변경사항]
  0) O2- 소이온 저장고 + 전계의존 이탈(detachment) -> 개시통계 모델
     - 대기 배경 소이온 n(O2-) = n(+) ~ 1e9 m^-3 (조정가능) 를 초기·상시 분포로 설정
       (배경 생성률 q [쌍/m^3/s, 우주선·자연방사선] 로 상시 보충)
     - 전계의존 충돌이탈 O2- + M -> e + O2 + M : k_det(E/N) 테이블 (log-log 보간,
       ~100 Td에서 이탈수명 ~10 ns 오더, Pancheshnyi PSST 2013 경향)
       -> 전계 인가 순간 고전계역의 O2- 에서 전자가 '풀려나' 방전을 개시
       -> 스트리머 채널 내 부착-이탈 재순환도 자동 반영
     - 개시 확률통계 (통계적 지연시간):
         유효 개시율 Lambda(t) = ∫_{alpha_eff>0} nu_det·n(O2-) dV  [1/s, 축대칭 환산]
         누적위험 H(t) = ∫Lambda dt,  개시확률 F(t) = 1 - e^{-H}
         역변환 표본추출로 지연시간 분포 히스토그램 생성 (N=2000)
     - 신규 "⑧ 개시통계" 탭: Lambda(t)·F(t) 곡선 + 지연시간 히스토그램
     - CSV 저장에 initiation_stats 포함
 [v12 변경사항]
  0) 타운젠트(Townsend) 개시 기구의 명시적 반영
     - 배경 전자밀도/seed 밀도를 0으로 설정 가능 (기본값 0)
     - 음극 표면 미소 전자방출 전류밀도 J0 [A/m^2] 입력 (우주선/광전자 등
       배경 방출에 해당; 침이 음극이면 침 표면, 양극이면 평판 표면에서 방출)
       -> 방출 전자가 양극으로 가속되며 alpha 증배(e^{int alpha dz}) 형성,
          양이온이 음극으로 되돌아와 gamma 2차전자방출로 후속 세대 생성
     - 자립방전 판정지표 실시간 표시: M = e^{int alpha_eff dz} (축상 적분),
       gamma*(M-1) >= 1 이면 Townsend 자립조건 충족
     - 주의: gamma 되먹임의 세대 주기는 이온 주행시간(수 us@cm급 갭)이므로
       ns 펄스에서는 침 음극(짧은 이온 귀환거리) 또는 DC 조건에서 관찰 권장
 [v11 변경사항]
  1) 그래프 높이 표기를 전 탭에서 y -> z 로 통일 (z=0: 평판, z=gap: 침 선단)
  2) 그래프 내 모든 캡션(제목/축라벨/컬러바/범례)을 영어로 통일
     -> 한글 폰트 미탑재 환경에서의 글자 깨짐(□) 원천 차단
  3) 축상 순전하밀도(symlog) 음수 눈금: ASCII 포맷터 적용 (유니코드 마이너스
     U+2212 가 폰트에 없어 깨지던 문제 수정)
  4) 축상 기록간격 기본값 5 ns -> 1 ns
 [v10 변경사항]
  0) ⑦탭 z축을 실제 전극배치 물리좌표로 정합 (아래쪽이 0)
     - z = (격자 y) - (평판 표면 y) : 평판 표면이 z=0(하단),
       위로 증가하여 침 선단이 z=gap [mm](상단)
     - 축 범위를 정확히 z = 0 ~ gap 으로 고정
       (기존: 격자 절단으로 표시 갭이 실제보다 1셀만큼 작게 나타남)
 [v9 변경사항]
  0) 발광(Optical Emission) 스트릭 이미지 추가
     - 발광강도 ∝ 전자충돌 N2(C3Πu) 여기율: S_em = ne · k_exc_N2(E/N) · [N2]
       (대기압에서 N2(C) 유효수명 ~1 ns 이하(소광지배) -> 발광 ≈ 순시 여기율)
     - "⑦ 스트릭/이력" 탭 (c)패널: z-t 발광 스트릭 (스트릭 카메라 측정 대응)
     - CSV 저장에 streak_emission 행렬 포함
 [v8 변경사항]
  0) 여기상태 분리: O(3P)/O(1D), N(4S)/N(2D)  -> 총 8종 활성종
     [상태별 해리채널]
       e+O2 --(6.0eV)--> O(P)+O(P)      : k_diss_O2_6eV
       e+O2 --(8.4eV)--> O(P)+O(D)      : k_diss_O2_8eV
       e+N2 --(~13eV)--> N(S)+N(D)      : k_diss_N2
       e+H2O           -> OH + H        : k_diss_H2O
     [소광(quenching) 및 상태별 반응]
       O(D)+N2 -> O(P)+N2 (3.1e-17)     O(D)+O2 -> O(P)+O2 (4.0e-17)
       O(D)+H2O -> 2OH     (2.2e-16)  * 잔광 OH의 주요 생성경로
       N(D)+O2 -> NO+O(P)  (6.0e-18)  * NO의 주요 생성경로
       N(D)+N2 -> N(S)+N2  (1.7e-20)
       N(S)+O2 -> NO+O(P)  (9e-23, 상온 Zeldovich: 느림)
       N(S)+NO -> N2+O(P)  (3.0e-17)    N(S)+OH -> NO+H (4.7e-17)
     - 여기종(O(D),N(D))·H 의 1차 소광은 지수형 해석 갱신(무조건 안정)
       -> ns급 소광수명이 잔광 dt(us)를 제한하지 않음
     - 반응률 CSV: k_diss_O2_6eV/_8eV 열 인식. 구버전 k_diss_O2 열은
       6eV:8.4eV = 40:60 으로 자동 분배(하위호환)
 [v7 변경사항]
  0) 논문형 그래프 탭 "⑦ 스트릭/이력" 추가 (첨부 그림 재현)
     (a) E/N [Td] 2D 맵 (z축: 선단=0 -> 평판, r축 [mm])
     (b) 전자밀도 [cm^-3] 2D 맵 (10^9~10^16, log)
     (c) 스트릭 이미지: z-t 다이어그램 (축상 log ne, 스트리머 진전 궤적)
     (d) 방전 후 활성종 총개수 N(t) [개]: O,N,OH,NO,O3,H (log-log, 축대칭 환산)
     - 라디칼 H 추가 (e+H2O->OH+H), H+O2+M->HO2 소멸
     - 시간이력/스트릭은 매 프레임 자동 기록, CSV 저장에 포함
 [v6 변경사항]
  1) HMI 레이아웃 개선: 시작/정지/1프레임/잔광/CSV저장/반응률로드 버튼을
     좌측 패널 "최상단 고정 제어바"로 이동 -> 화면 크기와 무관하게 항상 표시.
     파라미터 입력 영역은 마우스휠 스크롤 지원.
  2) BOLSIG+ 반응률 CSV 자동 연동:
     - 바탕화면/HighVoltage 폴더를 자동 탐색(없으면 생성)
     - 프로그램 시작 시 폴더 내 최신 CSV(헤더 EN_Td,...)를 자동 로드
     - "📂 반응률 로드" 버튼의 기본 폴더 = HighVoltage
     - "💾 CSV 저장" 기본 폴더 = HighVoltage
     CSV 형식: 1행 헤더 EN_Td,k_ion_N2,k_ion_O2,... (부분 교체 허용, SI: m3/s)
  * v5 기능(라디칼 5종 O/N/OH/NO/O3, 단면적 기반 반응률 테이블, 잔광모드,
    부유/유전체 평판, 순전하밀도, SEE, 펄스파형 등) 모두 유지

 [실행 - VS Code]  pip install numpy matplotlib  ->  python 이파일.py
==============================================================================
"""

import os
import glob
import time
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib import cm
from matplotlib.ticker import FuncFormatter
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = ["Times New Roman", "Nimbus Roman", "Liberation Serif",
                               "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"          # Times 계열 수식 글꼴
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams.update({"font.size": 16, "axes.titlesize": 17, "axes.labelsize": 16,
                     "xtick.labelsize": 16, "ytick.labelsize": 16,
                     "legend.fontsize": 16})

FONT_UI   = ("Times New Roman", 12)                # GUI 위젯 글꼴
FONT_UI_B = ("Times New Roman", 12, "bold")

# 흑백 인쇄 대응 선종/마커 순환
LSTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1)),
           (0, (1, 1)), (0, (5, 2, 1, 2))]
MARKERS = ["o", "s", "^", "v", "D", "P", "X", "*", "<", ">"]

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


# =============================================================================
#  HMI (tkinter GUI)
# =============================================================================
class DischargeHMI:
    MAX_OVERLAY = 10

    def __init__(self, root):
        self.root = root
        root.title("HV Discharge Simulator v17 | Rogowski Pair, 4-Gas Mixture")
        root.geometry("1460x950")
        self.sim = None
        self.running = False
        self.frame_count = 0
        self.rate_override = None
        self.hv_dir = get_hv_dir()            # 바탕화면\HighVoltage
        self._build_panel()
        self._build_tabs()
        self._auto_load_rates()               # 폴더 내 최신 BOLSIG+ CSV 자동 로드
        self.apply_params()

    # ------------------------------------------------------------------
    def _auto_load_rates(self):
        """HighVoltage 폴더에서 최신 반응률 CSV(헤더 EN_Td) 자동 로드"""
        try:
            files = sorted(glob.glob(os.path.join(self.hv_dir, "*.csv")),
                           key=os.path.getmtime, reverse=True)
            for p in files:
                try:
                    ov, loaded = parse_rate_csv(p)
                except Exception:
                    continue                   # 반응률 형식이 아닌 CSV는 건너뜀
                self.rate_override = ov
                self.lbl_rate.config(
                    text=f"Rates: {os.path.basename(p)} ({len(loaded)} auto-loaded)")
                return
            self.lbl_rate.config(text="Rates: built-in tables (no CSV in folder)")
        except Exception:
            self.lbl_rate.config(text="Rates: built-in tables")

    # ------------------------------------------------------------------
    def _build_panel(self):
        panel = ttk.Frame(self.root, padding=4)
        panel.pack(side=tk.LEFT, fill=tk.Y)
        panel.pack_propagate(False)
        panel.configure(width=360)

        # ================= (A) 최상단 고정 제어바: 항상 표시 =================
        ctrl = ttk.LabelFrame(panel, text="Controls (always visible)", padding=4)
        ctrl.pack(side=tk.TOP, fill=tk.X)
        r1 = ttk.Frame(ctrl); r1.pack(fill=tk.X, pady=1)
        ttk.Button(r1, text="⟳ Apply/Reset",
                   command=self.apply_params).pack(side=tk.LEFT, expand=True, fill=tk.X)
        self.btn_run = ttk.Button(r1, text="▶ Run", command=self.toggle_run)
        self.btn_run.pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(r1, text="⏭ 1 frame",
                   command=self.step_once).pack(side=tk.LEFT, expand=True, fill=tk.X)
        r2 = ttk.Frame(ctrl); r2.pack(fill=tk.X, pady=1)
        self.btn_ag = ttk.Button(r2, text="🌙 Afterglow", command=self.toggle_afterglow)
        self.btn_ag.pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(r2, text="💾 Save CSV",
                   command=self.export_csv).pack(side=tk.LEFT, expand=True, fill=tk.X)
        ttk.Button(r2, text="📂 Load rates",
                   command=self.load_rate_csv).pack(side=tk.LEFT, expand=True, fill=tk.X)
        self.lbl_rate = ttk.Label(ctrl, text="Rates: built-in tables",
                                  foreground="#0055aa", font=("Times New Roman", 12))
        self.lbl_rate.pack(anchor="w")
        ttk.Label(ctrl, text=f"Work folder: {self.hv_dir}", foreground="#888",
                  font=("Times New Roman", 12), wraplength=330).pack(anchor="w")

        # ================= (B) 스크롤 가능한 파라미터/상태 영역 =================
        body = ttk.Frame(panel)
        body.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        cv = tk.Canvas(body, highlightthickness=0, width=335)
        sb = ttk.Scrollbar(body, orient="vertical", command=cv.yview)
        inner = ttk.Frame(cv)
        inner.bind("<Configure>",
                   lambda e: cv.configure(scrollregion=cv.bbox("all")))
        cv.create_window((0, 0), window=inner, anchor="nw")
        cv.configure(yscrollcommand=sb.set)
        cv.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb.pack(side=tk.RIGHT, fill=tk.Y)

        def _wheel(e):                        # 마우스휠 스크롤 (Win/Linux)
            if getattr(e, "num", None) == 4 or getattr(e, "delta", 0) > 0:
                cv.yview_scroll(-2, "units")
            else:
                cv.yview_scroll(2, "units")
        def _bind_wheel(_):
            for ev in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                cv.bind_all(ev, _wheel)
        def _unbind_wheel(_):
            for ev in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                cv.unbind_all(ev)
        inner.bind("<Enter>", _bind_wheel)
        inner.bind("<Leave>", _unbind_wheel)

        self.vars = {}
        def add_entry(parent, label, key, default, unit=""):
            f = ttk.Frame(parent); f.pack(fill=tk.X, pady=1)
            ttk.Label(f, text=label, width=17, font=FONT_UI).pack(side=tk.LEFT)
            v = tk.StringVar(value=str(default))
            ttk.Entry(f, textvariable=v, width=7, font=FONT_UI).pack(side=tk.LEFT)
            ttk.Label(f, text=unit, font=FONT_UI).pack(side=tk.LEFT)
            self.vars[key] = v

        ttk.Label(inner, text="Electrodes / Gas", font=FONT_UI_B).pack(anchor="w")
        cfF = ttk.Frame(inner); cfF.pack(fill=tk.X, pady=1)
        ttk.Label(cfF, text="Configuration", width=16,
                  font=FONT_UI).pack(side=tk.LEFT)
        self.cfg = tk.StringVar(value="Needle - Plane")
        ttk.Combobox(cfF, textvariable=self.cfg, width=18, state="readonly",
                     values=list(CONFIG_MAP), font=FONT_UI).pack(side=tk.LEFT)
        add_entry(inner, "Gap (electrode gap)", "gap", 10.0, "mm")
        add_entry(inner, "Needle shank radius", "rmain", 1.0, "mm")
        add_entry(inner, "Needle tip radius", "rtip", 0.15, "mm")
        add_entry(inner, "Taper length", "taper", 3.0, "mm")
        add_entry(inner, "Needle spacing d (3-needle)", "nsep", 4.0, "mm")
        add_entry(inner, "Sphere radius R", "srad", 5.0, "mm")
        add_entry(inner, "Rogowski half-width", "rogw", 6.0, "mm")
        add_entry(inner, "Rogowski edge radius", "rogr", 2.0, "mm")
        mf = ttk.Frame(inner); mf.pack(fill=tk.X, pady=1)
        ttk.Label(mf, text="Plate material", width=17, font=FONT_UI).pack(side=tk.LEFT)
        self.mat = tk.StringVar(value="Metal (grounded)")
        ttk.Combobox(mf, textvariable=self.mat, width=12, state="readonly",
                     values=list(PLATE_PRESET.keys()), font=FONT_UI).pack(side=tk.LEFT)
        add_entry(inner, "Plate thickness", "pthick", 1.0, "mm")
        add_entry(inner, "Custom εr", "pepsr", 0.0, "")
        add_entry(inner, "Gas pressure", "press", 1.0, "atm")
        add_entry(inner, "N₂ fraction", "n2", 78.0, "%")
        add_entry(inner, "O₂ fraction", "o2", 21.0, "%")
        add_entry(inner, "Ar fraction", "ar", 1.0, "%")
        add_entry(inner, "H₂O fraction", "h2o", 0.0, "%")
        self.lbl_N = ttk.Label(inner, text="", foreground="#0055aa", font=FONT_UI)
        self.lbl_N.pack(anchor="w", padx=4)

        ttk.Separator(inner).pack(fill=tk.X, pady=3)
        ttk.Label(inner, text="Source / Rate model", font=FONT_UI_B).pack(anchor="w")
        add_entry(inner, "Peak voltage Vp", "volt", 30.0, "kV")
        pf = ttk.Frame(inner); pf.pack(fill=tk.X, pady=1)
        ttk.Label(pf, text="Needle polarity", width=17, font=FONT_UI).pack(side=tk.LEFT)
        self.pol = tk.IntVar(value=+1)
        ttk.Radiobutton(pf, text="＋", variable=self.pol, value=+1).pack(side=tk.LEFT)
        ttk.Radiobutton(pf, text="－", variable=self.pol, value=-1).pack(side=tk.LEFT)
        wf = ttk.Frame(inner); wf.pack(fill=tk.X, pady=1)
        ttk.Label(wf, text="Waveform", width=17, font=FONT_UI).pack(side=tk.LEFT)
        self.wave = tk.StringVar(value="IMPULSE")
        ttk.Radiobutton(wf, text="Impulse", variable=self.wave, value="IMPULSE").pack(side=tk.LEFT)
        ttk.Radiobutton(wf, text="DC", variable=self.wave, value="DC").pack(side=tk.LEFT)
        add_entry(inner, "Front time T_f", "tfront", 0.01, "µs")
        add_entry(inner, "Tail time T_t", "ttail", 0.1, "µs")
        rf = ttk.Frame(inner); rf.pack(fill=tk.X, pady=1)
        ttk.Label(rf, text="Rate model", width=17, font=FONT_UI).pack(side=tk.LEFT)
        self.rmodel = tk.StringVar(value="TABLE")
        ttk.Radiobutton(rf, text="Cross-section", variable=self.rmodel, value="TABLE").pack(side=tk.LEFT)
        ttk.Radiobutton(rf, text="Townsend", variable=self.rmodel, value="TOWNSEND").pack(side=tk.LEFT)

        ttk.Separator(inner).pack(fill=tk.X, pady=3)
        ttk.Label(inner, text="Discharge physics / Numerics", font=FONT_UI_B).pack(anchor="w")
        add_entry(inner, "SEE coefficient γ", "gamma", 0.05, "")
        add_entry(inner, "Cathode emission J₀", "j0", "1e-2", "A/m²")
        add_entry(inner, "Background O₂⁻", "nionbg", "1e9", "m⁻³")
        add_entry(inner, "Background rate q", "qbg", "1e7", "m⁻³s⁻¹")
        add_entry(inner, "Manual dt (0=auto)", "dtns", 0.0, "ns")
        add_entry(inner, "Snapshot interval", "snap", 1.0, "ns")
        add_entry(inner, "Grid points nx", "nx", 161, "")
        add_entry(inner, "Background nₑ", "nbg", "0", "m⁻³")
        add_entry(inner, "Seed density", "seed", "0", "m⁻³")
        add_entry(inner, "Substeps per frame", "nsub", 6, "")

        ttk.Separator(inner).pack(fill=tk.X, pady=3)
        ttk.Label(inner, text="Status", font=FONT_UI_B).pack(anchor="w")
        self.status = tk.Text(inner, width=44, height=20, font=("Times New Roman", 12),
                              state="disabled", relief="flat", background="#f4f4f4")
        self.status.pack(fill=tk.X, pady=2)

    # ------------------------------------------------------------------
    def _build_tabs(self):
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        def make_tab(title, figsize):
            frame = ttk.Frame(self.nb)
            self.nb.add(frame, text=title)
            fig = Figure(figsize=figsize, dpi=100)
            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
            return frame, fig, canvas

        _, self.figE, self.cvE = make_tab("① 전계", (9.8, 8.8))
        _, self.figN, self.cvN = make_tab("② 전자밀도", (9.8, 8.8))
        _, self.figZ, self.cvZ = make_tab("③ 선단확대", (9.8, 8.8))
        _, self.figP, self.cvP = make_tab("④ 축상", (9.8, 9.8))
        frR, self.figR, self.cvR = make_tab("⑥ 라디칼", (9.8, 9.4))
        _, self.figW, self.cvW = make_tab("⑤ 파형", (9.8, 8.8))
        _, self.figS, self.cvS = make_tab("⑦ 스트릭/이력", (10.4, 9.6))
        _, self.figT, self.cvT = make_tab("⑧ 개시통계", (9.8, 9.0))

        # 라디칼 탭 상단: 종 선택
        top = ttk.Frame(frR); top.pack(side=tk.TOP, fill=tk.X)
        ttk.Label(top, text="Species:", font=FONT_UI).pack(side=tk.LEFT, padx=6)
        self.rad_sel = tk.StringVar(value="O(P)")
        cb = ttk.Combobox(top, textvariable=self.rad_sel, width=6, state="readonly",
                          values=RAD_SPECIES, font=FONT_UI)
        cb.pack(side=tk.LEFT)
        cb.bind("<<ComboboxSelected>>", lambda e: self._refresh_plots())

        self.axE = self.figE.add_subplot(111)
        self.axN = self.figN.add_subplot(111)
        self.axZ = self.figZ.add_subplot(111)
        gs = self.figP.add_gridspec(4, 1, hspace=0.55)
        self.axP1 = self.figP.add_subplot(gs[0])
        self.axP2 = self.figP.add_subplot(gs[1])
        self.axP3 = self.figP.add_subplot(gs[2])
        self.axP4 = self.figP.add_subplot(gs[3])
        gr = self.figR.add_gridspec(2, 1, height_ratios=[1.5, 1.0], hspace=0.35)
        self.axR1 = self.figR.add_subplot(gr[0])   # 2D 맵
        self.axR2 = self.figR.add_subplot(gr[1])   # 축상 5종 중첩
        self.axW = self.figW.add_subplot(111)
        gsS = self.figS.add_gridspec(2, 3, height_ratios=[1.35, 1.0],
                                     wspace=0.55, hspace=0.40)
        self.axS1 = self.figS.add_subplot(gsS[0, 0])   # (a) E/N [Td]
        self.axS2 = self.figS.add_subplot(gsS[0, 1])   # (b) ne [cm^-3]
        self.axS3 = self.figS.add_subplot(gsS[0, 2])   # (c) 스트릭 z-t
        self.axS4 = self.figS.add_subplot(gsS[1, :])   # (d) 활성종 총개수
        self.cbS1 = self.cbS2 = self.cbS3 = None
        gsT = self.figT.add_gridspec(2, 1, hspace=0.42)
        self.axT1 = self.figT.add_subplot(gsT[0])   # Lambda(t), F(t)
        self.axT2 = self.figT.add_subplot(gsT[1])   # 지연시간 히스토그램
        self.imE = self.imN = self.imZ = self.imR = None
        self.cbE = self.cbN = self.cbZ = self.cbR = None
        self.ctr = None

    # ------------------------------------------------------------------
    def _read(self, key, cast=float):
        return cast(float(self.vars[key].get()))

    def load_rate_csv(self):
        """BOLSIG+ 반응률 CSV 로드 (기본 폴더: 바탕화면/HighVoltage)"""
        p = filedialog.askopenfilename(title="Select rate CSV (header: EN_Td,...)",
                                       initialdir=self.hv_dir,
                                       filetypes=[("CSV", "*.csv"), ("All", "*.*")])
        if not p:
            return
        try:
            ov, loaded = parse_rate_csv(p)
            self.rate_override = ov
            self.lbl_rate.config(
                text=f"Rates: {os.path.basename(p)} ({len(loaded)} replaced)")
            messagebox.showinfo("Loaded",
                f"Replaced rates: {', '.join(loaded)}\nApplied on next Apply/Reset.")
        except Exception as e:
            messagebox.showerror("CSV error", str(e))

    def apply_params(self):
        try:
            self.running = False
            self.btn_run.config(text="▶ Run")
            self.btn_ag.config(text="🌙 Afterglow")
            mode, epsr, color = PLATE_PRESET[self.mat.get()]
            if epsr is None:
                epsr = self._read("pepsr")
            self.plate_color = color
            self.sim = NeedlePlaneDischarge(
                gap_mm=self._read("gap"), voltage_kV=self._read("volt"),
                polarity=self.pol.get(), pressure_atm=self._read("press"),
                r_main_mm=self._read("rmain"), r_tip_mm=self._read("rtip"),
                taper_mm=self._read("taper"), nx=self._read("nx", int),
                n_background=self._read("nbg"), seed_density=self._read("seed"),
                wave_mode=self.wave.get(), t_front_us=self._read("tfront"),
                t_tail_us=self._read("ttail"),
                n2_pct=self._read("n2"), o2_pct=self._read("o2"),
                ar_pct=self._read("ar"),
                h2o_pct=self._read("h2o"), gamma_see=self._read("gamma"),
                dt_user_ns=self._read("dtns"), snap_ns=self._read("snap"),
                plate_mode=mode, plate_thick_mm=self._read("pthick"),
                plate_epsr=epsr, rate_model=self.rmodel.get(),
                rate_override=self.rate_override,
                cath_J0=self._read("j0"),
                n_ion_bg=self._read("nionbg"),
                q_bg=self._read("qbg"),
                config=CONFIG_MAP[self.cfg.get()],
                needle_sep_mm=self._read("nsep"),
                sphere_R_mm=self._read("srad"),
                rog_halfw_mm=self._read("rogw"),
                rog_edge_mm=self._read("rogr"))
            s = self.sim
            self.lbl_N.config(text=(
                f"N = {s.Ngas:.3e} m⁻³\n"
                f"   N₂/O₂/Ar/H₂O = {s.x_n2*100:.1f}/{s.x_o2*100:.1f}/"
                f"{s.x_ar*100:.1f}/{s.x_h2o*100:.1f} % (normalized)\n"
                f"   Plate: {self.mat.get()}, t={s.t_solid*1e3:.2f} mm"))
            if s.rtip < s.h:
                messagebox.showwarning("Grid resolution",
                    f"Tip radius ({s.rtip*1e3:.3f} mm) < grid spacing ({s.h*1e3:.3f} mm): "
                    f"tip approximated by one cell. Increase nx for accuracy.")
            self.frame_count = 0
            self._init_plots()
            self._update_status()
        except Exception as e:
            messagebox.showerror("Input error", str(e))

    # ------------------------------------------------------------------
    def _fmt_t(self, t):
        return f"{t*1e9:.2f} ns" if t < 1e-6 else f"{t*1e6:.3f} µs"

    def _init_plots(self):
        s = self.sim
        self.ext = [s.x[0]*1e3, s.x[-1]*1e3, s.y[0]*1e3, s.y[-1]*1e3]
        for cb in (self.cbE, self.cbN, self.cbZ, self.cbR):
            if cb:
                cb.remove()
        self.axE.clear(); self.axN.clear(); self.axZ.clear(); self.axR1.clear()

        Ekv = s.Emag_masked() / 1e5
        vmaxE = max(np.nanmax(Ekv), s.E0max_peak/1e5) * 1.05
        self.imE = self.axE.imshow(Ekv, origin="lower", extent=self.ext,
                                   cmap="inferno", vmin=0, vmax=vmaxE, aspect="equal")
        self.cbE = self.figE.colorbar(self.imE, ax=self.axE, fraction=0.046,
                                      label="|E| [kV/cm]")
        self.ctr = None
        self._draw_solids(self.axE)
        self.axE.set_title("Electric field |E| and equipotential lines")
        self.axE.set_xlabel("x [mm]"); self.axE.set_ylabel("z [mm]  (0 = bottom electrode)")

        self.imN = self.axN.imshow(np.maximum(s.n_e, 1e6), origin="lower",
                                   extent=self.ext, cmap="viridis",
                                   norm=LogNorm(vmin=max(s.n_bg, 1e8), vmax=1e21),
                                   aspect="equal")
        self.cbN = self.figN.colorbar(self.imN, ax=self.axN, fraction=0.046,
                                      label="$n_e$ [m$^{-3}$] (log)")
        self._draw_solids(self.axN)
        self.axN.set_title("Electron density $n_e$")
        self.axN.set_xlabel("x [mm]"); self.axN.set_ylabel("z [mm]  (0 = bottom electrode)")

        self.zx0, self.zx1 = -0.30*s.gap, 0.30*s.gap
        self.zy0 = max(s.t_solid, s.y_tip - 0.45*s.gap)
        self.zy1 = min(s.height, s.y_tip + 0.20*s.gap)
        self.zi0 = np.searchsorted(s.x, self.zx0)
        self.zi1 = np.searchsorted(s.x, self.zx1) + 1
        self.zj0 = np.searchsorted(s.y, self.zy0)
        self.zj1 = np.searchsorted(s.y, self.zy1) + 1
        zext = [s.x[self.zi0]*1e3, s.x[self.zi1-1]*1e3,
                s.y[self.zj0]*1e3, s.y[self.zj1-1]*1e3]
        zne = np.maximum(s.n_e[self.zj0:self.zj1, self.zi0:self.zi1], 1e6)
        self.imZ = self.axZ.imshow(zne, origin="lower", extent=zext,
                                   cmap="plasma",
                                   norm=LogNorm(vmin=max(s.n_bg, 1e8), vmax=1e21),
                                   aspect="equal")
        self.cbZ = self.figZ.colorbar(self.imZ, ax=self.axZ, fraction=0.046,
                                      label="$n_e$ [m$^{-3}$] (log)")
        zn = s.needle[self.zj0:self.zj1, self.zi0:self.zi1]
        self.axZ.contourf(s.X[self.zj0:self.zj1, self.zi0:self.zi1]*1e3,
                          s.Y[self.zj0:self.zj1, self.zi0:self.zi1]*1e3,
                          zn.astype(float), levels=[0.5, 1.5], colors=["#909090"])
        self.axZ.set_title("Zoomed $n_e$ near needle tip (SEE included)")
        self.axZ.set_xlabel("x [mm]"); self.axZ.set_ylabel("z [mm]  (0 = bottom electrode)")

        # 라디칼 2D 맵
        self.imR = self.axR1.imshow(np.maximum(s.rad["O(P)"], 1e6), origin="lower",
                                    extent=self.ext, cmap="cividis",
                                    norm=LogNorm(vmin=1e14, vmax=1e23),
                                    aspect="equal")
        self.cbR = self.figR.colorbar(self.imR, ax=self.axR1, fraction=0.046,
                                      label="Radical density [m$^{-3}$] (log)")
        self._draw_solids(self.axR1)
        self.axR1.set_xlabel("x [mm]"); self.axR1.set_ylabel("z [mm]  (0 = bottom electrode)")

        self._refresh_plots()

    def _draw_solids(self, ax):
        s = self.sim
        ax.contourf(s.X*1e3, s.Y*1e3, s.needle.astype(float),
                    levels=[0.5, 1.5], colors=["#909090"])
        if s.gnd.any():
            ax.contourf(s.X*1e3, s.Y*1e3, s.gnd.astype(float),
                        levels=[0.5, 1.5], colors=["#6a6a6a"])
        if s.n_solid >= 1:
            ax.axhspan(0, s.t_solid*1e3, color=self.plate_color, alpha=0.85)

    # ------------------------------------------------------------------
    def _refresh_paperfig(self):
        """⑦ 탭: 첨부 논문 그림형 (a)E/N (b)ne (c)스트릭 (d)활성종 이력"""
        s = self.sim
        for cb in (self.cbS1, self.cbS2, self.cbS3):
            if cb:
                try: cb.remove()
                except Exception: pass
        self.cbS1 = self.cbS2 = self.cbS3 = None
        for ax in (self.axS1, self.axS2, self.axS3, self.axS4):
            ax.clear()

        # 기체 갭 영역: z = (격자 y) - (평판 표면 y)  [물리좌표, 전극배치와 일치]
        #   z = 0 : 평판 표면(하단),  z = gap : 바늘 선단(상단)
        #   평판 두께 0(접지)일 때는 표면이 0행이므로 0행부터 포함 -> z=0 정확 일치
        j0 = s.j_ax_lo if (s.gnd.any() or s.n_solid >= 1) else 0
        j1 = max(s.j_ax_hi, j0 + 1)
        gapmm = s.gap * 1e3                                   # 실제 갭 [mm]
        hm = s.h * 1e3
        z_hi = (s.y[j1] - s.y_bot) * 1e3                    # 최상단 행 실제 z(≤gap)
        # imshow extent는 픽셀 '가장자리' 기준 (반 셀 보정) -> 각 행이 실제 z에 위치
        ext = [s.x[0]*1e3, s.x[-1]*1e3, -0.5*hm, z_hi + 0.5*hm]
        xr = min(0.20 * s.gap * 1e3, s.x[-1]*1e3)
        sub = slice(j0, j1 + 1)

        # ---- (a) E/N [Td] ----
        EN = (s.Emag / (s.Ngas * TD))[sub]
        EN = np.where(s.needle[sub] | s.gnd[sub], np.nan, EN)
        im1 = self.axS1.imshow(EN, origin="lower", extent=ext, aspect="auto",
                               cmap="jet", vmin=0,
                               vmax=max(np.nanmax(EN)*1.02, 10.0))
        self.cbS1 = self.figS.colorbar(im1, ax=self.axS1, fraction=0.10)
        self.axS1.set_title("(a) E/N [Td]", fontsize=16)
        self.axS1.set_xlim(-xr, xr)
        self.axS1.set_ylim(0.0, gapmm)                        # z: 0(평판)~gap(선단)
        self.axS1.set_xlabel("r [mm]", fontsize=16)
        self.axS1.set_ylabel("z [mm]  (0 = bottom electrode)", fontsize=16)

        # ---- (b) 전자밀도 [cm^-3] ----
        nec = np.maximum(s.n_e[sub] * 1e-6, 1e6)
        im2 = self.axS2.imshow(nec, origin="lower", extent=ext, aspect="auto",
                               cmap="jet", norm=LogNorm(vmin=1e9, vmax=1e16))
        self.cbS2 = self.figS.colorbar(im2, ax=self.axS2, fraction=0.10)
        self.axS2.set_title("(b) $n_e$ [cm$^{-3}$]", fontsize=16)
        self.axS2.set_xlim(-xr, xr)
        self.axS2.set_ylim(0.0, gapmm)
        self.axS2.set_xlabel("r [mm]", fontsize=16)
        self.axS2.set_ylabel("z [mm]  (0 = bottom electrode)", fontsize=16)

        # ---- (c) 발광 스트릭 (z-t, N2 SPS 여기율 ∝ 발광강도) ----
        if len(s.streak_em) >= 3:
            tt = np.array([p[0] for p in s.streak_em]) * 1e9     # [ns]
            M = np.array([np.maximum(p[1][sub], 1e-30) for p in s.streak_em])
            M = np.log10(np.maximum(M.T, 1e-30))                 # (z, t)
            vhi = M.max()
            zed = (s.y[j0:j1 + 1] - s.y_bot) * 1e3             # 행별 실제 z [mm]
            pm = self.axS3.pcolormesh(*np.meshgrid(tt, zed),
                                      M, cmap="hot", shading="nearest",
                                      vmin=vhi - 7.0, vmax=vhi)
            self.cbS3 = self.figS.colorbar(pm, ax=self.axS3, fraction=0.10,
                                           label="log$_{10}$ emission [a.u.]")
            self.axS3.set_ylim(0.0, gapmm)                       # z: 0(평판)~gap(선단)
        self.axS3.set_title("(c) Optical emission streak (Simulation)", fontsize=16)
        self.axS3.set_xlabel("Time [ns]", fontsize=16)
        self.axS3.set_ylabel("z [mm]  (0 = bottom electrode)", fontsize=16)

        # ---- (d) 활성종 총개수 이력 (log-log) ----
        if len(s.rad_hist) >= 3:
            tt = np.array([p[0] for p in s.rad_hist]) * 1e9
            m = tt > 0
            for sp in RAD_SPECIES:
                vals = np.array([p[1][sp] for p in s.rad_hist])
                c_, ls_, mk_ = RAD_STYLE[sp]
                npt = int(m.sum())
                self.axS4.plot(tt[m], np.maximum(vals[m], 1e6),
                               color=c_, ls=ls_, marker=mk_,
                               markevery=max(npt//10, 1), ms=6, mfc="none",
                               lw=1.8, label=sp)
            self.axS4.set_xscale("log"); self.axS4.set_yscale("log")
            self.axS4.set_xlim(max(tt[m][0], 1.0), max(tt[-1]*1.1, 10.0))
            vmax_all = max(max(p[1][sp] for sp in RAD_SPECIES)
                           for p in s.rad_hist)
            self.axS4.set_ylim(max(vmax_all*1e-5, 1e8), vmax_all*3)
            self.axS4.legend(loc="lower right", ncol=6, framealpha=0.6)
        gases = (f"H₂O({s.nH2O/s.Ngas*100:.0f}%)/" if s.nH2O > 0 else "") + "O₂/N₂"
        self.axS4.set_title(f"(d) Total number of radicals (axisymmetric estimate)  [{gases}]",
                            fontsize=16)
        self.axS4.set_xlabel("Time [ns]", fontsize=16)
        self.axS4.set_ylabel("Total number of radicals", fontsize=16)
        self.axS4.grid(alpha=0.3, which="both")
        self.cvS.draw_idle()

    # ------------------------------------------------------------------
    def _refresh_plots(self):
        s = self.sim
        if s is None:
            return
        pol = "+" if s.polarity > 0 else "−"
        tstr = self._fmt_t(s.t)

        Ekv = s.Emag_masked() / 1e5
        self.imE.set_data(Ekv)
        self.imE.set_clim(0, max(np.nanmax(Ekv)*1.02, 1e-3))
        if self.ctr is not None:
            try:
                self.ctr.remove()
            except Exception:
                for c in getattr(self.ctr, "collections", []):
                    try: c.remove()
                    except Exception: pass
        Vk = s.V.copy(); Vk[s.needle] = np.nan
        self.ctr = self.axE.contour(s.X*1e3, s.Y*1e3, Vk/1e3, 12,
                                    colors="w", linewidths=0.5, alpha=0.65)
        extra = " [afterglow: E frozen]" if s.afterglow else ""
        if s.plate_mode == "FLOAT_METAL":
            extra += f", V_plate={s.V_plate/1e3:.2f} kV"
        elif s.plate_mode == "DIELECTRIC":
            extra += f", σ_max={np.abs(s.sigma).max()*1e6:.2f} µC/m²"
        self.figE.suptitle(f"V(t)={s.Va_now/1e3:+.2f} kV (Vp={pol}{abs(s.Vpeak)/1e3:.1f} kV), "
                           f"t={tstr}{extra}", fontsize=16)

        self.imN.set_data(np.maximum(s.n_e, 1e6))
        self.figN.suptitle(f"max $n_e$={s.n_e.max():.2e} m$^{{-3}}$, t={tstr}", fontsize=16)
        self.imZ.set_data(np.maximum(s.n_e[self.zj0:self.zj1, self.zi0:self.zi1], 1e6))
        self.figZ.suptitle(f"gamma={s.gamma_see:g}, t={tstr}", fontsize=16)

        # ---- 라디칼 탭 ----
        sp = self.rad_sel.get()
        self.imR.set_data(np.maximum(s.rad[sp], 1e6))
        self.axR1.set_title(f"2D distribution of {sp}  (max {s.rad[sp].max():.2e} m$^{{-3}}$)")
        self.figR.suptitle(f"Radical distributions, t={tstr}"
                           f"{' [afterglow]' if s.afterglow else ''}", fontsize=16)
        ic = s.nx // 2
        yy = s.y * 1e3
        xmax = min((s.t_solid + 1.2 * s.gap) * 1e3, s.height * 1e3)
        self.axR2.clear()
        for name in RAD_SPECIES:
            c_, ls_, mk_ = RAD_STYLE[name]
            self.axR2.semilogy(yy, np.maximum(s.rad[name][:, ic], 1e6),
                               color=c_, ls=ls_, marker=mk_,
                               markevery=max(len(yy)//8, 1), ms=6, mfc="none",
                               lw=1.8, label=name)
        self.axR2.set_xlim(0, xmax)
        self.axR2.set_ylim(1e12, 1e24)
        self.axR2.set_xlabel("Axial position z [mm]  (0 = bottom electrode)")
        self.axR2.set_ylabel("Density [m$^{-3}$]")
        self.axR2.set_title("On-axis radical density profiles (current time)")
        self.axR2.grid(alpha=0.3, which="both")
        self.axR2.legend(loc="upper left", ncol=5)
        self.axR2.axvline(s.y_tip*1e3, color="#888", ls="--", lw=1.0)
        if s.n_solid >= 1:
            self.axR2.axvspan(0, s.t_solid*1e3, color=self.plate_color, alpha=0.4)

        # ---- 축상 프로파일 (하전종) ----
        snaps = s.profiles[-self.MAX_OVERLAY:]
        colors = cm.turbo(np.linspace(0.1, 0.95, len(snaps)))
        for ax in (self.axP1, self.axP2, self.axP3, self.axP4):
            ax.clear()
        for i, ((tsn, E_ax, ne_ax, nip_ax, rho_ax), c) in enumerate(zip(snaps, colors)):
            lbl = f"t={tsn*1e9:.1f} ns"
            ls = LSTYLES[i % len(LSTYLES)]
            mk = MARKERS[i % len(MARKERS)]
            kw = dict(color=c, lw=1.6, ls=ls, marker=mk,
                      markevery=max(len(yy)//9, 1), ms=6, mfc="none")
            self.axP1.plot(yy, E_ax/1e5, label=lbl, **kw)
            self.axP2.semilogy(yy, np.maximum(ne_ax, 1e6), **kw)
            self.axP3.semilogy(yy, np.maximum(nip_ax, 1e6), **kw)
            self.axP4.plot(yy, rho_ax, **kw)
        self.axP1.axhline(30*s.p_atm, color="k", ls=":", lw=1.0)
        self.axP1.set_ylabel("|E| [kV/cm]")
        self.axP1.set_title(f"On-axis electric field (overlaid every {s.snap_dt*1e9:g} ns; legend shared)")
        self.axP1.legend(loc="upper left", ncol=2, framealpha=0.6)
        self.axP2.set_ylabel("$n_e$ [m$^{-3}$]")
        self.axP2.set_title("On-axis electron density")
        self.axP3.set_ylabel("$n_+$ [m$^{-3}$]")
        self.axP3.set_title("On-axis positive-ion density")
        self.axP4.set_ylabel(r"$\rho$ [C/m$^3$]")
        self.axP4.set_yscale("symlog", linthresh=1e-3)
        # symlog 음수 눈금의 유니코드 마이너스(U+2212) 깨짐 방지: ASCII 포맷터
        def _ascii_fmt(v, _p):
            if v == 0:
                return "0"
            m, e = f"{v:.0e}".split("e")
            return f"{m}e{int(e)}"
        self.axP4.yaxis.set_major_formatter(FuncFormatter(_ascii_fmt))
        self.axP4.axhline(0, color="k", lw=0.7)
        self.axP4.set_title(r"On-axis net charge density  $\rho = q(n_+ - n_e - n_-)$")
        self.axP4.set_xlabel("Axial position z [mm]  (0 = plate, needle tip at z = gap)")
        for ax in (self.axP1, self.axP2, self.axP3, self.axP4):
            ax.set_xlim(0, xmax)
            ax.grid(alpha=0.3)
            ax.axvline(s.y_tip*1e3, color="#888", ls="--", lw=1.0)
            if s.n_solid >= 1:
                ax.axvspan(0, s.t_solid*1e3, color=self.plate_color, alpha=0.4)

        # ---- 파형 ----
        self.axW.clear()
        tt, vv = s.waveform_curve()
        unit, sc = ("µs", 1e6) if tt[-1] > 3e-6 else ("ns", 1e9)
        self.axW.plot(tt*sc, vv/1e3, "g-", lw=2.0, label="Applied voltage")
        if s.plate_mode == "FLOAT_METAL":
            self.axW.axhline(s.V_plate/1e3, color="m", ls="--", lw=1.5,
                             label=f"Floating-plate potential {s.V_plate/1e3:.2f} kV")
            self.axW.legend()
        tm = min(s.t, tt[-1])
        self.axW.plot(tm*sc, s.Va_now/1e3, "ro", ms=9)
        if s.wave_mode != "DC":
            self.axW.set_title(f"Applied voltage waveform, front {s.T_front*1e6:g} / "
                               f"tail {s.T_tail*1e6:g} us")
        else:
            self.axW.set_title("Applied voltage waveform (DC)")
        self.axW.set_xlabel(f"Time [{unit}]"); self.axW.set_ylabel("Voltage [kV]")
        self.axW.grid(alpha=0.3)

        self._refresh_paperfig()
        self._refresh_initstats()
        for cv in (self.cvE, self.cvN, self.cvZ, self.cvP, self.cvR, self.cvW):
            cv.draw_idle()

    # ------------------------------------------------------------------
    def _refresh_initstats(self):
        """⑧ 탭: 개시율 Lambda(t)/개시확률 F(t) 및 통계적 지연시간 분포"""
        s = self.sim
        self.axT1.clear(); self.axT2.clear()
        if len(s.init_stats) >= 3:
            tt = np.array([p[0] for p in s.init_stats]) * 1e9
            Lam = np.array([p[1] for p in s.init_stats])
            F = np.array([p[3] for p in s.init_stats])
            self.axT1.semilogy(tt, np.maximum(Lam, 1e-3), "b-", lw=1.8,
                               label=r"$\Lambda(t)$ [1/s]")
            self.axT1.set_ylabel(r"$\Lambda(t)$ [1/s]", color="b")
            ax2 = self.axT1.twinx()
            ax2.plot(tt, F, "r-", lw=1.8, label="F(t)")
            ax2.set_ylabel("Initiation probability F(t)", color="r")
            ax2.set_ylim(0, 1.02)
            self.axT1.set_title(
                r"Effective initiation rate $\Lambda=\int_{\alpha_{eff}>0}"
                r"\nu_{det}\,n_{O_2^-}\,dV$  and  $F(t)=1-e^{-\int\Lambda dt}$")
        self.axT1.set_xlabel("Time [ns]")
        self.axT1.grid(alpha=0.3, which="both")

        samples, cens = s.sample_delays(2000)
        if samples.size >= 5:
            self.axT2.hist(samples * 1e9, bins=40, color="#4477aa",
                           edgecolor="k", linewidth=0.4)
            med = np.median(samples) * 1e9
            self.axT2.axvline(med, color="r", ls="--", lw=1.5)
            self.axT2.set_title(
                f"Statistical delay-time distribution "
                f"(N=2000, median={med:.2f} ns, censored beyond t_end: {cens*100:.1f}%)")
        else:
            self.axT2.set_title(
                "Statistical delay-time distribution "
                f"(insufficient hazard yet; F(t_end)={0 if not s.init_stats else s.init_stats[-1][3]*100:.2f}%)")
        self.axT2.set_xlabel("Delay time [ns]")
        self.axT2.set_ylabel("Counts")
        self.axT2.grid(alpha=0.3)
        self.cvT.draw_idle()

    # ------------------------------------------------------------------
    def export_csv(self):
        if self.sim is None:
            return
        d = filedialog.askdirectory(title="Select folder to save CSV",
                                    initialdir=self.hv_dir)
        if not d:
            return
        s = self.sim
        stamp = time.strftime("%Y%m%d_%H%M%S")
        try:
            path1 = f"{d}/axial_profiles_{stamp}.csv"
            with open(path1, "w", encoding="utf-8-sig") as f:
                f.write("t_ns,y_mm,E_kV_per_cm,ne_m-3,nip_m-3,rho_C_per_m3\n")
                for (tsn, E_ax, ne_ax, nip_ax, rho_ax) in s.profiles:
                    for j in range(s.ny):
                        e = "" if np.isnan(E_ax[j]) else f"{E_ax[j]/1e5:.6e}"
                        f.write(f"{tsn*1e9:.4f},{s.y[j]*1e3:.4f},{e},"
                                f"{ne_ax[j]:.6e},{nip_ax[j]:.6e},{rho_ax[j]:.6e}\n")

            ic = s.nx // 2
            with open(f"{d}/axial_radicals_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                f.write("y_mm," + ",".join(f"{sp}_m-3" for sp in RAD_SPECIES) + "\n")
                for j in range(s.ny):
                    f.write(f"{s.y[j]*1e3:.4f}," +
                            ",".join(f"{s.rad[sp][j, ic]:.6e}" for sp in RAD_SPECIES) + "\n")

            def save2d(name, arr, fmt="%.6e"):
                with open(f"{d}/{name}_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                    f.write("y_mm\\x_mm," + ",".join(f"{v*1e3:.4f}" for v in s.x) + "\n")
                    for j in range(s.ny):
                        f.write(f"{s.y[j]*1e3:.4f}," +
                                ",".join(fmt % arr[j, i] for i in range(s.nx)) + "\n")
            save2d("field2D_Emag_Vpm", np.nan_to_num(s.Emag_masked()))
            save2d("field2D_ne_m-3",  s.n_e)
            save2d("field2D_nip_m-3", s.n_ip)
            save2d("field2D_rho_Cpm3", QE*(s.n_ip - s.n_e - s.n_in))
            for sp in RAD_SPECIES:
                save2d(f"field2D_{sp}_m-3", s.rad[sp])

            with open(f"{d}/radical_totals_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                f.write("t_ns," + ",".join(RAD_SPECIES) + "\n")
                for (tsn, tot) in s.rad_hist:
                    f.write(f"{tsn*1e9:.5f}," +
                            ",".join(f"{tot[sp]:.6e}" for sp in RAD_SPECIES) + "\n")
            with open(f"{d}/streak_ne_axis_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                f.write("y_mm\\t_ns," +
                        ",".join(f"{p[0]*1e9:.4f}" for p in s.streak) + "\n")
                for j in range(s.ny):
                    f.write(f"{s.y[j]*1e3:.4f}," +
                            ",".join(f"{p[1][j]:.4e}" for p in s.streak) + "\n")
            with open(f"{d}/streak_emission_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                f.write("y_mm\\t_ns," +
                        ",".join(f"{p[0]*1e9:.4f}" for p in s.streak_em) + "\n")
                for j in range(s.ny):
                    f.write(f"{s.y[j]*1e3:.4f}," +
                            ",".join(f"{p[1][j]:.4e}" for p in s.streak_em) + "\n")

            with open(f"{d}/initiation_stats_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                f.write("t_ns,Lambda_per_s,H_cumhazard,F_probability\n")
                for (tsn, lam, hh, ff) in s.init_stats:
                    f.write(f"{tsn*1e9:.5f},{lam:.6e},{hh:.6e},{ff:.6e}\n")

            tt, vv = s.waveform_curve(800)
            with open(f"{d}/waveform_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                f.write("t_s,V_applied_V\n")
                for a_, b_ in zip(tt, vv):
                    f.write(f"{a_:.6e},{b_:.6e}\n")

            with open(f"{d}/conditions_{stamp}.csv", "w", encoding="utf-8-sig") as f:
                f.write("parameter,value,unit\n")
                rows = [("gap", s.gap*1e3, "mm"), ("Vpeak", s.Vpeak/1e3, "kV"),
                        ("polarity", s.polarity, "-"), ("pressure", s.p_atm, "atm"),
                        ("N_gas", s.Ngas, "m-3"),
                        ("x_N2_auto", s.x_n2, "-"), ("x_O2", s.x_o2, "-"),
                        ("x_Ar", s.x_ar, "-"), ("x_H2O", s.x_h2o, "-"),
                        ("rate_model", s.rate_model, "-"),
                        ("r_main", s.r_main*1e3, "mm"), ("r_tip", s.rtip*1e3, "mm"),
                        ("plate_mode", s.plate_mode, "-"),
                        ("plate_thickness", s.t_solid*1e3, "mm"),
                        ("plate_epsr", s.plate_epsr, "-"),
                        ("V_plate_now", s.V_plate, "V"),
                        ("sigma_max", float(np.abs(s.sigma).max()), "C/m2"),
                        ("gamma_SEE", s.gamma_see, "-"),
                        ("wave_mode", s.wave_mode, "-"),
                        ("T_front", s.T_front*1e6, "us"), ("T_tail", s.T_tail*1e6, "us"),
                        ("t_now", s.t, "s"), ("afterglow", s.afterglow, "-"),
                        ("nx", s.nx, "-"), ("ny", s.ny, "-"),
                        ("h_grid", s.h*1e3, "mm")]
                for sp in RAD_SPECIES:
                    rows.append((f"max_{sp}", float(s.rad[sp].max()), "m-3"))
                for k, v, u in rows:
                    f.write(f"{k},{v},{u}\n")

            messagebox.showinfo("Saved",
                f"CSV files saved to: {d}\n"
                f"- axial_profiles / axial_radicals\n"
                f"- field2D_(Emag/ne/nip/rho/{'/'.join(RAD_SPECIES)})\n"
                f"- radical_totals / streak_ne_axis / streak_emission\n"
                f"- initiation_stats\n"
                f"- waveform / conditions")
        except Exception as e:
            messagebox.showerror("Save error", str(e))

    # ------------------------------------------------------------------
    def _update_status(self):
        s = self.sim
        dtu = s.dt_last*1e9
        dtinfo = f"{dtu:9.4f} ns" if dtu < 1e3 else f"{dtu/1e3:9.4f} µs"
        if s.dt_user > 0 and not s.afterglow:
            dtinfo += " (limited)" if s.dt_limited else " (manual)"
        plate = {"GND_METAL": "Metal (grounded)", "FLOAT_METAL": "Metal (floating)",
                 "DIELECTRIC": f"Dielectric εr={s.plate_epsr:g}"}[s.plate_mode]
        mode = "Afterglow" if s.afterglow else "Discharge (streamer)"
        lam_txt = (f" Init. rate Λ(t): {s.init_stats[-1][1]:10.3e} /s\n"
                   f" Init. prob F(t): {s.init_stats[-1][3]*100:10.4f} %\n"
                   if s.init_stats else "")
        cfg_lbl = {v: k for k, v in CONFIG_MAP.items()}.get(s.config, s.config)
        txt = (f" Config        : {cfg_lbl}\n"
               f" Mode          : {mode}\n"
               f" Time t        : {self._fmt_t(s.t):>12s}\n"
               f" dt applied    : {dtinfo}\n"
               f" V(t) applied  : {s.Va_now/1e3:10.2f} kV\n"
               f" Max |E|       : {np.nanmax(s.Emag_masked())/1e5:10.2f} kV/cm\n"
               f" Max nₑ        : {s.n_e.max():10.3e} m⁻³\n"
               f" Max n₊        : {s.n_ip.max():10.3e} m⁻³\n"
               f" Max n₋(O₂⁻)   : {s.n_in.max():10.3e} m⁻³\n"
               + "".join(f" Max {sp:<9s}: {s.rad[sp].max():10.3e} m⁻³\n"
                          for sp in RAD_SPECIES) +
               f" Plate         : {plate}, t={s.t_solid*1e3:.2f} mm\n"
               f" V_plate/σmax  : {s.V_plate/1e3:.3f} kV / "
               f"{np.abs(s.sigma).max()*1e6:.3f} µC/m²\n"
               f" SEE γ / J₀    : {s.gamma_see:g} / {s.J0:g} A/m²\n"
               f" M = e^∫αdz    : {s.M_townsend:10.3e}\n"
               f" γ(M−1)        : {s.gamma_MK:10.3e}"
               f"{' [self-sust ≥1]' if s.gamma_MK >= 1 else ''}\n"
               + lam_txt +
               f" Snapshots     : {len(s.profiles)} (every {s.snap_dt*1e9:g} ns)\n"
               f" Rate model    : "
               f"{'Cross-section tables' if s.rate_model=='TABLE' else 'Townsend approx.'}\n"
               f" Pressure      : {s.p_atm:g} atm  (N={s.Ngas:.3e} m⁻³)\n"
               f" Composition   : N₂ {s.x_n2*100:.0f} / O₂ {s.x_o2*100:.0f} / "
               f"Ar {s.x_ar*100:.0f} / H₂O {s.x_h2o*100:.0f} %\n"
               f" Grid          : {s.nx}x{s.ny} (h={s.h*1e3:.3f} mm)")
        self.status.config(state="normal")
        self.status.delete("1.0", tk.END)
        self.status.insert("1.0", txt)
        self.status.config(state="disabled")

    # ------------------------------------------------------------------
    def toggle_run(self):
        if self.sim is None:
            return
        self.running = not self.running
        self.btn_run.config(text="⏸ Pause" if self.running else "▶ Run")
        if self.running:
            self._loop()

    def toggle_afterglow(self):
        if self.sim is None:
            return
        self.sim.afterglow = not self.sim.afterglow
        self.btn_ag.config(text="⚡ Back to discharge" if self.sim.afterglow
                           else "🌙 Afterglow")
        self._update_status()

    def step_once(self):
        if self.sim is None:
            return
        if self.sim.afterglow:
            self.sim.step_afterglow(n_sub=self._read("nsub", int))
        else:
            self.sim.step(n_sub=self._read("nsub", int))
        self.frame_count += 1
        self._refresh_plots()
        self._update_status()

    def _loop(self):
        if not self.running:
            return
        self.step_once()
        self.root.after(15, self._loop)


# =============================================================================
if __name__ == "__main__":
    root = tk.Tk()
    try:
        st = ttk.Style(); st.theme_use("clam")
        for w in ("TLabel", "TButton", "TRadiobutton", "TEntry",
                  "TCombobox", "TLabelframe.Label", "TLabelframe"):
            st.configure(w, font=FONT_UI)
        root.option_add("*Font", "{Times New Roman} 12")
    except Exception:
        pass
    app = DischargeHMI(root)
    root.mainloop()
