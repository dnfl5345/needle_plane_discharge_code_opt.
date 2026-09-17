# -*- coding: utf-8 -*-
"""
==============================================================================
 고전압공학 방전현상 시뮬레이터 v18 (HMI) - needle_plane_discharge_hmi.py
==============================================================================
 v17 HMI(8개 탭, 다중 전극 구성, 4성분 기체, BOLSIG+ CSV 연동 등)를 그대로 유지하고
 물리 코어를 discharge_core.py(v18, 희소행렬 직접해법)로 교체한 버전.

 [v18 변경사항 - HMI]
  1) 물리 코어 분리: discharge_core.NeedlePlaneDischarge (GUI 없이 run_batch.py 로도 실행)
     - 포아송 직접해법(서브스텝당 65 ms -> ~1 ms), 반응률 LUT, 물리적 유전완화 dt
     - 초기화 6.5 s -> 0.1 s
  2) 화면 갱신 분리: 'Compute per tick'(기본 0.25 s) 동안 연속 계산하고,
     'Redraw interval'(기본 2 s)마다 한 번만 8개 탭을 다시 그림
     (v17: 6 서브스텝마다 전 탭 재렌더링 -> 계산보다 그리기가 오래 걸림)
  3) 갭 브리징 자동 감지: 스트리머가 대향전극에 닿으면 상태창에 표시하고
     'Auto afterglow' 체크 시 잔광 모드로 자동 전환 (v17 은 이 시점 이후 dt 가
     fs 대로 붕괴하여 사실상 진행 불가였음)
  4) 외부 회로 직렬저항 R (Ω) 입력: V_needle = V_src(t) - R·I(t), Sato 전류 I(t)
     -> 파형 탭에 전극전압/방전전류 곡선 추가 (R=0 이면 v17 과 동일한 이상 전압원)
  5) 'Run until t_end' 입력 (0 = 무한), 상태창에 벽시계/속도(ns per s)/dt 제한요인/
     서브스텝 수/전류 표시, CSV 저장은 discharge_export.py 공용 함수 사용
  6) v17 HMI 파일의 나머지 기능·표기는 동일

 [실행]  pip install numpy scipy matplotlib  ->  python needle_plane_discharge_hmi.py
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

from discharge_core import (NeedlePlaneDischarge, parse_rate_csv, PLATE_PRESET, CONFIG_MAP,
                            RAD_SPECIES, RAD_STYLE, TD, QE, HAVE_SCIPY)
from discharge_export import export_csv_files

def _available_font_families(prefs=("Times New Roman", "Nimbus Roman", "Liberation Serif",
                                    "DejaVu Serif")):
    """설치된 글꼴만 rcParams 에 등록 (미설치 글꼴은 matplotlib 이 텍스트마다 경고를 출력해
    Linux 등에서 렌더링이 매우 느려짐; Windows 에서는 Times New Roman 이 그대로 사용됨)"""
    from matplotlib import font_manager as _fm
    names = {f.name for f in _fm.fontManager.ttflist}
    fams = [p for p in prefs if p in names]
    return fams if fams else ["serif"]


plt.rcParams["font.family"] = _available_font_families()
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


# =============================================================================
#  HMI (tkinter GUI)
# =============================================================================
class DischargeHMI:
    MAX_OVERLAY = 10

    def __init__(self, root):
        self.root = root
        root.title("HV Discharge Simulator v18 | fast core (sparse Poisson), circuit R, auto-afterglow")
        root.geometry("1460x950")
        self.sim = None
        self.running = False
        self.frame_count = 0
        self._last_draw = 0.0                 # 마지막 그래프 갱신 벽시계
        self._t_draw = 0.0                    # 마지막 전체 재렌더링 소요시간 [s]
        self._wall_run = 0.0                  # 누적 실행 벽시계 [s]
        self._run_t0 = None
        self._bridge_notified = False
        self._msg = ""                        # 상태창 알림 문구
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
        add_entry(inner, "Density cap n_max", "ncap", "1e22", "m⁻³")

        ttk.Separator(inner).pack(fill=tk.X, pady=3)
        ttk.Label(inner, text="Circuit / Run control (v18)", font=FONT_UI_B).pack(anchor="w")
        add_entry(inner, "Series resistance R", "rser", 0.0, "Ω")
        add_entry(inner, "Run until t_end (0=∞)", "tend", 0.0, "ns")
        add_entry(inner, "Redraw interval", "redraw", 2.0, "s")
        add_entry(inner, "Compute per tick", "tick", 0.25, "s")
        self.auto_ag = tk.BooleanVar(value=True)
        ttk.Checkbutton(inner, text="Auto afterglow when gap is bridged",
                        variable=self.auto_ag).pack(anchor="w")
        ttk.Label(inner, text="(R=0: ideal source as v17; R>0: V_needle = V_src - R·I, "
                  "voltage collapse after bridging)", foreground="#666",
                  font=("Times New Roman", 10), wraplength=320).pack(anchor="w")

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
                rog_edge_mm=self._read("rogr"),
                R_series_ohm=self._read("rser"),
                n_cap=self._read("ncap"))
            s = self.sim
            self._last_draw = 0.0
            self._wall_run = 0.0
            self._bridge_notified = False
            self._msg = ("Poisson: sparse direct (scipy)" if s.poisson_used == "direct"
                         else "Poisson: SOR fallback (install scipy for speed)")
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
        self.axW.plot(tt*sc, vv/1e3, "g-", lw=2.0, label="Source voltage $V_{src}(t)$")
        if s.plate_mode == "FLOAT_METAL":
            self.axW.axhline(s.V_plate/1e3, color="m", ls="--", lw=1.5,
                             label=f"Floating-plate potential {s.V_plate/1e3:.2f} kV")
        if len(s.circ_hist) >= 2:
            C = np.array(s.circ_hist)
            if s.R_series > 0:
                self.axW.plot(C[:, 0]*sc, C[:, 2]/1e3, "r--", lw=1.6,
                              label="Needle voltage $V_{needle}=V_{src}-R\\,I$")
            axI = getattr(self, "axWI", None)
            if axI is None:
                self.axWI = axI = self.axW.twinx()
            axI.clear()
            axI.plot(C[:, 0]*sc, C[:, 3], "b-", lw=1.2, alpha=0.8, label="Discharge current I (Sato)")
            axI.set_ylabel("Current [A]", color="b")
            h2, l2 = axI.get_legend_handles_labels()
        else:
            h2, l2 = [], []
        if s.t_bridge is not None:
            self.axW.axvline(s.t_bridge*sc, color="k", ls=":", lw=1.2,
                             label=f"gap bridged {s.t_bridge*1e9:.1f} ns")
        h1, l1 = self.axW.get_legend_handles_labels()
        self.axW.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=12)
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
        try:
            files = export_csv_files(self.sim, d)
            messagebox.showinfo("Saved",
                f"{len(files)} CSV files saved to: {d}\n"
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
        if not s.afterglow:
            dtinfo += f" [{s.dt_limiter}]"
        wall = self._wall_run + ((time.perf_counter() - self._run_t0) if self._run_t0 else 0.0)
        speed = (s.t * 1e9 / wall) if wall > 0 else 0.0
        draw_txt = f" (redraw {self._t_draw:.2f} s)" if self._t_draw > 0 else ""
        bridge_txt = (f" Gap bridged   : t = {s.t_bridge*1e9:.2f} ns\n" if s.bridged else "")
        circ_txt = (f" V_src / V_ndl : {s.V_src_now/1e3:8.2f} / {s.Va_now/1e3:8.2f} kV\n"
                    if s.R_series > 0 else "")
        msg_txt = (f" >> {self._msg}\n" if self._msg else "")
        txt = (msg_txt +
               f" Config        : {cfg_lbl}\n"
               f" Mode          : {mode}\n"
               f" Time t        : {self._fmt_t(s.t):>12s}\n"
               f" Wall / speed  : {wall:8.1f} s / {speed:.3f} ns per s{draw_txt}\n"
               f" Substeps      : {s.n_substeps}\n"
               f" dt applied    : {dtinfo}\n"
               f" V(t) applied  : {s.Va_now/1e3:10.2f} kV\n"
               + circ_txt +
               f" I discharge   : {s.I_dis:10.3e} A (Sato, R={s.R_series:g} Ω)\n"
               + bridge_txt +
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
            self._run_t0 = time.perf_counter()
            self._loop()
        else:
            self._pause_clock()

    def _pause_clock(self):
        if self._run_t0 is not None:
            self._wall_run += time.perf_counter() - self._run_t0
            self._run_t0 = None

    def _stop(self, msg):
        self.running = False
        self.btn_run.config(text="▶ Run")
        self._pause_clock()
        self._msg = msg

    def toggle_afterglow(self):
        if self.sim is None:
            return
        self.sim.afterglow = not self.sim.afterglow
        self.btn_ag.config(text="⚡ Back to discharge" if self.sim.afterglow
                           else "🌙 Afterglow")
        self._update_status()

    def _advance(self, nsub):
        """한 프레임 전진 + 브리징/정체 자동 처리"""
        s = self.sim
        if s.afterglow:
            s.step_afterglow(n_sub=nsub)
        else:
            s.step(n_sub=nsub)
        self.frame_count += 1
        if s.bridged and not self._bridge_notified:
            self._bridge_notified = True
            self._msg = (f"Gap bridged at {s.t_bridge*1e9:.2f} ns (streamer reached counter electrode)")
            if self.auto_ag.get() and not s.afterglow:
                s.afterglow = True
                self.btn_ag.config(text="⚡ Back to discharge")
                self._msg += " -> auto afterglow (field frozen)"
        if s.stalled and not s.afterglow:
            self._stop(f"Paused: dt collapsed below {s.dt_floor*1e12:g} ps at "
                       f"{s.t*1e9:.2f} ns (spark transition; use Afterglow or R > 0)")

    def step_once(self):
        if self.sim is None:
            return
        self._advance(self._read("nsub", int))
        self._refresh_plots()
        self._last_draw = time.perf_counter()
        self._update_status()

    def _loop(self):
        """벽시계 예산(tick)만큼 계산한 뒤, redraw 간격이 지났을 때만 그래프를 갱신.
        계산 예산은 최소 '직전 재렌더링 시간의 2배'로 자동 확대 -> 그리기가 느린 PC 에서도
        벽시계의 2/3 이상이 계산에 쓰이도록 보장 (v17: 6 서브스텝마다 전 탭 재렌더링)"""
        if not self.running:
            return
        budget = max(self._read("tick"), 0.02, 2.0 * self._t_draw)
        nsub = self._read("nsub", int)
        t_end = self._read("tend") * 1e-9
        t0 = time.perf_counter()
        while time.perf_counter() - t0 < budget:
            self._advance(nsub)
            if not self.running:
                break
            if t_end > 0 and self.sim.t >= t_end:
                self._stop(f"t_end = {t_end*1e9:g} ns reached")
                break
        if (time.perf_counter() - self._last_draw >= max(self._read("redraw"), 0.05)
                or not self.running):
            td = time.perf_counter()
            self._refresh_plots()
            self.root.update_idletasks()          # 캔버스 그리기를 지금 수행(시간 측정)
            self._t_draw = time.perf_counter() - td
            self._last_draw = time.perf_counter()
        self._update_status()
        if self.running:
            self.root.after(10, self._loop)       # 짧은 휴지: Tk 이벤트(버튼/스크롤) 처리 기회

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
    if not HAVE_SCIPY:
        messagebox.showwarning("scipy not found",
            "scipy is not installed: falling back to the slow v17 SOR Poisson solver.\n"
            "Install it for the fast direct solver:  pip install scipy")
    app = DischargeHMI(root)
    root.mainloop()
