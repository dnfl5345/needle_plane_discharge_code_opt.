# -*- coding: utf-8 -*-
"""
run_batch.py - 방전 시뮬레이터 v18 헤드리스(GUI 없음) 배치 실행기
==============================================================================
GUI 와 동일한 물리 코어(discharge_core.NeedlePlaneDischarge)를 화면 갱신 없이
지정 시각(t_end)까지 실행하고, GUI 각 탭에 해당하는 그림(PNG)·CSV·소요시간 보고서를
저장한다. 기본값은 GUI 기본값(침-평판, 10 mm, +30 kV, 0.01/0.1 us 임펄스, nx=161).

사용 예)
  python run_batch.py                               # 300 ns(=3*T_tail) 전체 임펄스
  python run_batch.py --t_end_ns 50 --out results/quick
  python run_batch.py --nx 241 --afterglow_us 20    # 격자 세분 + 잔광 20 us
  python run_batch.py --config NEEDLE_ROGOWSKI --volt 40 --polarity -1

출력(out 폴더):
  progress.csv        벽시계 로그 (t, wall, dt, 제한요인, ne_max, E_max, 스트리머 선단 z ...)
  timing.json         총 소요시간/서브스텝/스텝당 ms/dt 통계
  fig1_field_ne.png   ① 전계 + ② 전자밀도 (방전 단계 종료 시점)
  fig2_axial.png      ④ 축상 프로파일 (E, ne, n+, rho) 시간중첩
  fig3_radicals.png   ⑥ 라디칼 2D + 축상 8종        (_afterglow.png: 잔광 종료 시점)
  fig4_paper.png      ⑦ (a)E/N (b)ne (c)발광 스트릭 (d)활성종 총개수 이력 (_afterglow.png)
  fig5_initiation.png ⑧ 개시율/개시확률/지연시간 분포
  fig6_timing.png     파형·벽시계·dt·ne_max·선단 위치 이력
  fig7_ne_montage.png 전자밀도 2D 시간 몽타주
  csv/                GUI 'Save CSV' 와 동일 형식의 CSV 묶음
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm

from discharge_core import (NeedlePlaneDischarge, CONFIG_MAP, PLATE_PRESET,
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
plt.rcParams["mathtext.fontset"] = "stix"
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams.update({"font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
                     "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10})
LSTYLES = ["-", "--", "-.", ":", (0, (5, 1)), (0, (3, 1, 1, 1)), (0, (1, 1)), (0, (5, 2, 1, 2))]
MARKERS = ["o", "s", "^", "v", "D", "P", "X", "*", "<", ">"]


# ----------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(description="HV discharge simulator v18 - headless batch run")
    a = p.add_argument
    a("--out", default=None, help="output folder (default results/run_<config>_nx<nx>_<stamp>)")
    a("--t_end_ns", type=float, default=None, help="discharge-phase end time [ns] (default 3*T_tail, DC: 100)")
    a("--afterglow_us", type=float, default=0.0, help="afterglow phase duration after t_end [us]")
    a("--nsub", type=int, default=6, help="substeps per frame (frame = history record)")
    a("--log_s", type=float, default=30.0, help="progress log interval [wall seconds]")
    a("--snap2d_ns", type=float, default=0.0, help="2D ne snapshot interval for montage [ns] (0=auto)")
    a("--wall_limit_h", type=float, default=24.0, help="abort discharge phase after this wall time [h]")
    a("--no_csv", action="store_true", help="skip CSV export")
    a("--dt_rule", default="v18", choices=["v17", "v18"])
    a("--R_series", type=float, default=0.0, help="external series resistance [Ohm] (0 = ideal source, v17)")
    a("--n_cap", type=float, default=1e22, help="charged-species density cap [m^-3]")
    a("--dt_floor_ps", type=float, default=0.01, help="dt below this for 20 substeps => stalled")
    a("--on_bridge", default="auto", choices=["auto", "stop", "continue"],
      help="after gap bridging: stop discharge phase (auto: stop if R_series==0, else continue)")
    a("--post_bridge_ns", type=float, default=0.0, help="extra time simulated after bridging when stopping [ns]")
    a("--poisson", default="direct", choices=["direct", "sor"])
    # ---- 물리/기하 (GUI 기본값) ----
    a("--config", default="NEEDLE_PLANE", choices=list(CONFIG_MAP.values()))
    a("--gap", type=float, default=10.0);      a("--rmain", type=float, default=1.0)
    a("--rtip", type=float, default=0.15);     a("--taper", type=float, default=3.0)
    a("--nsep", type=float, default=4.0);      a("--srad", type=float, default=5.0)
    a("--rogw", type=float, default=6.0);      a("--rogr", type=float, default=2.0)
    a("--plate", default="Metal (grounded)", choices=list(PLATE_PRESET))
    a("--pthick", type=float, default=1.0);    a("--pepsr", type=float, default=0.0)
    a("--press", type=float, default=1.0)
    a("--n2", type=float, default=78.0);       a("--o2", type=float, default=21.0)
    a("--ar", type=float, default=1.0);        a("--h2o", type=float, default=0.0)
    a("--volt", type=float, default=30.0);     a("--polarity", type=int, default=+1, choices=[1, -1])
    a("--wave", default="IMPULSE", choices=["IMPULSE", "DC"])
    a("--tfront", type=float, default=0.01);   a("--ttail", type=float, default=0.1)
    a("--rmodel", default="TABLE", choices=["TABLE", "TOWNSEND"])
    a("--gamma", type=float, default=0.05);    a("--j0", type=float, default=1e-2)
    a("--nionbg", type=float, default=1e9);    a("--qbg", type=float, default=1e7)
    a("--dtns", type=float, default=0.0);      a("--snap", type=float, default=1.0)
    a("--nx", type=int, default=161)
    a("--nbg", type=float, default=0.0);       a("--seed", type=float, default=0.0)
    a("--rates_csv", default=None, help="BOLSIG+ rate CSV (header EN_Td,...)")
    return p


def make_sim(args):
    mode, epsr, _ = PLATE_PRESET[args.plate]
    if epsr is None:
        epsr = args.pepsr
    ov = None
    if args.rates_csv:
        from discharge_core import parse_rate_csv
        ov, loaded = parse_rate_csv(args.rates_csv)
        print("rates replaced:", ", ".join(loaded))
    return NeedlePlaneDischarge(
        gap_mm=args.gap, voltage_kV=args.volt, polarity=args.polarity, pressure_atm=args.press,
        r_main_mm=args.rmain, r_tip_mm=args.rtip, taper_mm=args.taper, nx=args.nx,
        n_background=args.nbg, seed_density=args.seed, wave_mode=args.wave,
        t_front_us=args.tfront, t_tail_us=args.ttail,
        n2_pct=args.n2, o2_pct=args.o2, ar_pct=args.ar, h2o_pct=args.h2o,
        gamma_see=args.gamma, dt_user_ns=args.dtns, snap_ns=args.snap,
        plate_mode=mode, plate_thick_mm=args.pthick, plate_epsr=epsr,
        rate_model=args.rmodel, rate_override=ov,
        cath_J0=args.j0, n_ion_bg=args.nionbg, q_bg=args.qbg,
        config=args.config, needle_sep_mm=args.nsep, sphere_R_mm=args.srad,
        rog_halfw_mm=args.rogw, rog_edge_mm=args.rogr,
        poisson=args.poisson, dt_rule=args.dt_rule,
        R_series_ohm=args.R_series, n_cap=args.n_cap, dt_floor_ps=args.dt_floor_ps)


# ----------------------------------------------------------------------------
def head_position(s, thr=1e18):
    """축 근방 ne > thr 인 가장 낮은 z [mm] (스트리머 선단; 없으면 nan) - 코어의 head_z 사용"""
    return s.head_z(thr) * 1e3


def fmt_hms(sec):
    sec = int(round(sec))
    return f"{sec//3600:d}h {(sec%3600)//60:02d}m {sec%60:02d}s"


# ----------------------------------------------------------------------------
class Recorder:
    MAX_SNAPS = 48                              # 2D 스냅샷 상한 (초과 시 절반 솎아내고 간격 2배)

    def __init__(self, s, out, log_s, snap2d_ns, t_end):
        self.s, self.out, self.log_s = s, out, log_s
        self.t0 = time.perf_counter()
        self.next_log = 0.0
        self.row_s = max(min(log_s / 30.0, 2.0), 0.2)   # 이력 행 기록 간격(벽시계) - 그림용
        self.next_row = 0.0
        self.snap2d_dt = snap2d_ns * 1e-9        # 초기 간격 (자동 확장)
        self.next_snap2d = 0.0
        self.snaps2d = []                       # (t, ne float32, Emag float32)
        self.rows = []
        self.t_end = t_end
        self.f = open(os.path.join(out, "progress.csv"), "w")
        self.f.write("t_ns,wall_s,substeps,dt_ps,limiter,Va_kV,ne_max,nip_max,E_max_kVcm,"
                     "z_head_mm,N_OP,N_O3,N_NO,F_init\n")

    def __call__(self, s):
        w = time.perf_counter() - self.t0
        if s.t >= self.next_snap2d:
            self.snaps2d.append((s.t, s.n_e.astype(np.float32), s.Emag.astype(np.float32)))
            self.next_snap2d = s.t + self.snap2d_dt
            if len(self.snaps2d) > self.MAX_SNAPS:
                self.snaps2d = self.snaps2d[::2]
                self.snap2d_dt *= 2.0
        if w >= self.next_log or s.t >= self.t_end:
            self.next_log = w + self.log_s
            self.next_row = w + self.row_s
            self.record(s, w, print_line=True)
        elif w >= self.next_row:
            self.next_row = w + self.row_s
            self.record(s, w, print_line=False)

    def record(self, s, w, print_line=False):
        Emax = float(np.nanmax(s.Emag_masked())) / 1e5
        tot = s.rad_hist[-1][1]
        F = s.init_stats[-1][3] if s.init_stats else 0.0
        zh = head_position(s)
        row = (s.t * 1e9, w, s.n_substeps, s.dt_last * 1e12, s.dt_limiter, s.Va_now / 1e3,
               float(s.n_e.max()), float(s.n_ip.max()), Emax, zh,
               tot["O(P)"], tot["O3"], tot["NO"], F)
        self.rows.append(row)
        if print_line:
            self.f.write(",".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in row) + "\n")
            self.f.flush()
            frac = min(s.t / self.t_end, 1.0) if self.t_end > 0 else 1.0
            eta = (w / frac - w) if frac > 0.002 else float("nan")
            print(f"  t={s.t*1e9:9.3f} ns ({frac*100:5.1f}%)  wall {fmt_hms(w)}  ETA {fmt_hms(eta) if eta==eta else '--'}"
                  f"  steps={s.n_substeps:7d}  dt={s.dt_last*1e12:7.3f} ps [{s.dt_limiter}]"
                  f"  V={s.Va_now/1e3:6.2f} kV  ne_max={s.n_e.max():.2e}  E_max={Emax:6.1f} kV/cm"
                  f"  head z={zh:5.2f} mm", flush=True)

    def close(self):
        self.f.close()


# ----------------------------------------------------------------------------
#  그림
# ----------------------------------------------------------------------------
def _draw_solids(ax, s, plate_color):
    ax.contourf(s.X*1e3, s.Y*1e3, s.needle.astype(float), levels=[0.5, 1.5], colors=["#909090"])
    if s.gnd.any():
        ax.contourf(s.X*1e3, s.Y*1e3, s.gnd.astype(float), levels=[0.5, 1.5], colors=["#6a6a6a"])
    if s.n_solid >= 1:
        ax.axhspan(0, s.t_solid*1e3, color=plate_color, alpha=0.85)


def _tstr(t):
    return f"{t*1e9:.2f} ns" if t < 1e-6 else f"{t*1e6:.3f} µs"


def fig_field_ne(s, out, plate_color, title):
    ext = [s.x[0]*1e3, s.x[-1]*1e3, s.y[0]*1e3, s.y[-1]*1e3]
    fig, ax = plt.subplots(1, 2, figsize=(13, 6.2))
    Ekv = s.Emag_masked() / 1e5
    im = ax[0].imshow(Ekv, origin="lower", extent=ext, cmap="inferno", vmin=0,
                      vmax=max(np.nanmax(Ekv)*1.02, 1e-3), aspect="equal")
    fig.colorbar(im, ax=ax[0], fraction=0.046, label="|E| [kV/cm]")
    Vk = s.V.copy(); Vk[s.needle] = np.nan
    ax[0].contour(s.X*1e3, s.Y*1e3, Vk/1e3, 12, colors="w", linewidths=0.5, alpha=0.65)
    _draw_solids(ax[0], s, plate_color)
    ax[0].set_title("Electric field |E| and equipotential lines")
    ax[0].set_xlabel("x [mm]"); ax[0].set_ylabel("z [mm]  (0 = bottom electrode)")
    im2 = ax[1].imshow(np.maximum(s.n_e, 1e6), origin="lower", extent=ext, cmap="viridis",
                       norm=LogNorm(vmin=max(s.n_bg, 1e8), vmax=1e21), aspect="equal")
    fig.colorbar(im2, ax=ax[1], fraction=0.046, label="$n_e$ [m$^{-3}$] (log)")
    _draw_solids(ax[1], s, plate_color)
    ax[1].set_title(f"Electron density $n_e$ (max {s.n_e.max():.2e} m$^{{-3}}$)")
    ax[1].set_xlabel("x [mm]"); ax[1].set_ylabel("z [mm]  (0 = bottom electrode)")
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig1_field_ne.png"), dpi=130); plt.close(fig)


def fig_axial(s, out, plate_color, max_overlay=10):
    yy = s.y * 1e3
    xmax = min((s.t_solid + 1.2 * s.gap) * 1e3, s.height * 1e3)
    snaps = s.profiles
    if len(snaps) > max_overlay:
        pick = np.unique(np.linspace(0, len(snaps) - 1, max_overlay).astype(int))
        snaps = [snaps[i] for i in pick]
    colors = cm.turbo(np.linspace(0.1, 0.95, len(snaps)))
    fig, ax = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
    for i, ((tsn, E_ax, ne_ax, nip_ax, rho_ax), c) in enumerate(zip(snaps, colors)):
        kw = dict(color=c, lw=1.5, ls=LSTYLES[i % len(LSTYLES)], marker=MARKERS[i % len(MARKERS)],
                  markevery=max(len(yy)//9, 1), ms=5, mfc="none", label=f"t={tsn*1e9:.1f} ns")
        ax[0].plot(yy, E_ax/1e5, **kw)
        ax[1].semilogy(yy, np.maximum(ne_ax, 1e6), **kw)
        ax[2].semilogy(yy, np.maximum(nip_ax, 1e6), **kw)
        ax[3].plot(yy, rho_ax, **kw)
    ax[0].axhline(30*s.p_atm, color="k", ls=":", lw=1.0)
    ax[0].set_ylabel("|E| [kV/cm]"); ax[0].set_title("On-axis electric field (snapshots)")
    ax[0].legend(loc="upper left", ncol=2, framealpha=0.6)
    ax[1].set_ylabel("$n_e$ [m$^{-3}$]"); ax[1].set_title("On-axis electron density"); ax[1].set_ylim(1e8, 1e22)
    ax[2].set_ylabel("$n_+$ [m$^{-3}$]"); ax[2].set_title("On-axis positive-ion density"); ax[2].set_ylim(1e8, 1e22)
    ax[3].set_ylabel(r"$\rho$ [C/m$^3$]"); ax[3].set_yscale("symlog", linthresh=1e-3)
    ax[3].axhline(0, color="k", lw=0.7)
    ax[3].set_title(r"On-axis net charge density  $\rho = q(n_+ - n_e - n_-)$")
    ax[3].set_xlabel("Axial position z [mm]  (0 = plate, needle tip at z = gap)")
    for a in ax:
        a.set_xlim(0, xmax); a.grid(alpha=0.3); a.axvline(s.y_tip*1e3, color="#888", ls="--", lw=1.0)
        if s.n_solid >= 1:
            a.axvspan(0, s.t_solid*1e3, color=plate_color, alpha=0.4)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig2_axial.png"), dpi=130); plt.close(fig)


def fig_radicals(s, out, plate_color, suffix=""):
    ext = [s.x[0]*1e3, s.x[-1]*1e3, s.y[0]*1e3, s.y[-1]*1e3]
    fig = plt.figure(figsize=(11, 10))
    gr = fig.add_gridspec(2, 1, height_ratios=[1.5, 1.0], hspace=0.35)
    ax1 = fig.add_subplot(gr[0]); ax2 = fig.add_subplot(gr[1])
    im = ax1.imshow(np.maximum(s.rad["O(P)"], 1e6), origin="lower", extent=ext, cmap="cividis",
                    norm=LogNorm(vmin=1e14, vmax=1e23), aspect="equal")
    fig.colorbar(im, ax=ax1, fraction=0.046, label="Radical density [m$^{-3}$] (log)")
    _draw_solids(ax1, s, plate_color)
    ax1.set_title(f"2D distribution of O(P)  (max {s.rad['O(P)'].max():.2e} m$^{{-3}}$), t={_tstr(s.t)}")
    ax1.set_xlabel("x [mm]"); ax1.set_ylabel("z [mm]  (0 = bottom electrode)")
    ic = s.nx // 2; yy = s.y * 1e3
    xmax = min((s.t_solid + 1.2 * s.gap) * 1e3, s.height * 1e3)
    for name in RAD_SPECIES:
        c_, ls_, mk_ = RAD_STYLE[name]
        ax2.semilogy(yy, np.maximum(s.rad[name][:, ic], 1e6), color=c_, ls=ls_, marker=mk_,
                     markevery=max(len(yy)//8, 1), ms=5, mfc="none", lw=1.6, label=name)
    ax2.set_xlim(0, xmax); ax2.set_ylim(1e12, 1e24)
    ax2.set_xlabel("Axial position z [mm]  (0 = bottom electrode)"); ax2.set_ylabel("Density [m$^{-3}$]")
    ax2.set_title("On-axis radical density profiles (final time)")
    ax2.grid(alpha=0.3, which="both"); ax2.legend(loc="upper left", ncol=4)
    ax2.axvline(s.y_tip*1e3, color="#888", ls="--", lw=1.0)
    fig.savefig(os.path.join(out, f"fig3_radicals{suffix}.png"), dpi=130); plt.close(fig)


def fig_paper(s, out, suffix=""):
    fig = plt.figure(figsize=(14, 10.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1.0], wspace=0.55, hspace=0.40)
    ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2]); ax4 = fig.add_subplot(gs[1, :])
    j0 = s.j_ax_lo if (s.gnd.any() or s.n_solid >= 1) else 0
    j1 = max(s.j_ax_hi, j0 + 1)
    gapmm = s.gap * 1e3; hm = s.h * 1e3
    z_hi = (s.y[j1] - s.y_bot) * 1e3
    ext = [s.x[0]*1e3, s.x[-1]*1e3, -0.5*hm, z_hi + 0.5*hm]
    xr = min(0.20 * s.gap * 1e3, s.x[-1]*1e3)
    sub = slice(j0, j1 + 1)
    EN = (s.Emag / (s.Ngas * TD))[sub]
    EN = np.where(s.needle[sub] | s.gnd[sub], np.nan, EN)
    im1 = ax1.imshow(EN, origin="lower", extent=ext, aspect="auto", cmap="jet", vmin=0,
                     vmax=max(np.nanmax(EN)*1.02, 10.0))
    fig.colorbar(im1, ax=ax1, fraction=0.10)
    ax1.set_title("(a) E/N [Td]"); ax1.set_xlim(-xr, xr); ax1.set_ylim(0.0, gapmm)
    ax1.set_xlabel("r [mm]"); ax1.set_ylabel("z [mm]  (0 = bottom electrode)")
    nec = np.maximum(s.n_e[sub] * 1e-6, 1e6)
    im2 = ax2.imshow(nec, origin="lower", extent=ext, aspect="auto", cmap="jet", norm=LogNorm(vmin=1e9, vmax=1e16))
    fig.colorbar(im2, ax=ax2, fraction=0.10)
    ax2.set_title("(b) $n_e$ [cm$^{-3}$]"); ax2.set_xlim(-xr, xr); ax2.set_ylim(0.0, gapmm)
    ax2.set_xlabel("r [mm]"); ax2.set_ylabel("z [mm]  (0 = bottom electrode)")
    if len(s.streak_em) >= 3:
        tt = np.array([p[0] for p in s.streak_em]) * 1e9
        M = np.array([np.maximum(p[1][sub], 1e-30) for p in s.streak_em])
        M = np.log10(np.maximum(M.T, 1e-30))
        vhi = M.max()
        zed = (s.y[j0:j1 + 1] - s.y_bot) * 1e3
        pm = ax3.pcolormesh(*np.meshgrid(tt, zed), M, cmap="hot", shading="nearest", vmin=vhi - 7.0, vmax=vhi)
        fig.colorbar(pm, ax=ax3, fraction=0.10, label="log$_{10}$ emission [a.u.]")
        ax3.set_ylim(0.0, gapmm)
    ax3.set_title("(c) Optical emission streak (Simulation)")
    ax3.set_xlabel("Time [ns]"); ax3.set_ylabel("z [mm]  (0 = bottom electrode)")
    if len(s.rad_hist) >= 3:
        tt = np.array([p[0] for p in s.rad_hist]) * 1e9
        m = tt > 0
        for sp in RAD_SPECIES:
            vals = np.array([p[1][sp] for p in s.rad_hist])
            c_, ls_, mk_ = RAD_STYLE[sp]
            npt = int(m.sum())
            ax4.plot(tt[m], np.maximum(vals[m], 1e6), color=c_, ls=ls_, marker=mk_,
                     markevery=max(npt//10, 1), ms=5, mfc="none", lw=1.6, label=sp)
        ax4.set_xscale("log"); ax4.set_yscale("log")
        ax4.set_xlim(max(tt[m][0], 1.0), max(tt[-1]*1.1, 10.0))
        vmax_all = max(max(p[1][sp] for sp in RAD_SPECIES) for p in s.rad_hist)
        ax4.set_ylim(max(vmax_all*1e-5, 1e8), vmax_all*3)
        ax4.legend(loc="lower right", ncol=8, framealpha=0.6)
    gases = (f"H₂O({s.nH2O/s.Ngas*100:.0f}%)/" if s.nH2O > 0 else "") + "O₂/N₂"
    ax4.set_title(f"(d) Total number of radicals (axisymmetric estimate)  [{gases}]")
    ax4.set_xlabel("Time [ns]"); ax4.set_ylabel("Total number of radicals"); ax4.grid(alpha=0.3, which="both")
    fig.suptitle(f"Streak / history  (t = {_tstr(s.t)}, V(t) = {s.Va_now/1e3:+.2f} kV)", fontsize=13)
    fig.savefig(os.path.join(out, f"fig4_paper{suffix}.png"), dpi=130); plt.close(fig)


def fig_initstats(s, out):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 9))
    if len(s.init_stats) >= 3:
        tt = np.array([p[0] for p in s.init_stats]) * 1e9
        Lam = np.array([p[1] for p in s.init_stats]); F = np.array([p[3] for p in s.init_stats])
        ax1.semilogy(tt, np.maximum(Lam, 1e-3), "b-", lw=1.8, label=r"$\Lambda(t)$ [1/s]")
        ax1.set_ylabel(r"$\Lambda(t)$ [1/s]", color="b")
        axb = ax1.twinx(); axb.plot(tt, F, "r-", lw=1.8); axb.set_ylabel("Initiation probability F(t)", color="r")
        axb.set_ylim(0, 1.02)
        ax1.set_title(r"Effective initiation rate $\Lambda=\int_{\alpha_{eff}>0}\nu_{det}\,n_{O_2^-}\,dV$"
                      r"  and  $F(t)=1-e^{-\int\Lambda dt}$")
    ax1.set_xlabel("Time [ns]"); ax1.grid(alpha=0.3, which="both")
    samples, cens = s.sample_delays(2000, rng=np.random.default_rng(0))
    if samples.size >= 5:
        ax2.hist(samples * 1e9, bins=40, color="#4477aa", edgecolor="k", linewidth=0.4)
        med = np.median(samples) * 1e9
        ax2.axvline(med, color="r", ls="--", lw=1.5)
        ax2.set_title(f"Statistical delay-time distribution (N=2000, median={med:.2f} ns, "
                      f"censored beyond t_end: {cens*100:.1f}%)")
    else:
        Fe = 0 if not s.init_stats else s.init_stats[-1][3]*100
        ax2.set_title(f"Statistical delay-time distribution (insufficient hazard; F(t_end)={Fe:.2f}%)")
    ax2.set_xlabel("Delay time [ns]"); ax2.set_ylabel("Counts"); ax2.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig5_initiation.png"), dpi=130); plt.close(fig)


def fig_timing(s, rec, out, title):
    R = rec.rows
    t = np.array([r[0] for r in R]); w = np.array([r[1] for r in R])
    ne = np.array([r[6] for r in R]); Em = np.array([r[8] for r in R]); zh = np.array([r[9] for r in R])
    D = np.array([(a, b) for a, b, _ in s.dt_hist]); lim = [c for _, _, c in s.dt_hist]
    fig, ax = plt.subplots(2, 3, figsize=(17, 9))
    # (a) 회로: 전원/전극 전압, 전류
    tt, vv = s.waveform_curve(800)
    ax[0, 0].plot(tt*1e9, vv/1e3, "g-", lw=2, label="Source voltage $V_{src}(t)$")
    h2, l2 = [], []
    if s.circ_hist:
        C = np.array(s.circ_hist)
        if s.R_series > 0:
            ax[0, 0].plot(C[:, 0]*1e9, C[:, 2]/1e3, "r--", lw=1.5, label="Needle voltage $V_{needle}=V_{src}-RI$")
        axc = ax[0, 0].twinx()
        axc.plot(C[:, 0]*1e9, C[:, 3], "b-", lw=1.2, alpha=0.8, label="Discharge current I (Sato)")
        axc.set_ylabel("Current [A]", color="b")
        h2, l2 = axc.get_legend_handles_labels()
    if s.t_bridge is not None:
        ax[0, 0].axvline(s.t_bridge*1e9, color="k", ls=":", lw=1.2, label=f"gap bridged {s.t_bridge*1e9:.2f} ns")
    ax[0, 0].set_xlim(0, max(s.t*1e9*1.05, 1.0))
    h1, l1 = ax[0, 0].get_legend_handles_labels()
    ax[0, 0].legend(h1 + h2, l1 + l2, loc="best", fontsize=9)
    ax[0, 0].set(xlabel="Time [ns]", ylabel="Voltage [kV]", title=f"(a) Circuit: R_series = {s.R_series:g} Ω")
    # (b) 벽시계
    ax[0, 1].plot(t, w/60.0, "k-", lw=2)
    ax[0, 1].set(xlabel="Simulated time [ns]", ylabel="Wall-clock time [min]",
                 title=f"(b) Wall time vs simulated time ({fmt_hms(w[-1])}, {s.n_substeps} substeps)")
    # (c) dt / 제한요인
    cols = {"CFL": "tab:blue", "ionization": "tab:red", "dielectric": "tab:green", "waveform": "tab:orange",
            "diffusion": "tab:purple", "manual": "k"}
    for name in sorted(set(lim)):
        m = np.array([l == name for l in lim])
        base = name.replace("semi:", "")
        ax[0, 2].semilogy(D[m, 0]*1e9, D[m, 1]*1e12, "." if not name.startswith("semi") else "x", ms=4,
                          color=cols.get(base, "gray"), label=f"{name} ({m.mean()*100:.0f}%)")
    ax[0, 2].set(xlabel="Simulated time [ns]", ylabel="dt [ps]", title="(c) Time step and active limiter (x: semi-implicit)")
    ax[0, 2].legend(fontsize=9)
    # (d) ne_max, (e) E_max, (f) head z
    ax[1, 0].semilogy(t, np.maximum(ne, 1e6), "b-", lw=1.8)
    ax[1, 0].set(xlabel="Simulated time [ns]", ylabel="max $n_e$ [m$^{-3}$]", title="(d) Peak electron density")
    ax[1, 1].plot(t, Em, "r-", lw=1.8)
    ax[1, 1].set(xlabel="Simulated time [ns]", ylabel="max |E| [kV/cm]", title="(e) Peak field (gas region)")
    ok = np.isfinite(zh)
    ax[1, 2].plot(t[ok], zh[ok], "g-", lw=1.8, label="head z ($n_e$ > 1e18 m$^{-3}$)")
    hv, lv = [], []
    if ok.sum() >= 3:
        v = -np.gradient(zh[ok], t[ok])          # mm/ns
        axv = ax[1, 2].twinx(); axv.plot(t[ok], v, "m:", lw=1.3, label="velocity [mm/ns]")
        axv.set_ylabel("Streamer velocity [mm/ns]", color="m"); axv.set_ylim(0, max(np.nanmax(v)*1.2, 0.1))
        hv, lv = axv.get_legend_handles_labels()
    ax[1, 2].set_ylim(0, s.gap*1e3*1.05)
    hh, ll = ax[1, 2].get_legend_handles_labels(); ax[1, 2].legend(hh + hv, ll + lv, fontsize=9)
    ax[1, 2].set(xlabel="Simulated time [ns]", ylabel="Head position z [mm]", title="(f) Streamer head position / velocity")
    for a in ax.flat: a.grid(alpha=0.3)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(out, "fig6_timing.png"), dpi=120); plt.close(fig)


def fig_montage(s, rec, out, plate_color, n=12):
    snaps = rec.snaps2d
    if len(snaps) < 2:
        return
    pick = np.unique(np.linspace(0, len(snaps) - 1, n).astype(int))
    snaps = [snaps[i] for i in pick]
    ext = [s.x[0]*1e3, s.x[-1]*1e3, s.y[0]*1e3, s.y[-1]*1e3]
    ncol = 4; nrow = int(np.ceil(len(snaps) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.0*ncol, 4.6*nrow), squeeze=False)
    xr = min(0.35 * s.gap * 1e3, s.x[-1]*1e3)
    for a, (t, ne, Em) in zip(axs.flat, snaps):
        im = a.imshow(np.maximum(ne, 1e6), origin="lower", extent=ext, cmap="viridis",
                      norm=LogNorm(vmin=1e12, vmax=1e21), aspect="equal")
        _draw_solids(a, s, plate_color)
        a.set_xlim(-xr, xr); a.set_ylim(0, (s.y_tip + 0.15*s.gap)*1e3)
        a.set_title(f"t = {t*1e9:.1f} ns, max $n_e$={ne.max():.1e}", fontsize=10)
        a.set_xlabel("x [mm]"); a.set_ylabel("z [mm]")
    for a in list(axs.flat)[len(snaps):]:
        a.axis("off")
    fig.colorbar(im, ax=axs, fraction=0.02, label="$n_e$ [m$^{-3}$] (log)")
    fig.suptitle("Electron density $n_e$ - time montage", fontsize=13)
    fig.savefig(os.path.join(out, "fig7_ne_montage.png"), dpi=120); plt.close(fig)


# ----------------------------------------------------------------------------
def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.t_end_ns is None:
        args.t_end_ns = 3.0 * max(args.ttail, 2*args.tfront) * 1e3 if args.wave == "IMPULSE" else 100.0
    stamp = time.strftime("%Y%m%d_%H%M%S")
    out = args.out or os.path.join("results", f"run_{args.config}_nx{args.nx}_{stamp}")
    os.makedirs(out, exist_ok=True)
    print(f"== HV discharge simulator v18 batch ==  out: {out}")
    print(f"   scipy: {'yes' if HAVE_SCIPY else 'NO (SOR fallback, slow)'}  poisson={args.poisson} dt_rule={args.dt_rule}")
    t0 = time.perf_counter()
    s = make_sim(args)
    t_init = time.perf_counter() - t0
    if s.rtip < s.h:
        print(f"   [warn] tip radius {s.rtip*1e3:.3f} mm < grid {s.h*1e3:.3f} mm: tip = one cell (increase nx)")
    print(f"   grid {s.nx}x{s.ny} (h={s.h*1e3:.3f} mm), unknowns={getattr(s,'n_unknowns','-')}, "
          f"init {t_init:.2f} s (LU factor {getattr(s,'t_wall_factor',0)*1e3:.0f} ms), poisson={s.poisson_used}")
    t_end = args.t_end_ns * 1e-9
    snap2d = args.snap2d_ns if args.snap2d_ns > 0 else 0.25       # 자동: 0.25 ns 로 시작, 48장 초과 시 간격 2배
    rec = Recorder(s, out, args.log_s, snap2d, t_end)
    print(f"-- discharge phase: 0 -> {args.t_end_ns:g} ns  (V(t): {args.wave}, {args.tfront:g}/{args.ttail:g} us, "
          f"Vp={args.polarity*args.volt:+g} kV, config={args.config})")
    rec.record(s, 0.0, print_line=True)
    stop_on_bridge = (args.on_bridge == "stop") or (args.on_bridge == "auto" and args.R_series <= 0)
    stat = s.run_until(t_end, n_sub=args.nsub, on_frame=rec, wall_limit=args.wall_limit_h*3600,
                       stop_on_bridge=stop_on_bridge, post_bridge=args.post_bridge_ns*1e-9)
    W_dis = stat["wall_s"]
    rec.record(s, W_dis, print_line=True)
    if s.bridged:
        print(f"   ** gap bridged (streamer reached counter electrode) at t = {s.t_bridge*1e9:.3f} ns **")
    if stat["stop_reason"] == "bridged":
        print(f"   discharge phase stopped {args.post_bridge_ns:g} ns after bridging (fluid-model validity ends: "
              f"streamer-to-spark transition; use --R_series > 0 or --on_bridge continue to go on)")
    elif stat["stop_reason"] == "dt_stalled":
        print(f"   [stop] dt collapsed below {args.dt_floor_ps} ps at t={s.t*1e9:.3f} ns (numerical stall guard)")
    elif stat["stop_reason"] == "wall_limit":
        print(f"   [abort] wall limit {args.wall_limit_h} h reached at t={s.t*1e9:.2f} ns")
    _, _, plate_color = PLATE_PRESET[args.plate]
    title = (f"{args.config}, gap {args.gap:g} mm, Vp={args.polarity*args.volt:+g} kV {args.wave}, "
             f"nx={s.nx}, t={_tstr(s.t)}, wall {fmt_hms(W_dis)}")
    tf = time.perf_counter()
    fig_field_ne(s, out, plate_color, title)
    fig_axial(s, out, plate_color)
    fig_radicals(s, out, plate_color, suffix="")
    fig_paper(s, out, suffix="")
    fig_initstats(s, out)
    fig_timing(s, rec, out, title)
    fig_montage(s, rec, out, plate_color)
    print(f"-- discharge-phase figures written in {time.perf_counter()-tf:.1f} s")
    W_ag = 0.0
    if args.afterglow_us > 0 and stat["stop_reason"] in ("t_end", "bridged"):
        print(f"-- afterglow phase: {args.afterglow_us:g} us (field frozen, radical diffusion + chemistry)")
        s.afterglow = True
        t_ag_end = s.t + args.afterglow_us * 1e-6
        ta = time.perf_counter(); nag = 0
        while s.t < t_ag_end:
            s.step_afterglow(n_sub=4); nag += 4
        W_ag = time.perf_counter() - ta
        print(f"   afterglow done: {nag} substeps, {fmt_hms(W_ag)}, t={s.t*1e6:.3f} us")
        fig_radicals(s, out, plate_color, suffix="_afterglow")
        fig_paper(s, out, suffix="_afterglow")
    rec.close()

    # ---- 결과 통계 ----
    dts = np.array([d[1] for d in s.dt_hist]); lims = [d[2] for d in s.dt_hist]
    lim_share = {k: float(np.mean([l == k for l in lims])) for k in sorted(set(lims))}
    W_tot = time.perf_counter() - t0
    timing = dict(
        out=out, config=args.config, nx=s.nx, ny=s.ny, h_mm=s.h*1e3, unknowns=getattr(s, "n_unknowns", None),
        poisson=s.poisson_used, dt_rule=s.dt_rule, t_end_ns=args.t_end_ns, t_reached_ns=s.t*1e9,
        t_end_reached=stat["t_end_reached"], stop_reason=stat["stop_reason"], afterglow_us=args.afterglow_us,
        t_bridge_ns=(s.t_bridge*1e9 if s.t_bridge is not None else None), R_series=args.R_series,
        I_max_A=float(max(abs(r[3]) for r in s.circ_hist)) if s.circ_hist else None,
        wall_init_s=t_init, wall_discharge_s=W_dis, wall_afterglow_s=W_ag, wall_total_s=W_tot,
        wall_total_hms=fmt_hms(W_tot), substeps=s.n_substeps, frames=stat["frames"],
        ms_per_substep=1e3*W_dis/max(s.n_substeps, 1), field_solve_share=s.t_wall_field/max(s.t_wall_step, 1e-9),
        dt_ps_mean=float(np.mean(dts))*1e12 if dts.size else None, dt_ps_min=float(dts.min())*1e12 if dts.size else None,
        dt_ps_max=float(dts.max())*1e12 if dts.size else None, dt_limiter_share=lim_share,
        ne_max=float(s.n_e.max()), E_max_kVcm=float(np.nanmax(s.Emag_masked()))/1e5,
        head_z_mm=head_position(s), F_init=float(s.init_stats[-1][3]) if s.init_stats else None,
        M_townsend=s.M_townsend, gamma_MK=s.gamma_MK,
        radicals_total={sp: s.rad_hist[-1][1][sp] for sp in RAD_SPECIES},
        params=vars(args), cpu_count=os.cpu_count(), python=sys.version.split()[0])
    json.dump(timing, open(os.path.join(out, "timing.json"), "w"), indent=2, default=float)

    # ---- CSV ----
    if not args.no_csv:
        tf = time.perf_counter()
        export_csv_files(s, os.path.join(out, "csv"), stamp)
        print(f"-- CSV written in {time.perf_counter()-tf:.1f} s")
    print("\n==== SUMMARY ====")
    print(f" simulated {s.t*1e9:.3f} ns  |  substeps {s.n_substeps}  |  {timing['ms_per_substep']:.2f} ms/substep  "
          f"|  field solve {timing['field_solve_share']*100:.0f}% of step time")
    print(f" wall: init {t_init:.2f} s, discharge {fmt_hms(W_dis)}, afterglow {fmt_hms(W_ag)}, TOTAL {fmt_hms(W_tot)}")
    print(f" dt [ps]: mean {timing['dt_ps_mean']:.3f} min {timing['dt_ps_min']:.3f} max {timing['dt_ps_max']:.3f} ; limiter share "
          + ", ".join(f"{k} {v*100:.0f}%" for k, v in lim_share.items()))
    print(f" final: max ne {timing['ne_max']:.3e} m^-3, max |E| {timing['E_max_kVcm']:.1f} kV/cm, head z {timing['head_z_mm']:.2f} mm, "
          f"F_init {timing['F_init']}")
    print(" radicals (total number):", ", ".join(f"{sp} {v:.2e}" for sp, v in timing["radicals_total"].items()))
    return timing


if __name__ == "__main__":
    main()
