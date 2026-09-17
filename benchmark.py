# -*- coding: utf-8 -*-
"""benchmark.py - v17 원본 코어 vs v18 최적화 코어: 속도 및 물리결과 비교
사용법: python benchmark.py [t_sim_ns=3.0] [nx=161]
  1) 동일 파라미터(GUI 기본값)로 t_sim 까지 세 가지를 실행
     A) v17 원본(legacy/core_v17.py)  B) v18, dt_rule="v17"  C) v18, dt_rule="v18"(기본)
  2) 서브스텝 수/벽시계/ms per substep, 최종 max ne, max |E|, 축상 프로파일 차이 출력
  3) results/benchmark_<nx>.json, results/benchmark_<nx>.png 저장
"""
import sys, os, json, time
sys.path[:0] = ["legacy", "."]
import numpy as np
import core_v17 as c17
import discharge_core as c18

GUI_DEFAULTS = dict(gap_mm=10.0, voltage_kV=30.0, polarity=+1, pressure_atm=1.0,
                    r_main_mm=1.0, r_tip_mm=0.15, taper_mm=3.0, nx=161,
                    n_background=0.0, seed_density=0.0, wave_mode="IMPULSE",
                    t_front_us=0.01, t_tail_us=0.1, n2_pct=78, o2_pct=21, ar_pct=1, h2o_pct=0,
                    gamma_see=0.05, dt_user_ns=0.0, snap_ns=1.0, plate_mode="GND_METAL",
                    plate_thick_mm=1.0, plate_epsr=3.8, rate_model="TABLE",
                    cath_J0=1e-2, n_ion_bg=1e9, q_bg=1e7, config="NEEDLE_PLANE")

def run(sim, t_end, label):
    t0 = time.perf_counter(); n0 = 0; frames = 0; hist = []
    while sim.t < t_end:
        sim.step(n_sub=6); frames += 1
        hist.append((sim.t, sim.dt_last, sim.n_e.max(), np.nanmax(sim.Emag_masked())))
    W = time.perf_counter() - t0
    nsub = 6 * frames
    ic = sim.nx // 2
    out = dict(label=label, wall_s=W, substeps=nsub, ms_per_substep=1e3 * W / nsub,
               t_ns=sim.t * 1e9, ne_max=float(sim.n_e.max()), nip_max=float(sim.n_ip.max()),
               E_max_kVcm=float(np.nanmax(sim.Emag_masked()) / 1e5),
               ne_axis=sim.n_e[:, ic].copy(), E_axis=np.nan_to_num(sim.Emag[:, ic]).copy(),
               y_mm=sim.y * 1e3, hist=np.array(hist),
               O_total=float(sim.rad_hist[-1][1]["O(P)"]))
    print(f"[{label:22s}] t={out['t_ns']:.3f} ns  substeps={nsub:6d}  wall={W:7.1f} s  "
          f"{out['ms_per_substep']:6.2f} ms/substep  ne_max={out['ne_max']:.3e}  E_max={out['E_max_kVcm']:.1f} kV/cm", flush=True)
    return out

if __name__ == "__main__":
    t_sim = float(sys.argv[1]) * 1e-9 if len(sys.argv) > 1 else 3.0e-9
    nx = int(sys.argv[2]) if len(sys.argv) > 2 else 161
    P = dict(GUI_DEFAULTS, nx=nx)
    os.makedirs("results", exist_ok=True)
    res = {}
    t0 = time.perf_counter(); s = c17.NeedlePlaneDischarge(**P); ti17 = time.perf_counter() - t0
    print(f"v17 init {ti17:.1f} s, grid {s.nx}x{s.ny}")
    res["v17"] = run(s, t_sim, "v17 original")
    t0 = time.perf_counter(); s = c18.NeedlePlaneDischarge(**P, dt_rule="v17"); ti18 = time.perf_counter() - t0
    print(f"v18 init {ti18:.2f} s (poisson={s.poisson_used})")
    res["v18_dt17"] = run(s, t_sim, "v18 (dt rule v17)")
    s = c18.NeedlePlaneDischarge(**P, dt_rule="v18")
    res["v18"] = run(s, t_sim, "v18 (dt rule v18)")
    base = res["v17"]
    summ = {}
    for k, r in res.items():
        dne = np.abs(np.log10(np.maximum(r["ne_axis"], 1e6)) - np.log10(np.maximum(base["ne_axis"], 1e6))).max()
        dE = np.abs(r["E_axis"] - base["E_axis"]).max() / max(base["E_axis"].max(), 1)
        summ[k] = dict(wall_s=r["wall_s"], substeps=r["substeps"], ms_per_substep=r["ms_per_substep"],
                       ne_max=r["ne_max"], E_max_kVcm=r["E_max_kVcm"], O_total=r["O_total"],
                       speedup_vs_v17=base["wall_s"] / r["wall_s"],
                       max_dlog10_ne_axis=float(dne), max_rel_dE_axis=float(dE))
    summ["init_s"] = dict(v17=ti17, v18=ti18)
    summ["params"] = dict(P, t_sim_ns=t_sim * 1e9)
    json.dump(summ, open(f"results/benchmark_nx{nx}.json", "w"), indent=2, default=float)
    print("\n==== summary ====")
    for k in ("v17", "v18_dt17", "v18"):
        r = summ[k]
        print(f"{k:9s} wall {r['wall_s']:7.1f} s | {r['substeps']:6d} substeps | {r['ms_per_substep']:6.2f} ms/substep | "
              f"speedup x{r['speedup_vs_v17']:5.1f} | ne_max {r['ne_max']:.3e} | E_max {r['E_max_kVcm']:6.1f} | "
              f"max|dlog10 ne|axis {r['max_dlog10_ne_axis']:.3f} | max rel dE axis {r['max_rel_dE_axis']:.2e}")
    # ---- 그림 ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 8.5))
    sty = {"v17": ("k", "-", "v17 original"), "v18_dt17": ("tab:blue", "--", "v18 (dt rule v17)"), "v18": ("tab:red", ":", "v18 (dt rule v18)")}
    for k, r in res.items():
        c, ls, lb = sty[k]; h = r["hist"]
        ax[0, 0].semilogy(r["y_mm"], np.maximum(r["ne_axis"], 1e6), color=c, ls=ls, lw=2, label=lb)
        ax[0, 1].plot(r["y_mm"], r["E_axis"] / 1e5, color=c, ls=ls, lw=2, label=lb)
        ax[1, 0].semilogy(h[:, 0] * 1e9, h[:, 2], color=c, ls=ls, lw=2, label=lb)
        ax[1, 1].semilogy(h[:, 0] * 1e9, h[:, 1] * 1e12, color=c, ls=ls, lw=2, label=lb)
    ax[0, 0].set(xlabel="z [mm]", ylabel="$n_e$ on axis [m$^{-3}$]", title=f"On-axis $n_e$ at t = {t_sim*1e9:.1f} ns", ylim=(1e8, 1e22))
    ax[0, 1].set(xlabel="z [mm]", ylabel="|E| on axis [kV/cm]", title=f"On-axis |E| at t = {t_sim*1e9:.1f} ns")
    ax[1, 0].set(xlabel="t [ns]", ylabel="max $n_e$ [m$^{-3}$]", title="Peak electron density history")
    ax[1, 1].set(xlabel="t [ns]", ylabel="dt [ps]", title="Time step history")
    for a in ax.flat: a.grid(alpha=0.3); a.legend()
    fig.suptitle(f"v17 vs v18 benchmark (nx={nx}, needle-plane, 30 kV impulse 0.01/0.1 us)  "
                 f"wall: v17 {res['v17']['wall_s']:.0f} s / v18 {res['v18']['wall_s']:.0f} s", fontsize=13)
    fig.tight_layout(); fig.savefig(f"results/benchmark_nx{nx}.png", dpi=120)
    print("saved results/benchmark_nx%d.png/.json" % nx)
