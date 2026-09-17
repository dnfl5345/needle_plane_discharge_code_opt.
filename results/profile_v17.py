"""v17 원본 코어 헤드리스 프로파일: GUI 기본값으로 일정 wall-time 동안 step() 반복"""
import sys, time, cProfile, pstats, io
sys.path.insert(0, "legacy")
import numpy as np
import core_v17 as c

GUI_DEFAULTS = dict(gap_mm=10.0, voltage_kV=30.0, polarity=+1, pressure_atm=1.0,
                    r_main_mm=1.0, r_tip_mm=0.15, taper_mm=3.0, nx=161,
                    n_background=0.0, seed_density=0.0, wave_mode="IMPULSE",
                    t_front_us=0.01, t_tail_us=0.1, n2_pct=78, o2_pct=21, ar_pct=1, h2o_pct=0,
                    gamma_see=0.05, dt_user_ns=0.0, snap_ns=1.0, plate_mode="GND_METAL",
                    plate_thick_mm=1.0, plate_epsr=3.8, rate_model="TABLE",
                    cath_J0=1e-2, n_ion_bg=1e9, q_bg=1e7, config="NEEDLE_PLANE")

wall_budget = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
t0 = time.perf_counter(); s = c.NeedlePlaneDischarge(**GUI_DEFAULTS); t_init = time.perf_counter() - t0
print(f"init: {t_init:.1f} s  grid {s.nx}x{s.ny}")
dts, nemax, tsim, wall = [], [], [], []
t0 = time.perf_counter(); nframes = 0
while time.perf_counter() - t0 < wall_budget:
    s.step(n_sub=6); nframes += 1
    dts.append(s.dt_last); nemax.append(s.n_e.max()); tsim.append(s.t); wall.append(time.perf_counter()-t0)
    if nframes % 50 == 0:
        print(f"frame {nframes:5d} t={s.t*1e9:8.4f} ns dt={s.dt_last*1e12:7.3f} ps ne_max={s.n_e.max():.2e} "
              f"E_max={np.nanmax(s.Emag_masked())/1e5:7.1f} kV/cm  {wall[-1]/(6*nframes)*1e3:6.1f} ms/substep", flush=True)
W = time.perf_counter() - t0
print(f"\n{nframes} frames = {6*nframes} substeps in {W:.1f} s -> {W/(6*nframes)*1e3:.1f} ms/substep")
print(f"sim time reached {s.t*1e9:.4f} ns; mean dt {np.mean(dts)*1e12:.3f} ps; last dt {dts[-1]*1e12:.3f} ps")
rate = s.t / W
print(f"=> {rate*1e9:.5f} ns of simulation per wall-second; 300 ns impulse would need {300e-9/rate/3600:.1f} h (at current dt)")
np.savez("results/profile_v17.npz", dt=dts, nemax=nemax, tsim=tsim, wall=wall)

# 함수별 분해(짧게)
pr = cProfile.Profile(); pr.enable()
for _ in range(3): s.step(n_sub=6)
pr.disable(); sio = io.StringIO(); pstats.Stats(pr, stream=sio).sort_stats("cumulative").print_stats(12); print(sio.getvalue()[:3000])
