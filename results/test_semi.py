"""반암시적(semi) vs 명시적(explicit) 비교: R=1 kOhm, 브리징 후 t_stop 까지"""
import sys, time, json; sys.path[:0]=["legacy","."]
import numpy as np, discharge_core as c, benchmark as B
t_stop=float(sys.argv[1])*1e-9; semi = sys.argv[2]=="semi"; out=sys.argv[3]
s=c.NeedlePlaneDischarge(**B.GUI_DEFAULTS, R_series_ohm=1000.0, semi_implicit=semi)
rows=[]; t0=time.perf_counter()
while s.t < t_stop and not s.stalled:
    s.step(6); rows.append((s.t, s.dt_last, s.Va_now, s.I_dis, s.n_e.max(), np.nanmax(s.Emag_masked()), s.head_z()))
    if s.n_substeps % 6000 < 6:
        print(f"[{sys.argv[2]}] t={s.t*1e9:7.3f} ns steps={s.n_substeps} semi={s.n_semi_steps} fact={s.n_semi_factor} dt={s.dt_last*1e12:.3f} ps [{s.dt_limiter}] Va={s.Va_now/1e3:.2f} kV I={s.I_dis:.2f} A ne={s.n_e.max():.2e} wall={time.perf_counter()-t0:.0f}s", flush=True)
W=time.perf_counter()-t0
print(f"[{sys.argv[2]}] DONE t={s.t*1e9:.3f} ns wall={W:.0f}s steps={s.n_substeps} semi={s.n_semi_steps} refactor={s.n_semi_factor} t_bridge={s.t_bridge*1e9 if s.t_bridge else None} ne_max={s.n_e.max():.3e} Va={s.Va_now:.1f} I={s.I_dis:.3f}")
np.savez(out, rows=np.array(rows), wall=W, steps=s.n_substeps, n_e=s.n_e, n_ip=s.n_ip, V=s.V, x=s.x, y=s.y, y_bot=s.y_bot)
