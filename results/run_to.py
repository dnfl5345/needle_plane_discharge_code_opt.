"""run_to.py <t_ns> <dt_rule> <out.npz> : GUI 기본값으로 t_ns 까지 실행, 이력+최종상태 저장"""
import sys, time; sys.path[:0]=["legacy","."]
import numpy as np, discharge_core as c18, benchmark as B
from run_batch import head_position
t_end=float(sys.argv[1])*1e-9; rule=sys.argv[2]; out=sys.argv[3]
s=c18.NeedlePlaneDischarge(**B.GUI_DEFAULTS, dt_rule=rule)
rows=[]; t0=time.perf_counter()
while s.t<t_end and not s.stalled and not (s.bridged and s.t>s.t_bridge+0.5e-9):
    s.step(6); rows.append((s.t, s.dt_last, s.n_e.max(), np.nanmax(s.Emag_masked()), head_position(s), s.Va_now))
W=time.perf_counter()-t0
E=s.Emag_masked(); j,i=np.unravel_index(np.nanargmax(E),E.shape)
print(f"{rule}: bridged={s.bridged} t_bridge={(s.t_bridge or 0)*1e9:.3f} ns | t={s.t*1e9:.3f} ns wall={W:.0f}s steps={s.n_substeps} dt={s.dt_last*1e12:.4f} ps ({s.dt_limiter}) ne_max={s.n_e.max():.2e} Emax={E[j,i]/1e5:.1f} kV/cm at (j={j},i={i}) z={(s.y[j]-s.y_bot)*1e3:.2f}mm x={s.x[i]*1e3:.2f}mm | ne,n+,n- there: {s.n_e[j,i]:.2e} {s.n_ip[j,i]:.2e} {s.n_in[j,i]:.2e}")
np.savez(out, rows=np.array(rows), n_e=s.n_e, n_ip=s.n_ip, n_in=s.n_in, V=s.V, Emag=s.Emag, x=s.x, y=s.y, y_bot=s.y_bot, t=s.t)
