"""검증 2: v17 을 서브스텝당 SOR 3000회(수렴)로 돌리면 v18(직접해법) 궤적과 일치하는가"""
import sys, time, json; sys.path[:0]=["legacy","."]
import numpy as np, core_v17 as c17, discharge_core as c18, benchmark as B
P = B.GUI_DEFAULTS; t_end = float(sys.argv[1])*1e-9 if len(sys.argv)>1 else 3e-9
rows = []
s17 = c17.NeedlePlaneDischarge(**P); s18 = c18.NeedlePlaneDischarge(**P, dt_rule="v17"); s60 = c17.NeedlePlaneDischarge(**P)
t0 = time.perf_counter()
while s17.t < t_end:
    s17.step(n_sub=1, field_iters=3000); s18.step(n_sub=1); s60.step(n_sub=1)
    rows.append((s17.t*1e9, s17.dt_last*1e12, s18.dt_last*1e12, s17.n_e.max(), s18.n_e.max(), s60.n_e.max(),
                 np.nanmax(s17.Emag_masked())/1e5, np.nanmax(s18.Emag_masked())/1e5, np.nanmax(s60.Emag_masked())/1e5,
                 np.abs(s17.V-s18.V).max(), c18.NeedlePlaneDischarge.poisson_residual.__wrapped__(s17) if hasattr(c18.NeedlePlaneDischarge.poisson_residual,'__wrapped__') else 0.0))
    if len(rows) % 10 == 0:
        r = rows[-1]
        print(f"t={r[0]:6.3f} ns dt17={r[1]:7.2f} dt18={r[2]:7.2f} ps | ne_max: v17conv {r[3]:.3e} v18 {r[4]:.3e} v17(60it) {r[5]:.3e} | "
              f"Emax: {r[6]:6.2f} {r[7]:6.2f} {r[8]:6.2f} kV/cm | max|V17conv-V18|={r[9]:.2e} V  [{time.perf_counter()-t0:.0f}s]", flush=True)
A = np.array(rows); np.save("results/verify_converged.npy", A)
print(f"\nfinal t={A[-1,0]:.3f} ns: ne_max v17(converged SOR)={A[-1,3]:.4e}  v18={A[-1,4]:.4e}  ratio={A[-1,4]/A[-1,3]:.4f} | v17(60 it)={A[-1,5]:.4e} ratio={A[-1,5]/A[-1,3]:.3f}")
print(f"Emax v17conv={A[-1,6]:.2f} v18={A[-1,7]:.2f} v17(60it)={A[-1,8]:.2f} kV/cm ; max|V17conv - V18| over run = {A[:,9].max():.3e} V")
