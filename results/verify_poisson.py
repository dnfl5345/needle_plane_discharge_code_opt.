"""검증 1: v18 직접해법 vs v17 SOR(고반복 수렴) 동일성 + 서브스텝 속도"""
import sys, time; sys.path[:0] = ["legacy", "."]
import numpy as np, core_v17 as c17, discharge_core as c18
P = dict(gap_mm=10.0, voltage_kV=30.0, polarity=+1, nx=161, n_background=0.0, seed_density=0.0,
         t_front_us=0.01, t_tail_us=0.1, snap_ns=1.0, plate_thick_mm=1.0, cath_J0=1e-2,
         n_ion_bg=1e9, q_bg=1e7, config="NEEDLE_PLANE")
def res17(s):
    src = (s.h**2) * (c18.QE*(s.n_ip-s.n_e-s.n_in) + (np.pad(s.sigma[None,:]/s.h, ((s.j_s,s.ny-s.j_s-1),(0,0))) if s.plate_mode=="DIELECTRIC" else 0)) / c18.EPS0
    V=s.V; nbs=(s.cN*np.roll(V,-1,0)+s.cS*np.roll(V,1,0)+s.cE*np.roll(V,-1,1)+s.cW*np.roll(V,1,1))
    m=(s.red|s.black)&(~s.gnd); res=s.csum*V-nbs-src
    return float(np.abs(res[m]).max()/max(np.abs(s.csum*V)[m].max(),1e-300))
for cfg, extra in [("NEEDLE_PLANE", {}), ("NEEDLE_PLANE", {"plate_mode": "DIELECTRIC", "plate_epsr": 3.8}),
                   ("NEEDLE_PLANE", {"plate_mode": "FLOAT_METAL"}), ("NEEDLE_ROGOWSKI", {}),
                   ("SPHERE_SPHERE", {}), ("THREE_NEEDLE", {}), ("NEEDLE_NEEDLE", {}), ("PLANE_PLANE", {})]:
    p = dict(P, config=cfg, **extra)
    t0 = time.perf_counter(); s18 = c18.NeedlePlaneDischarge(**p); ti = time.perf_counter() - t0
    # 같은 상태(전하분포 0에 가까움)에서 v18 직접해 vs v17 SOR 5000회
    s17 = c17.NeedlePlaneDischarge(**p)
    s18.V_plate = s17.V_plate = 1234.0 if extra.get("plate_mode") == "FLOAT_METAL" else 0.0
    s18.Va_now = s17.Va_now = s18.Vpeak      # 피크전압으로 통일
    # 인위적 공간전하: 침 선단 아래 가우시안 양전하 (직접해법의 소스항 검증)
    g = 1e19*np.exp(-((s18.X**2 + (s18.Y - s18.y_tip + 1e-3)**2)/(2*(0.3e-3)**2)))
    for s in (s18, s17):
        s.n_ip[:] = g; s.n_e[:] = 0.0; s.n_in[:] = 0.0; s._zero_in_solids()
        if extra.get("plate_mode") == "DIELECTRIC": s.sigma[:] = 1e-5*np.exp(-(s.x/2e-3)**2)
    s17._apply_bc(); s17.solve_field(iters=6000, w=1.95)
    t0 = time.perf_counter(); s18._apply_bc(); s18.solve_field(); tsolve = time.perf_counter() - t0
    dV = np.abs(s18.V - s17.V).max(); Vs = np.abs(s17.V).max()
    print(f"{cfg:16s} {extra.get('plate_mode','GND'):11s} grid {s18.nx}x{s18.ny} unk={s18.n_unknowns:6d} "
          f"init {ti:5.2f}s factor {s18.t_wall_factor*1e3:6.1f}ms solve {tsolve*1e3:5.2f}ms | "
          f"max|V18-V17(SOR6000)|/max|V| = {dV/Vs:.2e}  residual v18={s18.poisson_residual():.1e} v17={res17(s17):.1e}")
