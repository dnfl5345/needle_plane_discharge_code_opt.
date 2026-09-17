# -*- coding: utf-8 -*-
"""헤드리스 회귀 테스트:  python tests/test_core.py
 1) 직접해법 == v17 SOR 수렴해 (침-평판 3종 평판 모드)
 2) LUT 반응률 == v17 log-log 보간 (상대오차 < 0.5%)
 3) v18(dt 규칙 v17) 1 서브스텝 == v17 원본 코어 1 서브스텝 (dt 상대차 < 1e-3, 밀도 상대차 < 2e-3; LUT 보간오차 수준)
 4) 브리징/정체 가드 및 회로(R>0) 동작, 반암시적 스텝 동작
"""
import os, sys, time
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "legacy"), ROOT]
import numpy as np
import core_v17 as c17
import discharge_core as c18

P = dict(gap_mm=10.0, voltage_kV=30.0, polarity=+1, nx=121, n_background=0.0, seed_density=0.0,
         t_front_us=0.01, t_tail_us=0.1, snap_ns=1.0, plate_thick_mm=1.0, cath_J0=1e-2,
         n_ion_bg=1e9, q_bg=1e7, config="NEEDLE_PLANE")

def test_poisson_direct_equals_converged_sor():
    for extra in ({}, {"plate_mode": "DIELECTRIC", "plate_epsr": 3.8}, {"plate_mode": "FLOAT_METAL"}):
        p = dict(P, **extra)
        s18 = c18.NeedlePlaneDischarge(**p); s17 = c17.NeedlePlaneDischarge(**p)
        g = 1e19 * np.exp(-((s18.X**2 + (s18.Y - s18.y_tip + 1e-3)**2) / (2 * (0.3e-3)**2)))
        for s in (s18, s17):
            s.Va_now = s.Vpeak; s.V_plate = 500.0 if extra.get("plate_mode") == "FLOAT_METAL" else 0.0
            s.n_ip[:] = g; s.n_e[:] = 0; s.n_in[:] = 0; s._zero_in_solids()
            if extra.get("plate_mode") == "DIELECTRIC":
                s.sigma[:] = 1e-5 * np.exp(-(s.x / 2e-3)**2)
        s17._apply_bc(); s17.solve_field(iters=8000, w=1.95)
        s18._apply_bc(); s18.solve_field()
        err = np.abs(s18.V - s17.V).max() / np.abs(s17.V).max()
        assert err < 1e-8, (extra, err)
        assert s18.poisson_residual() < 1e-12
    print("[ok] direct Poisson == converged SOR")

def test_lut_matches_v17_rates():
    s17 = c17.NeedlePlaneDischarge(**P); s18 = c18.NeedlePlaneDischarge(**P)
    EN = np.logspace(-1, 3.5, 5000)
    for name in ["k_ion_N2", "k_ion_O2", "k_att_O2", "k_det_O2m", "k_exc_N2", "muN", "k_diss_N2", "k_diss_O2_8eV"]:
        err = np.abs(s18._k(name, EN) / s17._k(name, EN) - 1).max()
        assert err < 5e-3, (name, err)
    print("[ok] LUT rates within 0.5% of v17 interpolation")

def test_single_substep_matches_v17():
    s17 = c17.NeedlePlaneDischarge(**P); s18 = c18.NeedlePlaneDischarge(**P, dt_rule="v17", semi_implicit=False)
    s17.solve_field(iters=6000, w=1.95)                         # v17 초기장도 수렴시켜 동일 출발
    for a in ("V", "Ex", "Ey", "Emag"):
        getattr(s18, a)[:] = getattr(s17, a)
    for k in range(6):
        for a in ("V", "n_e", "n_ip", "n_in", "Ex", "Ey", "Emag"):
            getattr(s18, a)[:] = getattr(s17, a)
        for sp in c18.RAD_SPECIES:
            s18.rad[sp][:] = s17.rad[sp]
        s18.t = s17.t; s18.Va_now = s17.Va_now
        s17.step(n_sub=1, field_iters=4000); s18.step(n_sub=1)
        assert abs(s17.dt_last - s18.dt_last) / s17.dt_last < 1e-3, (s17.dt_last, s18.dt_last)
        rel = np.abs(s17.n_e - s18.n_e).max() / max(s17.n_e.max(), 1e-300)
        assert rel < 2e-3, (k, rel)
        assert np.abs(s17.V - s18.V).max() < 1.0, np.abs(s17.V - s18.V).max()   # SOR 잔차 수준(V)
    print("[ok] single-substep evolution matches v17 (converged field)")

def test_guards_circuit_semi():
    # 짧은 갭/고전압으로 빠르게 브리징: 가드·회로·반암시적 경로 동작 확인
    p = dict(P, gap_mm=3.0, voltage_kV=25.0, nx=81)
    s = c18.NeedlePlaneDischarge(**p, R_series_ohm=2000.0)
    st = s.run_until(40e-9, n_sub=6, stop_on_bridge=True, post_bridge=1e-9)
    assert s.bridged and st["stop_reason"] == "bridged", st
    assert s.n_semi_steps > 0 and s.n_semi_factor > 0
    assert 0.0 <= s.Va_now <= s.V_src_now and s.I_dis > 0
    assert s.n_e.max() <= s.n_cap * (1 + 1e-12)
    assert np.isfinite(s.V).all() and np.nanmax(s.Emag_masked()) < 1e9   # 전계 발산 없음
    assert not s.stalled
    print(f"[ok] guards/circuit/semi: bridged at {s.t_bridge*1e9:.2f} ns, semi steps {s.n_semi_steps}, "
          f"I={s.I_dis:.2f} A, Va={s.Va_now/1e3:.2f} kV")

if __name__ == "__main__":
    t0 = time.perf_counter()
    test_poisson_direct_equals_converged_sor()
    test_lut_matches_v17_rates()
    test_single_substep_matches_v17()
    test_guards_circuit_semi()
    print(f"ALL TESTS PASSED ({time.perf_counter()-t0:.0f} s)")
