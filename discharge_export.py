# -*- coding: utf-8 -*-
"""discharge_export.py - 시뮬레이션 결과 CSV 저장 (GUI 'Save CSV' 와 배치 실행기 공용)
v17 DischargeHMI.export_csv 의 파일 형식을 그대로 유지한다."""
import os
import time
import numpy as np
from discharge_core import QE, RAD_SPECIES


def export_csv_files(s, d, stamp=None):
    """s: NeedlePlaneDischarge, d: 저장 폴더. 반환: 저장한 파일 경로 목록"""
    os.makedirs(d, exist_ok=True)
    stamp = stamp or time.strftime("%Y%m%d_%H%M%S")
    files = []

    def _p(name):
        p = os.path.join(d, f"{name}_{stamp}.csv")
        files.append(p)
        return p

    with open(_p("axial_profiles"), "w", encoding="utf-8-sig") as f:
        f.write("t_ns,y_mm,E_kV_per_cm,ne_m-3,nip_m-3,rho_C_per_m3\n")
        for (tsn, E_ax, ne_ax, nip_ax, rho_ax) in s.profiles:
            for j in range(s.ny):
                e = "" if np.isnan(E_ax[j]) else f"{E_ax[j]/1e5:.6e}"
                f.write(f"{tsn*1e9:.4f},{s.y[j]*1e3:.4f},{e},"
                        f"{ne_ax[j]:.6e},{nip_ax[j]:.6e},{rho_ax[j]:.6e}\n")

    ic = s.nx // 2
    with open(_p("axial_radicals"), "w", encoding="utf-8-sig") as f:
        f.write("y_mm," + ",".join(f"{sp}_m-3" for sp in RAD_SPECIES) + "\n")
        for j in range(s.ny):
            f.write(f"{s.y[j]*1e3:.4f}," +
                    ",".join(f"{s.rad[sp][j, ic]:.6e}" for sp in RAD_SPECIES) + "\n")

    def save2d(name, arr, fmt="%.6e"):
        with open(_p(name), "w", encoding="utf-8-sig") as f:
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

    with open(_p("radical_totals"), "w", encoding="utf-8-sig") as f:
        f.write("t_ns," + ",".join(RAD_SPECIES) + "\n")
        for (tsn, tot) in s.rad_hist:
            f.write(f"{tsn*1e9:.5f}," +
                    ",".join(f"{tot[sp]:.6e}" for sp in RAD_SPECIES) + "\n")
    with open(_p("streak_ne_axis"), "w", encoding="utf-8-sig") as f:
        f.write("y_mm\\t_ns," +
                ",".join(f"{p[0]*1e9:.4f}" for p in s.streak) + "\n")
        for j in range(s.ny):
            f.write(f"{s.y[j]*1e3:.4f}," +
                    ",".join(f"{p[1][j]:.4e}" for p in s.streak) + "\n")
    with open(_p("streak_emission"), "w", encoding="utf-8-sig") as f:
        f.write("y_mm\\t_ns," +
                ",".join(f"{p[0]*1e9:.4f}" for p in s.streak_em) + "\n")
        for j in range(s.ny):
            f.write(f"{s.y[j]*1e3:.4f}," +
                    ",".join(f"{p[1][j]:.4e}" for p in s.streak_em) + "\n")

    with open(_p("initiation_stats"), "w", encoding="utf-8-sig") as f:
        f.write("t_ns,Lambda_per_s,H_cumhazard,F_probability\n")
        for (tsn, lam, hh, ff) in s.init_stats:
            f.write(f"{tsn*1e9:.5f},{lam:.6e},{hh:.6e},{ff:.6e}\n")

    tt, vv = s.waveform_curve(800)
    with open(_p("waveform"), "w", encoding="utf-8-sig") as f:
        f.write("t_s,V_applied_V\n")
        for a_, b_ in zip(tt, vv):
            f.write(f"{a_:.6e},{b_:.6e}\n")

    with open(_p("conditions"), "w", encoding="utf-8-sig") as f:
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
                ("h_grid", s.h*1e3, "mm"),
                ("poisson_solver", s.poisson_used, "-"),
                ("dt_rule", s.dt_rule, "-"),
                ("substeps_total", s.n_substeps, "-"),
                ("wall_time_step_total", s.t_wall_step, "s"),
                ("wall_time_field_total", s.t_wall_field, "s")]
        for sp in RAD_SPECIES:
            rows.append((f"max_{sp}", float(s.rad[sp].max()), "m-3"))
        for k, v, u in rows:
            f.write(f"{k},{v},{u}\n")
    return files
