# -*- coding: utf-8 -*-
"""
==============================================================================
 고전압 방전현상 시뮬레이터 v18 - 통합 실행 진입점 (main.py)
==============================================================================
 배포용 실행 파일(HVDischargeSim.exe)이 호출하는 진입점이며, 파이썬으로 직접
 실행할 때도 동일하게 동작한다.

   python main.py                 -> GUI 실행 (기본)
   python main.py --batch [옵션]  -> 헤드리스 배치 실행 (run_batch.py 의 모든 옵션)
   python main.py --selftest      -> 의존성/코어 자체점검 (실행 파일 검증용)
   python main.py --gui-selftest [초] -> GUI 를 실제로 띄워 N 초 계산 후 8개 탭 그림 저장
   python main.py --version       -> 버전·의존성 정보 출력

 실행 파일 빌드:  python build_exe.py            (build_exe.bat / build_exe.sh)
==============================================================================
"""
import os
import sys

APP_NAME = "HV Discharge Simulator"
APP_VERSION = "18.0"


def is_frozen():
    """PyInstaller 등으로 패키징된 실행 파일로 동작 중인가"""
    return getattr(sys, "frozen", False)


def app_dir():
    """실행 파일(또는 스크립트)이 놓인 폴더 - 결과 저장 기본 위치"""
    if is_frozen():
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


def _ensure_import_path():
    """스크립트 실행 시 같은 폴더의 모듈을 import 할 수 있게 보장"""
    d = os.path.dirname(os.path.abspath(__file__))
    if d not in sys.path:
        sys.path.insert(0, d)


def version_text():
    lines = [f"{APP_NAME} v{APP_VERSION}",
             f"  python      : {sys.version.split()[0]}",
             f"  frozen exe  : {'yes (' + os.path.basename(sys.executable) + ')' if is_frozen() else 'no'}",
             f"  app folder  : {app_dir()}"]
    for mod in ("numpy", "scipy", "matplotlib"):
        try:
            m = __import__(mod)
            lines.append(f"  {mod:<12}: {getattr(m, '__version__', 'ok')}")
        except Exception as e:                      # pragma: no cover
            lines.append(f"  {mod:<12}: MISSING ({e})")
    try:
        import tkinter
        lines.append(f"  {'tkinter':<12}: Tk {tkinter.TkVersion}")
    except Exception as e:                          # pragma: no cover
        lines.append(f"  {'tkinter':<12}: MISSING ({e})")
    return "\n".join(lines)


def selftest():
    """실행 파일이 제대로 묶였는지 확인: 의존성 + 코어 몇 스텝 + 그림 저장"""
    import time
    print(version_text())
    _ensure_import_path()
    import numpy as np
    from discharge_core import NeedlePlaneDischarge, HAVE_SCIPY
    print(f"  poisson     : {'sparse direct (scipy)' if HAVE_SCIPY else 'SOR fallback (scipy missing)'}")
    t0 = time.perf_counter()
    s = NeedlePlaneDischarge(gap_mm=10.0, voltage_kV=30.0, nx=81, n_background=0.0,
                             seed_density=0.0, snap_ns=1.0, plate_thick_mm=1.0)
    t_init = time.perf_counter() - t0
    t0 = time.perf_counter()
    for _ in range(5):
        s.step(n_sub=6)
    t_step = time.perf_counter() - t0
    assert np.isfinite(s.n_e).all() and s.t > 0, "core produced non-finite state"
    print(f"  core init   : {t_init:.2f} s (grid {s.nx}x{s.ny}, solver {s.poisson_used})")
    print(f"  30 substeps : {t_step:.2f} s -> t = {s.t*1e9:.3f} ns, max ne = {s.n_e.max():.2e} m^-3")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = os.path.join(app_dir(), "selftest_field.png")
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.imshow(np.nan_to_num(s.Emag_masked()) / 1e5, origin="lower", cmap="inferno")
    ax.set_title("selftest: |E| [kV/cm]")
    fig.tight_layout()
    fig.savefig(out, dpi=90)
    plt.close(fig)
    print(f"  figure      : {out}")
    print("SELFTEST OK")
    return 0


def gui_selftest(seconds=20.0):
    """실행 파일에 묶인 GUI 를 실제로 띄워 N 초 동안 계산시키고, 8개 탭 그림을 저장한다.
    (화면 없는 환경에서는  xvfb-run -a <exe> --gui-selftest 20  으로 확인)"""
    import time
    _ensure_import_path()
    import tkinter as tk
    import needle_plane_discharge_hmi as H

    print(version_text())
    out = os.path.join(app_dir(), "gui_selftest")
    os.makedirs(out, exist_ok=True)
    root = tk.Tk()
    app = H.DischargeHMI(root)
    root.update()
    s = app.sim
    assert s is not None, "simulation was not created"
    tabs = [("tab1_field", app.figE), ("tab2_ne", app.figN), ("tab3_zoom", app.figZ),
            ("tab4_axial", app.figP), ("tab5_waveform", app.figW), ("tab6_radicals", app.figR),
            ("tab7_streak", app.figS), ("tab8_initiation", app.figT)]
    print(f"  window      : {root.winfo_width()}x{root.winfo_height()}, tabs {len(tabs)}, "
          f"grid {s.nx}x{s.ny}, solver {s.poisson_used}")

    t0 = time.perf_counter()
    app.toggle_run()
    root.after(int(float(seconds) * 1000), root.quit)
    root.mainloop()
    if app.running:
        app.toggle_run()
    root.update()
    W = time.perf_counter() - t0
    print(f"  run {W:.1f} s : t = {s.t*1e9:.3f} ns, {s.n_substeps} substeps, "
          f"{app.frame_count} frames, redraw {app._t_draw:.2f} s, bridged = {s.bridged}")
    for name, fig in tabs:
        fig.savefig(os.path.join(out, name + ".png"), dpi=100)
    print(f"  figures     : {out} ({len(tabs)} PNG)")
    assert s.n_substeps > 0, "no substep was computed"
    root.destroy()
    print("GUI SELFTEST OK")
    return 0


def run_gui():
    _ensure_import_path()
    import tkinter as tk
    from tkinter import ttk, messagebox
    import needle_plane_discharge_hmi as H

    root = tk.Tk()
    try:
        st = ttk.Style()
        st.theme_use("clam")
        for w in ("TLabel", "TButton", "TRadiobutton", "TEntry",
                  "TCombobox", "TLabelframe.Label", "TLabelframe"):
            st.configure(w, font=H.FONT_UI)
        root.option_add("*Font", "{Times New Roman} 12")
    except Exception:
        pass
    if not H.HAVE_SCIPY:
        messagebox.showwarning(
            "scipy not found",
            "scipy is not installed: falling back to the slow v17 SOR Poisson solver.\n"
            "Install it for the fast direct solver:  pip install scipy")
    H.DischargeHMI(root)
    root.mainloop()
    return 0


def run_batch_cli(argv):
    _ensure_import_path()
    import run_batch
    # 실행 파일에서 켰을 때 결과가 exe 옆 results/ 로 가도록 작업 폴더를 고정
    if is_frozen():
        os.chdir(app_dir())
    run_batch.main(argv)
    return 0


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in ("--version", "-V"):
        print(version_text())
        return 0
    if argv and argv[0] == "--selftest":
        return selftest()
    if argv and argv[0] == "--gui-selftest":
        return gui_selftest(argv[1] if len(argv) > 1 else 20.0)
    if argv and argv[0] == "--batch":
        return run_batch_cli(argv[1:])
    if argv and argv[0] in ("--help", "-h"):
        print(__doc__)
        return 0
    if argv:
        print(f"unknown option: {argv[0]}\n")
        print(__doc__)
        return 2
    return run_gui()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception:                               # pragma: no cover
        import traceback
        msg = traceback.format_exc()
        sys.stderr.write(msg)
        # 콘솔 없는 GUI 실행 파일에서도 원인을 볼 수 있도록 로그 + 대화상자
        try:
            with open(os.path.join(app_dir(), "error_log.txt"), "w", encoding="utf-8") as f:
                f.write(msg)
        except Exception:
            pass
        try:
            import tkinter as tk
            from tkinter import messagebox
            r = tk.Tk()
            r.withdraw()
            messagebox.showerror(f"{APP_NAME} - error",
                                 msg[-2000:] + "\n\n(error_log.txt saved next to the program)")
            r.destroy()
        except Exception:
            pass
        sys.exit(1)
