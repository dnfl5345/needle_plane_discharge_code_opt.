# -*- coding: utf-8 -*-
"""GUI 스모크 테스트 (Xvfb 가상 디스플레이에서 실행):
   xvfb-run -a python3 tests/test_gui_smoke.py [run_seconds]
 - HMI 생성/Apply, 1프레임, Run(벽시계 N초)/Pause, 구성 변경, 잔광 토글, CSV 저장, 8개 탭 figure 확인"""
import os, sys, time, tempfile
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import tkinter as tk
import needle_plane_discharge_hmi as H

run_s = float(sys.argv[1]) if len(sys.argv) > 1 else 60.0
root = tk.Tk()
app = H.DischargeHMI(root)
root.update()
s = app.sim
assert s is not None and s.poisson_used == "direct", s.poisson_used
print(f"[gui] created; grid {s.nx}x{s.ny}, poisson={s.poisson_used}")
figs = [app.figE, app.figN, app.figZ, app.figP, app.figR, app.figW, app.figS, app.figT]
assert all(len(f.axes) > 0 for f in figs), "tabs without axes"
print(f"[gui] 8 tabs, axes per fig: {[len(f.axes) for f in figs]}")

# 1 frame
app.step_once(); root.update()
assert s.t > 0 and app.frame_count == 1
print(f"[gui] step_once OK: t={s.t*1e9:.3f} ns")

# Run for run_s seconds of wall time with redraw every 1 s (실제 Tk 이벤트루프 사용)
app.vars["redraw"].set("1.0"); app.vars["tick"].set("0.25")
t0 = time.perf_counter(); app.toggle_run(); assert app.running
root.after(int(run_s * 1000), root.quit)
root.mainloop()
if app.running:
    app.toggle_run()
root.update()
W = time.perf_counter() - t0
print(f"[gui] run {W:.1f} s: t={s.t*1e9:.3f} ns, substeps={s.n_substeps}, frames={app.frame_count}, "
      f"speed={s.t*1e9/W:.4f} ns/s, last redraw {app._t_draw:.2f} s, bridged={s.bridged} "
      f"(t_bridge={None if s.t_bridge is None else round(s.t_bridge*1e9,3)}), afterglow={s.afterglow}, msg='{app._msg}'")
assert s.n_substeps > 100

# 잔광 토글 + 1프레임
if not s.afterglow:
    app.toggle_afterglow()
assert s.afterglow
app.step_once(); root.update()
print(f"[gui] afterglow step OK: t={s.t*1e9:.3f} ns dt={s.dt_last*1e9:.3f} ns")
app.toggle_afterglow(); assert not s.afterglow

# CSV 저장 (다이얼로그 대체)
d = tempfile.mkdtemp()
H.filedialog.askdirectory = lambda **k: d
H.messagebox.showinfo = lambda *a, **k: None
app.export_csv()
n_csv = len([f for f in os.listdir(d) if f.endswith(".csv")])
assert n_csv >= 18, n_csv
print(f"[gui] export_csv OK: {n_csv} files in {d}")

# 구성 변경 + R>0 + Apply + 몇 프레임
app.cfg.set("Needle - Rogowski plane"); app.vars["rser"].set("500"); app.vars["nx"].set("121")
app.apply_params(); root.update()
s = app.sim
assert s.config == "NEEDLE_ROGOWSKI" and s.R_series == 500 and s.nx == 121
for _ in range(3):
    app.step_once()
root.update()
print(f"[gui] config switch OK: {s.config}, R={s.R_series}, grid {s.nx}x{s.ny}, t={s.t*1e9:.3f} ns, I={s.I_dis:.3e} A")
app.cfg.set("Sphere - Sphere"); app.apply_params(); app.step_once(); root.update()
print(f"[gui] sphere-sphere OK: t={app.sim.t*1e9:.3f} ns")
root.destroy()
print("[gui] ALL OK")
