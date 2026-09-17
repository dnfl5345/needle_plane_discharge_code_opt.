# -*- coding: utf-8 -*-
"""실행 진입점(main.py) 회귀 테스트:  python tests/test_entrypoint.py
 - --version / --selftest / --batch 가 하위 프로세스로 정상 동작하는지 확인
 - 빌드된 실행 파일도 같은 방법으로 점검할 수 있다:
       python tests/test_entrypoint.py dist/HVDischargeSim/HVDischargeSim
   (GUI 포함 점검은  xvfb-run -a <exe> --gui-selftest 15)
"""
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def build_cmd(target, args):
    """target 이 실행 파일이면 그대로, 아니면 'python main.py' 로 호출"""
    if target:
        return [target] + args
    return [sys.executable, os.path.join(ROOT, "main.py")] + args


def run(target, args, timeout=600):
    p = subprocess.run(build_cmd(target, args), cwd=ROOT, capture_output=True,
                       text=True, timeout=timeout)
    out = p.stdout + p.stderr
    assert p.returncode == 0, f"exit {p.returncode} for {args}\n{out[-3000:]}"
    return out


def main():
    target = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else None
    what = os.path.basename(target) if target else "main.py"

    out = run(target, ["--version"])
    for key in ("HV Discharge Simulator", "numpy", "scipy", "matplotlib", "tkinter"):
        assert key in out, f"--version missing {key}\n{out}"
    assert "MISSING" not in out, f"dependency missing:\n{out}"
    print(f"[ok] {what} --version")

    out = run(target, ["--selftest"])
    assert "SELFTEST OK" in out, out[-2000:]
    assert "sparse direct" in out, f"scipy solver not active:\n{out}"
    print(f"[ok] {what} --selftest")

    with tempfile.TemporaryDirectory() as d:
        out = run(target, ["--batch", "--t_end_ns", "2", "--nx", "81",
                           "--out", d, "--no_csv", "--log_s", "60"])
        assert "SUMMARY" in out, out[-2000:]
        figs = [f for f in os.listdir(d) if f.endswith(".png")]
        assert len(figs) >= 6, f"only {figs} written"
        assert "timing.json" in os.listdir(d)
    print(f"[ok] {what} --batch (figures: {len(figs)})")

    out = run(target, ["--help"])
    assert "--gui-selftest" in out
    print(f"[ok] {what} --help")
    print("ALL ENTRYPOINT TESTS PASSED")


if __name__ == "__main__":
    main()
