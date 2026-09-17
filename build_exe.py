# -*- coding: utf-8 -*-
"""
==============================================================================
 실행 파일(.exe) 빌드 스크립트 - build_exe.py
==============================================================================
 GUI 가 포함된 단독 실행 파일을 만든다. 파이썬·numpy·scipy·matplotlib·tkinter 가
 모두 안에 들어가므로 배포 대상 PC 에 파이썬을 설치하지 않아도 실행된다.

   pip install -r requirements.txt pyinstaller
   python build_exe.py              # 폴더형 (기본, 권장: 시작 2~3 초)
   python build_exe.py --onefile    # 단일 파일 (배포 간편, 시작 15~30 초)
   python build_exe.py --console    # 콘솔 창 유지 (오류 확인용)
   python build_exe.py --clean      # build/ dist/ 정리 후 빌드

 결과물
   폴더형 : dist/HVDischargeSim/HVDischargeSim(.exe)   <- 폴더 전체를 배포
   단일형 : dist/HVDischargeSim(.exe)                  <- 파일 하나만 배포

 * 실행 파일은 만들어지는 OS 용으로만 동작한다(Windows 용 .exe 는 Windows 에서 빌드).
 * 빌드 후 자체점검:  dist/.../HVDischargeSim --selftest
==============================================================================
"""
import os
import shutil
import sys
import time

NAME = "HVDischargeSim"
ENTRY = "main.py"
HERE = os.path.dirname(os.path.abspath(__file__))

# 실행 파일에 반드시 포함할 프로젝트 모듈 (main.py 가 지연 import 하므로 명시)
PROJECT_MODULES = ["discharge_core", "needle_plane_discharge_hmi", "run_batch", "discharge_export"]

# PyInstaller 가 자동으로 찾지 못할 수 있는 동적 import
HIDDEN = ["scipy.sparse.linalg", "scipy.sparse.csgraph", "scipy._lib.messagestream",
          "matplotlib.backends.backend_tkagg", "matplotlib.backends.backend_agg",
          "tkinter", "tkinter.filedialog", "tkinter.messagebox", "tkinter.ttk"]

# 쓰지 않는 무거운 패키지 제외 (용량·빌드시간 절약)
EXCLUDE = ["PyQt5", "PyQt6", "PySide2", "PySide6", "wx", "IPython", "jupyter", "notebook",
           "pytest", "pandas", "sphinx", "docutils", "PIL.ImageQt", "matplotlib.backends.backend_qt5agg",
           "matplotlib.backends.backend_webagg", "scipy.spatial.cKDTree"]


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    onefile = "--onefile" in argv
    console = "--console" in argv
    clean = "--clean" in argv
    for a in argv:
        if a not in ("--onefile", "--console", "--clean"):
            print(f"unknown option: {a}\n{__doc__}")
            return 2

    try:
        import PyInstaller.__main__ as pyi
        import PyInstaller
    except ImportError:
        print("PyInstaller 가 없습니다.  pip install pyinstaller  후 다시 실행하세요.")
        return 1
    for mod in ("numpy", "scipy", "matplotlib"):
        try:
            __import__(mod)
        except ImportError:
            print(f"{mod} 가 없습니다.  pip install -r requirements.txt  후 다시 실행하세요.")
            return 1
    try:
        import tkinter                                   # noqa: F401
    except ImportError:
        print("tkinter 가 없습니다. Windows 는 파이썬 설치 시 'tcl/tk' 옵션을, "
              "Linux 는 python3-tk 패키지를 설치하세요.")
        return 1

    if clean:
        for d in ("build", "dist"):
            p = os.path.join(HERE, d)
            if os.path.isdir(p):
                print(f"removing {p}")
                shutil.rmtree(p, ignore_errors=True)

    args = [os.path.join(HERE, ENTRY),
            "--name", NAME,
            "--noconfirm",
            "--distpath", os.path.join(HERE, "dist"),
            "--workpath", os.path.join(HERE, "build"),
            "--specpath", HERE,
            "--paths", HERE]
    args += ["--onefile"] if onefile else ["--onedir"]
    # GUI 프로그램: 콘솔 창 없이 실행 (--console 로 오류 메시지 확인 가능)
    args += ["--console"] if console else ["--windowed"]
    for m in PROJECT_MODULES + HIDDEN:
        args += ["--hidden-import", m]
    for m in EXCLUDE:
        args += ["--exclude-module", m]

    print(f"PyInstaller {PyInstaller.__version__} | mode: "
          f"{'onefile' if onefile else 'onedir'}, {'console' if console else 'windowed'}")
    print(" ".join(args))
    t0 = time.perf_counter()
    pyi.run(args)
    dt = time.perf_counter() - t0

    exe = os.path.join(HERE, "dist", NAME if onefile else os.path.join(NAME, NAME))
    if os.name == "nt":
        exe += ".exe"
    ok = os.path.exists(exe)
    if ok:
        if onefile:
            size = os.path.getsize(exe)
        else:
            size = sum(os.path.getsize(os.path.join(r, f))
                       for r, _, fs in os.walk(os.path.join(HERE, "dist", NAME)) for f in fs)
        print(f"\nBUILD OK ({dt:.0f} s)\n  {exe}\n  size: {size/1e6:.0f} MB"
              f"{' (folder total)' if not onefile else ''}")
        print(f"\n자체점검:  \"{exe}\" --selftest")
        print(f"GUI 실행 :  \"{exe}\"")
        print(f"배치 실행:  \"{exe}\" --batch --afterglow_us 10")
    else:
        print(f"\nBUILD FAILED: {exe} 가 생성되지 않았습니다 ({dt:.0f} s)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
