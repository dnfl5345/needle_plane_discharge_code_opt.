#!/usr/bin/env bash
# ============================================================================
#  Linux / macOS 실행 파일 빌드
#    ./build_exe.sh              폴더형 (기본)
#    ./build_exe.sh --onefile    단일 실행 파일
#  * tkinter 가 포함된 파이썬이 필요합니다 (Ubuntu: sudo apt install python3-tk)
# ============================================================================
set -e
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
"$PY" -c "import tkinter" 2>/dev/null || {
  echo "tkinter 가 없는 파이썬입니다: $PY"
  echo "  Ubuntu/Debian: sudo apt install python3-tk"
  echo "  다른 파이썬 사용: PYTHON=/usr/bin/python3.12 ./build_exe.sh"
  exit 1
}
echo "[1/2] 의존성 확인/설치 ..."
"$PY" -m pip install -r requirements.txt pyinstaller
echo "[2/2] 빌드 ..."
"$PY" build_exe.py "$@"
