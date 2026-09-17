@echo off
REM ===========================================================================
REM  Windows 실행 파일(.exe) 빌드 - 더블클릭하거나 명령창에서 실행
REM    build_exe.bat              폴더형 (기본, 권장)
REM    build_exe.bat --onefile    단일 exe
REM ===========================================================================
setlocal
cd /d "%~dp0"
echo [1/2] 의존성 확인/설치 ...
python -m pip install -r requirements.txt pyinstaller || goto :err
echo [2/2] 빌드 ...
python build_exe.py %* || goto :err
echo.
echo 완료. dist 폴더를 확인하세요.
pause
exit /b 0
:err
echo.
echo 빌드에 실패했습니다. 위 메시지를 확인하세요.
pause
exit /b 1
