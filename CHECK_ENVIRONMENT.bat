@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
set "VPY=%~dp0.venv\Scripts\python.exe"
if not exist "%VPY%" (
  echo Brak .venv. Najpierw uruchom START_XRD_ANALYZER.bat
  pause
  exit /b 1
)
"%VPY%" -c "import sys,PyQt6,numpy,scipy,matplotlib,openpyxl,gemmi; print('Python:',sys.executable); print('Version:',sys.version); print('gemmi:',gemmi.__version__); print('numpy:',numpy.__version__); print('scipy:',scipy.__version__); print('matplotlib:',matplotlib.__version__)"
pause
