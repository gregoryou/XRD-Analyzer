@echo off
setlocal EnableExtensions
cd /d "%~dp0"
set "PY=.venv\Scripts\python.exe"
if not exist "%PY%" (
  echo Brak .venv. Uruchom najpierw START_XRD_ANALYZER.bat, aby utworzyc srodowisko.
  pause
  exit /b 2
)
"%PY%" -m pytest --version >nul 2>&1
if errorlevel 1 (
  echo Instalowanie pytest w wirtualnym srodowisku XRD Analyzer...
  "%PY%" -m pip install pytest
  if errorlevel 1 goto ERROR
)
set QT_QPA_PLATFORM=offscreen
"%PY%" -m pytest -v tests\test_cod_offline_2035.py
if errorlevel 1 goto ERROR
echo.
echo Testy COD Offline 20.3.5: PASS
pause
exit /b 0
:ERROR
echo Testy COD Offline nie powiodly sie.
pause
exit /b 1
