@echo off
setlocal EnableExtensions EnableDelayedExpansion
chcp 65001 >nul
cd /d "%~dp0"

title XRD Analyzer - uruchamianie programu
set "VENV=%~dp0.venv"
set "VPY=%VENV%\Scripts\python.exe"

echo ============================================================
echo XRD Analyzer - uruchamianie w dedykowanym srodowisku Python
echo ============================================================
echo.

if exist "%VPY%" goto CHECK_ENV

echo [1/4] Tworzenie lokalnego srodowiska .venv...
set "BASEPY="

rem Prefer Python 3.13 because scientific packages are broadly supported.
py -3.13 -c "import sys; print(sys.executable)" >nul 2>&1
if not errorlevel 1 set "BASEPY=py -3.13"
if defined BASEPY goto CREATE_VENV

py -3.14 -c "import sys; print(sys.executable)" >nul 2>&1
if not errorlevel 1 set "BASEPY=py -3.14"
if defined BASEPY goto CREATE_VENV

py -3.12 -c "import sys; print(sys.executable)" >nul 2>&1
if not errorlevel 1 set "BASEPY=py -3.12"
if defined BASEPY goto CREATE_VENV

py -3.11 -c "import sys; print(sys.executable)" >nul 2>&1
if not errorlevel 1 set "BASEPY=py -3.11"
if defined BASEPY goto CREATE_VENV

python -c "import sys; print(sys.executable)" >nul 2>&1
if not errorlevel 1 set "BASEPY=python"

if not defined BASEPY (
    echo BLAD: Nie znaleziono Pythona.
    echo Zainstaluj Python 3.13 z python.org i zaznacz "Add Python to PATH".
    pause
    exit /b 1
)

:CREATE_VENV
echo Uzywam: %BASEPY%
%BASEPY% -m venv "%VENV%"
if errorlevel 1 (
    echo BLAD: Nie udalo sie utworzyc srodowiska .venv.
    pause
    exit /b 1
)

:CHECK_ENV
if not exist "%VPY%" (
    echo BLAD: Brak interpretera %VPY%
    pause
    exit /b 1
)

echo [2/4] Sprawdzanie bibliotek...
"%VPY%" -c "import PyQt6,numpy,scipy,matplotlib,openpyxl,gemmi" >nul 2>&1
if not errorlevel 1 goto RUN_APP

echo [3/4] Instalowanie / aktualizowanie wymaganych bibliotek...
"%VPY%" -m pip install --upgrade pip
if errorlevel 1 goto INSTALL_ERROR
"%VPY%" -m pip install --upgrade -r "%~dp0requirements.txt"
if errorlevel 1 goto INSTALL_ERROR

"%VPY%" -c "import PyQt6,numpy,scipy,matplotlib,openpyxl,gemmi; print('OK - gemmi', gemmi.__version__)"
if errorlevel 1 goto INSTALL_ERROR

:RUN_APP
echo [4/4] Uruchamianie XRD Analyzer...
"%VPY%" "%~dp0XRD_Analyzer_GUI.py"
set "RC=%ERRORLEVEL%"
if not "%RC%"=="0" (
    echo.
    echo Program zakonczyl sie kodem bledu %RC%.
    pause
)
exit /b %RC%

:INSTALL_ERROR
echo.
echo ============================================================
echo BLAD INSTALACJI BIBLIOTEK

echo Interpreter programu:
echo "%VPY%"
echo.
echo Sprobuj uruchomic REPAIR_XRD_ENVIRONMENT.bat.
echo ============================================================
pause
exit /b 1
