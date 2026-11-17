@echo off
setlocal
chcp 65001 >nul
cd /d "%~dp0"
title XRD Analyzer - naprawa srodowiska

if exist ".venv" (
  echo Usuwanie starego srodowiska .venv...
  rmdir /s /q ".venv"
)

echo Srodowisko usuniete. Launcher utworzy je od nowa.
call "%~dp0START_XRD_ANALYZER.bat"
