@echo off
title Protein Atlas
cd /d "%~dp0"

echo.
echo  ===========================================
echo    Protein Atlas - launching
echo  ===========================================
echo.

if not exist "api.py" (
  echo ERROR: api.py is missing from this folder. Extract the WHOLE zip first ^(right-click, Extract All^).
  pause
  exit /b 1
)
if not exist "frontend\index.html" (
  echo ERROR: the frontend folder is missing. Extract the WHOLE zip first ^(right-click, Extract All^).
  pause
  exit /b 1
)

set "PY="
where py >nul 2>nul && set "PY=py -3"
if not defined PY ( where python >nul 2>nul && set "PY=python" )
if not defined PY (
  echo ERROR: Python was not found. Install it from https://www.python.org/downloads/
  echo Tick "Add python.exe to PATH" during install, then run this file again.
  pause
  exit /b 1
)

if not exist ".venv\Scripts\python.exe" (
  echo Setting up - first run only, takes a minute or two...
  %PY% -m venv .venv
  if errorlevel 1 ( echo Could not create the virtual environment. & pause & exit /b 1 )
)

echo Checking packages...
".venv\Scripts\python.exe" -m pip install --quiet --disable-pip-version-check -r requirements.txt
if errorlevel 1 (
  echo ERROR: Package install failed. Check your internet connection and run again.
  pause
  exit /b 1
)

echo.
echo  Starting server. Chrome will open in a few seconds.
echo  First run only: the PROSITE database and the dataset statistics download in the
echo  background (a few minutes) - the site is usable straight away.
echo  Keep this window open while you use the site. Close it to stop.
echo.

start "" cmd /c "timeout /t 5 /nobreak >nul & start chrome http://127.0.0.1:8000/"

".venv\Scripts\python.exe" -m uvicorn api:app --host 127.0.0.1 --port 8000
pause
