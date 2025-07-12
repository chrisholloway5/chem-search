@echo off
echo Starting XAMPP services for Chem Search...
echo.

REM Start XAMPP Control Panel
echo Opening XAMPP Control Panel...
start "" "C:\xampp\xampp-control.exe"

echo.
echo XAMPP Control Panel opened.
echo.
echo Manual steps to complete:
echo 1. In XAMPP Control Panel, click "Start" for Apache
echo 2. If MySQL is not running, click "Start" for MySQL
echo 3. Wait for both services to show "Running" status
echo 4. Open browser and go to: http://localhost/chem-search
echo.
echo Junction created: C:\xampp\htdocs\chem-search
echo Points to: %~dp0
echo.
pause
