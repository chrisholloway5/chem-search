@echo off
echo.
echo ============================================================
echo 🧪 CHEMISTRY PLATFORM - STATUS CHECK
echo ============================================================
echo.

REM Check if XAMPP is running
echo [1/4] Checking XAMPP Services...
tasklist /FI "IMAGENAME eq httpd.exe" 2>NUL | find /I /N "httpd.exe" >NUL
if "%ERRORLEVEL%"=="0" (
    echo ✅ Apache is running
) else (
    echo ❌ Apache not running - Run start_xampp.bat
)

tasklist /FI "IMAGENAME eq mysqld.exe" 2>NUL | find /I /N "mysqld.exe" >NUL
if "%ERRORLEVEL%"=="0" (
    echo ✅ MySQL is running
) else (
    echo ❌ MySQL not running - Run start_xampp.bat
)

echo.
echo [2/4] Testing Web Access...
powershell -Command "(Invoke-WebRequest -Uri 'http://localhost/chem-search' -UseBasicParsing).StatusCode" >nul 2>&1
if "%ERRORLEVEL%"=="0" (
    echo ✅ Website accessible at http://localhost/chem-search
) else (
    echo ❌ Website not accessible
)

echo.
echo [3/4] Checking Python Packages...
cd project
.\.venv\Scripts\python.exe -c "import rdkit; print('✅ RDKit:', rdkit.__version__)" 2>nul || echo ❌ RDKit not available
.\.venv\Scripts\python.exe -c "import openai; print('✅ OpenAI:', openai.__version__)" 2>nul || echo ❌ OpenAI not available

echo.
echo [4/4] Running Full Status Check...
.\.venv\Scripts\python.exe platform_status.py

echo.
echo ============================================================
echo 🚀 QUICK LINKS:
echo    Main Site: http://localhost/chem-search
echo    Test Page: http://localhost/chem-search/test_setup.php
echo    Molecule Analyzer: http://localhost/chem-search/tools/molecule_analyzer.php
echo ============================================================
echo.
pause
