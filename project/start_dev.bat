@echo off
cd /d "C:\xampp\htdocs\chemcrow\project"
echo.
echo ===================================
echo   Chem Development Environment
echo ===================================
echo.
echo Current Directory: %CD%
echo Python Environment: C:\xampp\htdocs\chemcrow\project\.venv\Scripts\python.exe
echo Website URL: http://localhost/chemcrow/
echo.
echo Available Commands:
echo   python --version    - Check Python version
echo   pip list           - List installed packages
echo   dir                - List directory contents
echo.
cmd /k
