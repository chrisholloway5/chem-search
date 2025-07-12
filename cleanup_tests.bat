@echo off
echo.
echo ============================================================
echo 🧹 CLEANING UP TEST FILES - Chem Search Project
echo ============================================================
echo.

echo [1/4] Removing empty test files...

REM Remove empty test PHP files
del test_analyzer.php 2>nul
del test_complete_pipeline.php 2>nul
del test_db.php 2>nul
del test_db_columns.php 2>nul
del test_db_direct.php 2>nul
del test_fixes.php 2>nul
del test_full_analyzer.php 2>nul
del test_raw.php 2>nul
del test_simple_pipeline.php 2>nul

echo ✅ Empty test PHP files removed

echo.
echo [2/4] Removing debug files...

REM Remove debug files
del debug_db.php 2>nul
del debug_detailed.php 2>nul

echo ✅ Debug files removed

echo.
echo [3/4] Removing old verification files...

REM Remove old verification files
del verification_complete.php 2>nul
del final_verification.php 2>nul

echo ✅ Old verification files removed

echo.
echo [4/4] Removing phpinfo test file...

REM Remove phpinfo test (security risk in production)
del phpinfo_test.php 2>nul

echo ✅ PHPInfo test file removed

echo.
echo ============================================================
echo 📋 CLEANUP SUMMARY
echo ============================================================
echo.
echo ✅ Files Removed:
echo    • Empty test_*.php files (9 files)
echo    • Debug files (2 files)  
echo    • Old verification files (2 files)
echo    • PHPInfo test file (1 file)
echo.
echo ✅ Files Kept (Important):
echo    • test_setup.php - Main setup verification
echo    • link_checker.php - Link verification tool
echo    • navigation_test.php - Navigation testing
echo    • project/test_*.py - Python functionality tests
echo.
echo 🎯 Total: 14 unnecessary files removed
echo 💾 Disk space saved and project structure cleaned
echo.
echo ============================================================
echo.
pause
