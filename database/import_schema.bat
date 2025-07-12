@echo off
echo Importing Chem Search database schema...
cd /d "C:\xampp\mysql\bin"
mysql.exe -u root chem_search < "C:\xampp\htdocs\Subchem-search\database\schema.sql"
if %errorlevel% equ 0 (
    echo ✓ Database schema imported successfully!
    mysql.exe -u root -e "USE chem_search; SHOW TABLES;"
) else (
    echo ❌ Error importing schema
)
pause
