# Database Setup Guide for Subchems Search

## Prerequisites

1. **XAMPP with MySQL/MariaDB** - Make sure MySQL is running
2. **Python with mysql-connector-python** (optional, for automated setup)

## Quick Setup

### Option 1: Automated Setup (Recommended)

1. **Install mysql-connector-python:**
   ```bash
   pip install mysql-connector-python
   ```

2. **Run the setup script:**
   ```bash
   cd C:\xampp\htdocs\chemcrow\database
   python setup.py
   ```

### Option 2: Manual Setup

1. **Open phpMyAdmin** (http://localhost/phpmyadmin)

2. **Create the database:**
   ```sql
   CREATE DATABASE subchems_search 
   CHARACTER SET utf8mb4 
   COLLATE utf8mb4_unicode_ci;
   ```

3. **Import the schema:**
   - Select the `subchems_search` database
   - Go to Import tab
   - Choose file: `C:\xampp\htdocs\chemcrow\database\schema.sql`
   - Click Go

### Option 3: Command Line

1. **Open Command Prompt as Administrator**

2. **Navigate to MySQL bin directory:**
   ```cmd
   cd C:\xampp\mysql\bin
   ```

3. **Create database and import schema:**
   ```cmd
   mysql -u root -p -e "CREATE DATABASE subchems_search CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;"
   mysql -u root -p subchems_search < "C:\xampp\htdocs\chemcrow\database\schema.sql"
   ```

## Configuration

1. **Copy environment file:**
   ```bash
   copy .env.example .env
   ```

2. **Edit .env file** with your database credentials:
   ```
   DB_HOST=localhost
   DB_NAME=subchems_search
   DB_USER=root
   DB_PASS=your_mysql_password
   USE_DATABASE=true
   ```

## Verify Installation

1. **Test database connection:**
   - Visit: http://localhost/chemcrow/tools/database_manager.php
   - You should see the database dashboard

2. **Test molecule analysis with storage:**
   - Visit: http://localhost/chemcrow/tools/molecule_analyzer.php
   - Analyze a molecule (e.g., caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`)
   - Check database manager to see stored data

## Database Structure

The database contains 8 main tables:

### Core Tables
- **`chemicals`** - Basic chemical identifiers (SMILES, InChI, names)
- **`molecular_properties`** - Computed molecular properties
- **`drug_properties`** - Drug-likeness and ADMET properties
- **`physical_properties`** - Physical and thermodynamic properties
- **`safety_properties`** - Safety and toxicology data

### Supporting Tables
- **`chemical_synonyms`** - Alternative names and synonyms
- **`analysis_cache`** - Cached analysis results
- **`search_history`** - User search and activity logs

## Features

### ✅ **Comprehensive Chemical Profiles**
- IUPAC names, InChI, InChI Keys, SMILES notations
- CAS numbers and external database IDs
- Complete molecular property calculations

### ✅ **Advanced Property Calculations**
- Molecular Weight, Exact Mass, Monoisotopic Mass
- XLogP3, CLogP (lipophilicity)
- Hydrogen Bond Donors/Acceptors
- Rotatable Bond Count, TPSA
- Heavy Atom Count, Formal Charge
- Complexity scores and stereochemistry data

### ✅ **Drug-Like Property Assessment**
- Lipinski's Rule of Five compliance
- Veber's Rule violations
- QED (Quantitative Estimate of Drug-likeness)
- ADMET predictions

### ✅ **Caching & Performance**
- Analysis result caching (24-hour default)
- Fast lookups by SMILES, InChI Key, CAS number
- Search history and popular searches tracking

### ✅ **Data Management**
- Export capabilities (JSON, CSV)
- Database statistics and analytics
- Recent analyses tracking
- Search pattern analysis

## Usage Examples

### PHP API Usage
```php
require_once 'includes/chemical_model.php';
$chemical = new Chemical();

// Find or create chemical
$result = $chemical->findOrCreateBySmiles('CCO'); // Ethanol

// Store comprehensive properties
$chemical->storeMolecularProperties($chemicalId, $properties);
$chemical->storeDrugProperties($chemicalId, $drugProps);

// Get full profile
$profile = $chemical->getFullProfile($chemicalId);
```

### Database Queries
```sql
-- Get all drug-like compounds
SELECT c.smiles, c.common_name, mp.molecular_weight, dp.lipinski_violations
FROM chemicals c
JOIN molecular_properties mp ON c.id = mp.chemical_id
JOIN drug_properties dp ON c.id = dp.chemical_id
WHERE dp.drug_like = 1;

-- Search by molecular weight range
SELECT * FROM chemical_full_profile 
WHERE molecular_weight BETWEEN 200 AND 500;
```

## Maintenance

### Regular Tasks
1. **Monitor log files** in `/logs` directory
2. **Clean old cache entries** (optional):
   ```sql
   DELETE FROM analysis_cache 
   WHERE analysis_date < DATE_SUB(NOW(), INTERVAL 30 DAY);
   ```
3. **Backup database** regularly
4. **Update statistics** in database manager

### Performance Optimization
- Indexes are automatically created for fast lookups
- Consider partitioning for very large datasets (>1M compounds)
- Regular OPTIMIZE TABLE commands for large tables

## Troubleshooting

### Common Issues

1. **"Database connection failed"**
   - Check MySQL is running in XAMPP
   - Verify credentials in .env file
   - Check database exists

2. **"Missing database tables"**
   - Re-run schema.sql import
   - Check for SQL errors in import process

3. **"Analysis caching not working"**
   - Verify `analysis_cache` table exists
   - Check storage permissions
   - Review error logs

### Support
- Check `/logs/error.log` for detailed error information
- Verify Python environment and RDKit installation
- Test individual components (database, Python analysis) separately

## Security Notes

- Database credentials are stored in `.env` file (not version controlled)
- All user inputs are sanitized before database queries
- Prepared statements used for all SQL operations
- Rate limiting applied to prevent abuse
- Access logs maintained for security monitoring
