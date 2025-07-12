#!/usr/bin/env python3
"""
Database Setup Script for Subchems Search
Sets up the MySQL database and tables
"""

import mysql.connector
import sys
import os

def create_database():
    """Create the database and tables"""
    
    # Database configuration
    config = {
        'host': 'localhost',
        'user': 'root',
        'password': '',  # Update with your MySQL password
        'charset': 'utf8mb4',
        'collation': 'utf8mb4_unicode_ci'
    }
    
    try:
        # Connect to MySQL server
        connection = mysql.connector.connect(**config)
        cursor = connection.cursor()
        
        print("✓ Connected to MySQL server")
        
        # Create database
        cursor.execute("CREATE DATABASE IF NOT EXISTS subchems_search CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci")
        print("✓ Database 'subchems_search' created/verified")
        
        # Use the database
        cursor.execute("USE subchems_search")
        
        # Read and execute schema file
        schema_file = os.path.join(os.path.dirname(__file__), '..', 'database', 'schema.sql')
        
        if not os.path.exists(schema_file):
            print("❌ Schema file not found:", schema_file)
            return False
        
        with open(schema_file, 'r', encoding='utf-8') as f:
            schema_sql = f.read()
        
        # Split and execute SQL statements
        statements = schema_sql.split(';')
        
        for statement in statements:
            statement = statement.strip()
            if statement and not statement.startswith('--'):
                try:
                    cursor.execute(statement)
                except mysql.connector.Error as e:
                    if "already exists" not in str(e):
                        print(f"Warning: {e}")
        
        connection.commit()
        print("✓ Database schema created successfully")
        
        # Verify tables
        cursor.execute("SHOW TABLES")
        tables = cursor.fetchall()
        
        expected_tables = [
            'chemicals', 'molecular_properties', 'drug_properties',
            'physical_properties', 'safety_properties', 'chemical_synonyms',
            'analysis_cache', 'search_history'
        ]
        
        existing_tables = [table[0] for table in tables]
        
        print(f"\n📋 Tables created: {len(existing_tables)}")
        for table in existing_tables:
            print(f"   ✓ {table}")
        
        missing_tables = set(expected_tables) - set(existing_tables)
        if missing_tables:
            print(f"\n❌ Missing tables: {', '.join(missing_tables)}")
            return False
        
        print("\n✅ Database setup completed successfully!")
        print(f"Database: subchems_search")
        print(f"Tables: {len(existing_tables)}")
        print(f"Host: {config['host']}")
        
        return True
        
    except mysql.connector.Error as e:
        print(f"❌ MySQL Error: {e}")
        return False
    except Exception as e:
        print(f"❌ Error: {e}")
        return False
    finally:
        if 'connection' in locals() and connection.is_connected():
            cursor.close()
            connection.close()

def main():
    print("🔧 Setting up Subchems Search Database...")
    print("=" * 50)
    
    success = create_database()
    
    if success:
        print("\n🎉 Setup completed! You can now:")
        print("1. Update your .env file with database credentials")
        print("2. Test the molecule analyzer with database storage")
        print("3. View stored chemical profiles in the database")
        sys.exit(0)
    else:
        print("\n❌ Setup failed. Please check the errors above.")
        sys.exit(1)

if __name__ == '__main__':
    main()
