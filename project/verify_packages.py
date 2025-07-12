#!/usr/bin/env python3
"""
Chemistry Package Verification Script
Checks if all required packages for Chem Search are installed and working.
"""

import sys
import importlib
from typing import Dict, List, Tuple

def check_package(package_name: str, import_name: str = None) -> Tuple[bool, str, str]:
    """
    Check if a package is installed and get its version.
    
    Args:
        package_name: Name of the package for display
        import_name: Name to use for import (if different from package_name)
    
    Returns:
        Tuple of (is_installed, version, error_message)
    """
    if import_name is None:
        import_name = package_name.lower().replace('-', '_')
    
    try:
        module = importlib.import_module(import_name)
        
        # Try to get version
        version = "Unknown"
        for attr in ['__version__', 'version', 'VERSION']:
            if hasattr(module, attr):
                version = getattr(module, attr)
                break
        
        return True, str(version), ""
    
    except ImportError as e:
        return False, "", str(e)
    except Exception as e:
        return False, "", f"Error loading {package_name}: {str(e)}"

def main():
    """Main verification function."""
    print("🧪 Chem Search - Python Package Verification")
    print("=" * 50)
    
    # Core chemistry packages
    core_packages = [
        ("RDKit", "rdkit"),
        ("ChemCrow", "chemcrow"),
        ("OpenAI", "openai"),
        ("LangChain", "langchain"),
    ]
    
    # Supporting packages
    support_packages = [
        ("Python-dotenv", "dotenv"),
        ("Requests", "requests"),
        ("NumPy", "numpy"),
        ("Pandas", "pandas"),
        ("Tiktoken", "tiktoken"),
    ]
    
    # Optional packages
    optional_packages = [
        ("Streamlit", "streamlit"),
        ("IPython", "IPython"),
        ("Paper-QA", "paperqa"),
        ("Molbloom", "molbloom"),
    ]
    
    all_good = True
    
    def check_category(packages: List[Tuple[str, str]], category: str, required: bool = True):
        nonlocal all_good
        print(f"\n📦 {category} Packages:")
        print("-" * 30)
        
        category_good = True
        for pkg_name, import_name in packages:
            is_installed, version, error = check_package(pkg_name, import_name)
            
            if is_installed:
                print(f"✅ {pkg_name:<15} v{version}")
            else:
                status = "❌" if required else "⚠️"
                print(f"{status} {pkg_name:<15} MISSING")
                if required:
                    category_good = False
                    print(f"   Error: {error}")
        
        if required and not category_good:
            all_good = False
        
        return category_good
    
    # Check all categories
    core_good = check_category(core_packages, "Core Chemistry", required=True)
    support_good = check_category(support_packages, "Supporting", required=True)
    check_category(optional_packages, "Optional", required=False)
    
    # Final summary
    print("\n" + "=" * 50)
    print("📋 Installation Summary:")
    print("-" * 25)
    
    if all_good:
        print("🎉 All required packages are installed!")
        print("✅ Your Chem Search environment is ready to use.")
        print("\n🚀 Next steps:")
        print("   1. Start XAMPP (Apache + MySQL)")
        print("   2. Visit: http://localhost/chem-search")
        print("   3. Add your OpenAI API key to .env file")
    else:
        print("❌ Some required packages are missing.")
        print("\n🔧 To fix missing packages:")
        print("   pip install -r requirements.txt")
        print("\n   Or install individually:")
        if not core_good:
            print("   pip install rdkit chemcrow openai langchain")
        if not support_good:
            print("   pip install python-dotenv requests numpy pandas")
    
    print("\n📖 For detailed setup instructions:")
    print("   See README.md or XAMPP_SETUP.md")
    
    return 0 if all_good else 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
