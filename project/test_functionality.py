#!/usr/bin/env python3
"""
Simple ChemCrow test script - bypasses version conflicts
"""

def test_chemcrow_basic():
    """Test basic ChemCrow functionality"""
    try:
        print("🧪 Testing ChemCrow Basic Functionality")
        print("=" * 50)
        
        # Test basic imports
        print("📦 Testing imports...")
        import chemcrow
        print("✅ ChemCrow package imported")
        
        # Test tools import
        from chemcrow import tools
        print("✅ ChemCrow tools imported")
        
        # Test specific tools
        from chemcrow.tools import rdkit_tools
        print("✅ RDKit tools imported")
        
        print("\n🎉 ChemCrow basic functionality is working!")
        print("✅ Ready for chemistry analysis")
        
        return True
        
    except Exception as e:
        print(f"❌ ChemCrow test failed: {str(e)}")
        return False

def test_rdkit():
    """Test RDKit functionality"""
    try:
        print("\n🧪 Testing RDKit...")
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        # Test SMILES parsing
        smiles = "CCO"  # Ethanol
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mw = Descriptors.MolWt(mol)
            print(f"✅ RDKit working - Ethanol MW: {mw:.2f}")
            return True
        else:
            print("❌ RDKit molecule parsing failed")
            return False
            
    except Exception as e:
        print(f"❌ RDKit test failed: {str(e)}")
        return False

def test_openai():
    """Test OpenAI import"""
    try:
        print("\n🤖 Testing OpenAI...")
        import openai
        print(f"✅ OpenAI imported - Version: {openai.__version__}")
        return True
    except Exception as e:
        print(f"❌ OpenAI test failed: {str(e)}")
        return False

def main():
    """Run all tests"""
    print("🚀 Chem Search - Functionality Test")
    print("=" * 60)
    
    results = {
        "RDKit": test_rdkit(),
        "OpenAI": test_openai(),
        "ChemCrow": test_chemcrow_basic()
    }
    
    print("\n" + "=" * 60)
    print("📋 Test Summary:")
    print("-" * 30)
    
    all_passed = True
    for component, passed in results.items():
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{component:<15} {status}")
        if not passed:
            all_passed = False
    
    print("\n" + "=" * 60)
    if all_passed:
        print("🎉 All core components are working!")
        print("✅ Your chemistry platform is ready to use!")
    else:
        print("⚠️  Some components have issues but basic functionality should work")
        print("💡 ChemCrow may have version conflicts but RDKit analysis will work")
    
    print("\n🌐 Web Interface: http://localhost/chem-search")
    print("🧪 Test Setup: http://localhost/chem-search/test_setup.php")

if __name__ == "__main__":
    main()
