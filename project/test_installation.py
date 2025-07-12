"""
Test Chem Installation (without API calls)

This script tests that Chem is properly installed and can be imported
without making any API calls.
"""

def test_imports():
    """Test that all Chem components can be imported."""
    print("Testing Chem installation...")
    
    try:
        import chemcrow
        print("✅ Chem base package imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import chemcrow: {e}")
        return False
    
    try:
        from chemcrow.agents import ChemCrow
        print("✅ Chem agents module imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import Chem agents: {e}")
        return False
    
    try:
        import chemcrow.tools
        print("✅ Chem tools imported successfully")
    except ImportError as e:
        print(f"❌ Failed to import Chem tools: {e}")
        return False
    
    return True

def test_rdkit():
    """Test that RDKit (chemistry toolkit) is working."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
        
        # Simple test: calculate molecular weight of water (H2O)
        mol = Chem.MolFromSmiles("O")
        if mol is not None:
            mw = Descriptors.MolWt(mol)
            print(f"✅ RDKit working - Water (H2O) molecular weight: {mw:.2f}")
            return True
        else:
            print("❌ RDKit failed to parse simple molecule")
            return False
    except ImportError as e:
        print(f"❌ Failed to import RDKit: {e}")
        return False

def main():
    print("=" * 50)
    print("Chem Installation Test")
    print("=" * 50)
    
    success = True
    success &= test_imports()
    success &= test_rdkit()
    
    print("\n" + "=" * 50)
    if success:
        print("🎉 All tests passed! Chem is properly installed.")
        print("\nTo use Chem with API calls, you need:")
        print("1. Valid OpenAI API key with available quota")
        print("2. Set environment variable: $env:OPENAI_API_KEY='your-key'")
    else:
        print("❌ Some tests failed. Please check the installation.")
    print("=" * 50)

if __name__ == "__main__":
    main()
