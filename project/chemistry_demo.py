"""
Chem Local Demo (No API Required)

This demo shows some of Chem' chemistry capabilities that work locally
without requiring API calls.
"""

def demo_rdkit_calculations():
    """Demonstrate local chemistry calculations using RDKit."""
    print("🧪 RDKit Chemistry Calculations Demo")
    print("-" * 40)
    
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    
    # Common molecules to analyze
    molecules = {
        "Water": "O",
        "Methane": "C",
        "Ethanol": "CCO",
        "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "Tylenol (Acetaminophen)": "CC(=O)NC1=CC=C(C=C1)O"
    }
    
    for name, smiles in molecules.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = rdMolDescriptors.CalcNumHBD(mol)
            hba = rdMolDescriptors.CalcNumHBA(mol)
            
            print(f"\n{name}:")
            print(f"  SMILES: {smiles}")
            print(f"  Molecular Weight: {mw:.2f} g/mol")
            print(f"  LogP: {logp:.2f}")
            print(f"  H-Bond Donors: {hbd}")
            print(f"  H-Bond Acceptors: {hba}")

def demo_molecule_conversion():
    """Demonstrate molecule format conversions."""
    print("\n\n🔄 Molecule Format Conversion Demo")
    print("-" * 40)
    
    from rdkit import Chem
    
    # Example: Convert caffeine from SMILES to different formats
    caffeine_smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    mol = Chem.MolFromSmiles(caffeine_smiles)
    
    if mol is not None:
        print(f"Caffeine conversions:")
        print(f"  SMILES: {caffeine_smiles}")
        print(f"  InChI: {Chem.MolToInchi(mol)}")
        print(f"  InChI Key: {Chem.MolToInchiKey(mol)}")
        
        # Add hydrogens for molecular formula
        mol_with_h = Chem.AddHs(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol_with_h)
        print(f"  Molecular Formula: {formula}")

def main():
    print("=" * 50)
    print("Chem Local Chemistry Demo")
    print("=" * 50)
    
    demo_rdkit_calculations()
    demo_molecule_conversion()
    
    print("\n" + "=" * 50)
    print("🎯 Local chemistry calculations completed!")
    print("\nFor AI-powered chemistry questions, you'll need:")
    print("1. Valid OpenAI API key with available credits")
    print("2. Set: $env:OPENAI_API_KEY='your-api-key'")
    print("3. Then use: python example_usage.py")
    print("=" * 50)

if __name__ == "__main__":
    main()
