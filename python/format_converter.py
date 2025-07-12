#!/usr/bin/env python3
"""
Chemical Format Converter Script for Chem Search
Converts between different chemical structure formats
"""

import sys
import json
import warnings
warnings.filterwarnings('ignore')

def convert_format(input_str, input_type, output_type):
    """Convert between chemical formats using RDKit"""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        
        # Parse input based on type
        mol = None
        if input_type.lower() == 'smiles':
            mol = Chem.MolFromSmiles(input_str.strip())
        elif input_type.lower() == 'inchi':
            mol = Chem.MolFromInchi(input_str.strip())
        elif input_type.lower() == 'mol':
            mol = Chem.MolFromMolBlock(input_str.strip())
        else:
            return {"error": f"Unsupported input format: {input_type}"}
        
        if mol is None:
            return {"error": f"Invalid {input_type.upper()} format or could not parse structure"}
        
        # Convert to output format
        output = None
        if output_type.lower() == 'smiles':
            output = Chem.MolToSmiles(mol)
        elif output_type.lower() == 'inchi':
            output = Chem.MolToInchi(mol)
        elif output_type.lower() == 'inchi_key':
            output = Chem.MolToInchiKey(mol)
        elif output_type.lower() == 'formula':
            mol_with_h = Chem.AddHs(mol)
            output = rdMolDescriptors.CalcMolFormula(mol_with_h)
        elif output_type.lower() == 'mol':
            output = Chem.MolToMolBlock(mol)
        else:
            return {"error": f"Unsupported output format: {output_type}"}
        
        if output is None:
            return {"error": f"Failed to convert to {output_type.upper()} format"}
        
        # Additional information
        additional_info = {}
        try:
            additional_info['molecular_weight'] = f"{Chem.Descriptors.MolWt(mol):.2f} g/mol"
            additional_info['num_atoms'] = str(mol.GetNumAtoms())
            additional_info['num_bonds'] = str(mol.GetNumBonds())
            mol_with_h = Chem.AddHs(mol)
            additional_info['molecular_formula'] = rdMolDescriptors.CalcMolFormula(mol_with_h)
        except:
            pass
        
        result = {
            "input": input_str,
            "input_type": input_type,
            "output": output,
            "output_type": output_type,
            "additional_info": additional_info
        }
        
        return result
        
    except ImportError as e:
        return {"error": f"Missing required library: {str(e)}"}
    except Exception as e:
        return {"error": f"Conversion failed: {str(e)}"}

def main():
    if len(sys.argv) != 4:
        print(json.dumps({"error": "Usage: python format_converter.py <input> <input_type> <output_type>"}))
        sys.exit(1)
    
    input_str = sys.argv[1]
    input_type = sys.argv[2]
    output_type = sys.argv[3]
    
    result = convert_format(input_str, input_type, output_type)
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
