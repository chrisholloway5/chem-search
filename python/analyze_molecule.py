#!/usr/bin/env python3
"""
Enhanced Molecule Analysis Script for Chem Search
Computes comprehensive molecular properties and stores in database
"""

import sys
import json
import traceback
from typing import Dict, Any, Optional

try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski, QED
    from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds, CalcTPSA
    from rdkit.Chem.GraphDescriptors import BertzCT
    import rdkit
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

def calculate_comprehensive_properties(smiles: str) -> Dict[str, Any]:
    """Calculate comprehensive molecular properties from SMILES"""
    
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is not available")
    
    try:
        # Parse the molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")
        
        # Get canonical SMILES
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        
        # Basic identifiers
        properties = {
            'smiles_canonical': canonical_smiles,
            'molecular_formula': rdMolDescriptors.CalcMolFormula(mol),
        }
        
        # Try to generate InChI and InChI Key
        try:
            inchi = Chem.MolToInchi(mol)
            inchi_key = Chem.MolToInchiKey(mol)
            properties.update({
                'inchi': inchi,
                'inchi_key': inchi_key
            })
        except:
            properties.update({
                'inchi': None,
                'inchi_key': None
            })
        
        # Molecular properties
        molecular_props = {
            'molecular_weight': round(Descriptors.MolWt(mol), 4),
            'exact_mass': round(Descriptors.ExactMolWt(mol), 8),
            'monoisotopic_mass': round(Descriptors.ExactMolWt(mol), 8),
            
            # Lipophilicity
            'xlogp3': round(Descriptors.MolLogP(mol), 4),
            'clogp': round(Crippen.MolLogP(mol), 4),
            
            # Hydrogen bonding
            'hbd_count': Descriptors.NumHDonors(mol),
            'hba_count': Descriptors.NumHAcceptors(mol),
            
            # Structural properties
            'rotatable_bond_count': CalcNumRotatableBonds(mol),
            'heavy_atom_count': mol.GetNumHeavyAtoms(),
            'aromatic_ring_count': rdMolDescriptors.CalcNumAromaticRings(mol),
            'aliphatic_ring_count': rdMolDescriptors.CalcNumAliphaticRings(mol),
            
            # Surface area and topology
            'tpsa': round(CalcTPSA(mol), 3),
            
            # Formal properties
            'formal_charge': Chem.GetFormalCharge(mol),
            
            # Complexity
            'bertz_complexity': round(BertzCT(mol), 3),
            
            # Stereochemistry
            'defined_atom_stereocenter_count': len(Chem.FindMolChiralCenters(mol, includeUnassigned=False)),
            'undefined_atom_stereocenter_count': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)) - len(Chem.FindMolChiralCenters(mol, includeUnassigned=False)),
            
            # Component information
            'covalently_bonded_unit_count': len(Chem.GetMolFrags(mol)),
            
            # Validation
            'compound_is_canonicalized': True,
            'has_stereochemistry': len(Chem.FindMolChiralCenters(mol, includeUnassigned=True)) > 0
        }
        
        # Additional descriptors
        try:
            molecular_props.update({
                'radical_electron_count': Descriptors.NumRadicalElectrons(mol),
                'isotope_atom_count': sum(1 for atom in mol.GetAtoms() if atom.GetIsotope() != 0),
                'total_charge': sum(atom.GetFormalCharge() for atom in mol.GetAtoms()),
            })
        except:
            pass
        
        # Drug-like properties (Lipinski's Rule of Five)
        drug_props = calculate_drug_properties(mol)
        
        # Physical properties estimates
        physical_props = calculate_physical_properties(mol)
        
        # Additional molecular descriptors
        additional_descriptors = calculate_additional_descriptors(mol)
        
        return {
            'basic_properties': properties,
            'molecular_properties': molecular_props,
            'drug_properties': drug_props,
            'physical_properties': physical_props,
            'additional_descriptors': additional_descriptors,
            'analysis_metadata': {
                'rdkit_version': rdkit.__version__,
                'canonical_smiles': canonical_smiles,
                'analysis_success': True
            }
        }
        
    except Exception as e:
        raise Exception(f"Error calculating properties: {str(e)}")

def calculate_drug_properties(mol) -> Dict[str, Any]:
    """Calculate drug-like properties and ADMET predictions"""
    
    # Lipinski's Rule of Five
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    lipinski_violations = 0
    lipinski_violations += 1 if mw > 500 else 0
    lipinski_violations += 1 if logp > 5 else 0
    lipinski_violations += 1 if hbd > 5 else 0
    lipinski_violations += 1 if hba > 10 else 0
    
    # Veber's Rule
    rotatable_bonds = CalcNumRotatableBonds(mol)
    tpsa = CalcTPSA(mol)
    
    veber_violations = 0
    veber_violations += 1 if rotatable_bonds > 10 else 0
    veber_violations += 1 if tpsa > 140 else 0
    
    # QED (Quantitative Estimate of Drug-likeness)
    try:
        qed_score = QED.qed(mol)
    except:
        qed_score = None
    
    return {
        'lipinski_violations': lipinski_violations,
        'lipinski_mw_violation': mw > 500,
        'lipinski_logp_violation': logp > 5,
        'lipinski_hbd_violation': hbd > 5,
        'lipinski_hba_violation': hba > 10,
        
        'veber_violations': veber_violations,
        'veber_rotatable_bonds_violation': rotatable_bonds > 10,
        'veber_tpsa_violation': tpsa > 140,
        
        'qed_score': round(qed_score, 4) if qed_score else None,
        
        # Basic ADMET predictions (simplified)
        'drug_like': lipinski_violations <= 1 and veber_violations == 0,
        'oral_bioavailability_likely': lipinski_violations <= 1 and rotatable_bonds <= 7,
    }

def calculate_physical_properties(mol) -> Dict[str, Any]:
    """Calculate estimated physical properties"""
    
    # Basic estimates - these would ideally come from experimental data or ML models
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    
    # Very basic estimates - replace with ML models in production
    estimated_props = {
        'estimated_logS': round(-0.5 * logp - 0.01 * mw + 0.5, 3),  # Simplified solubility estimate
        'estimated_melting_point': None,  # Would need ML model
        'estimated_boiling_point': None,  # Would need ML model
    }
    
    return estimated_props

def calculate_additional_descriptors(mol) -> Dict[str, Any]:
    """Calculate additional molecular descriptors"""
    
    try:
        additional = {
            # Kappa descriptors
            'kappa1': round(rdMolDescriptors.CalcKappa1(mol), 4),
            'kappa2': round(rdMolDescriptors.CalcKappa2(mol), 4),
            'kappa3': round(rdMolDescriptors.CalcKappa3(mol), 4),
            
            # Connectivity indices
            'chi0v': round(rdMolDescriptors.CalcChi0v(mol), 4),
            'chi1v': round(rdMolDescriptors.CalcChi1v(mol), 4),
            'chi2v': round(rdMolDescriptors.CalcChi2v(mol), 4),
            
            # Ring information
            'ring_count': rdMolDescriptors.CalcNumRings(mol),
            'saturated_ring_count': rdMolDescriptors.CalcNumSaturatedRings(mol),
            'heterocycle_count': rdMolDescriptors.CalcNumHeterocycles(mol),
            
            # Atom counts by type
            'carbon_count': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C'),
            'nitrogen_count': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N'),
            'oxygen_count': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O'),
            'sulfur_count': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'S'),
            'phosphorus_count': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'P'),
            'halogen_count': sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in ['F', 'Cl', 'Br', 'I']),
        }
        
        return additional
        
    except Exception as e:
        return {'error': f"Error calculating additional descriptors: {str(e)}"}

def main():
    if len(sys.argv) != 2:
        print(json.dumps({
            'error': 'Usage: python analyze_molecule.py <SMILES>',
            'success': False
        }))
        sys.exit(1)
    
    smiles = sys.argv[1].strip()
    
    if not smiles:
        print(json.dumps({
            'error': 'Empty SMILES provided',
            'success': False
        }))
        sys.exit(1)
    
    try:
        result = calculate_comprehensive_properties(smiles)
        result['success'] = True
        print(json.dumps(result, indent=2))
        
    except Exception as e:
        error_result = {
            'error': str(e),
            'success': False,
            'traceback': traceback.format_exc() if '--debug' in sys.argv else None
        }
        print(json.dumps(error_result))
        sys.exit(1)

if __name__ == '__main__':
    main()
