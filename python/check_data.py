#!/usr/bin/env python3
"""
Quick test to see what data structure is being returned
"""
import json
import sys
import os

# Add the current directory to path
sys.path.insert(0, os.path.dirname(__file__))

# Import and run the analysis
from analyze_molecule import calculate_comprehensive_properties

smiles = 'CCO'
result = calculate_comprehensive_properties(smiles)

print("=== PHYSICAL PROPERTIES ===")
print(json.dumps(result['physical_properties'], indent=2))

print("\n=== DRUG PROPERTIES ===")
print(json.dumps(result['drug_properties'], indent=2))
