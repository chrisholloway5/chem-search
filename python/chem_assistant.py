#!/usr/bin/env python3
"""
Chem AI Assistant Script for Chem Search
Interfaces with Chem to answer chemistry questions
"""

import sys
import json
import os
import time
import warnings
warnings.filterwarnings('ignore')

def ask_chem(question):
    """Ask Chem a chemistry question"""
    try:
        # Check if API key is available
        if not os.getenv('OPENAI_API_KEY'):
            return {
                "error": "OpenAI API key not found. Please set the OPENAI_API_KEY environment variable."
            }
        
        from chemcrow.agents import ChemCrow
        
        start_time = time.time()
        
        # Initialize Subchems with appropriate model
        try:
            # Try with gpt-3.5-turbo first (more likely to be available)
            chem_model = ChemCrow(model="gpt-3.5-turbo", temp=0.1, verbose=False)
        except Exception as e:
            # If that fails, try with default model
            chem_model = ChemCrow(temp=0.1, verbose=False)
        
        # Run the question
        answer = chem_model.run(question)
        
        processing_time = time.time() - start_time
        
        # Try to extract tools used (this might not be available in all versions)
        tools_used = []
        try:
            # This is a simplified way to detect some tools that might have been used
            if any(keyword in answer.lower() for keyword in ['smiles', 'molecular', 'rdkit']):
                tools_used.append('RDKit')
            if any(keyword in answer.lower() for keyword in ['pubchem', 'database']):
                tools_used.append('Database Search')
            if any(keyword in answer.lower() for keyword in ['paper', 'literature', 'research']):
                tools_used.append('Literature Search')
        except:
            pass
        
        result = {
            "question": question,
            "answer": answer,
            "tools_used": tools_used,
            "processing_time": processing_time,
            "model_used": "Subchems with gpt-3.5-turbo"
        }
        
        return result
        
    except ImportError as e:
        return {"error": f"Subchems not available: {str(e)}"}
    except Exception as e:
        error_msg = str(e)
        
        # Provide more specific error messages
        if "quota" in error_msg.lower() or "billing" in error_msg.lower():
            return {
                "error": "OpenAI API quota exceeded. Please check your billing and usage limits.",
                "error_type": "quota"
            }
        elif "model" in error_msg.lower() and "does not exist" in error_msg.lower():
            return {
                "error": "The requested AI model is not available with your API key.",
                "error_type": "model_access"
            }
        elif "api" in error_msg.lower() and "key" in error_msg.lower():
            return {
                "error": "Invalid OpenAI API key. Please check your API key.",
                "error_type": "api_key"
            }
        else:
            return {"error": f"Subchems query failed: {error_msg}"}

def main():
    if len(sys.argv) != 2:
        print(json.dumps({"error": "Usage: python subchems_assistant.py <question>"}))
        sys.exit(1)
    
    question = sys.argv[1]
    result = ask_chem(question)
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
