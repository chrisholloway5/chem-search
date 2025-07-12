"""
Example usage of Chem

Before running this example, make sure to set your API keys:
- Set OPENAI_API_KEY environment variable with your OpenAI API key
- Optionally set SERP_API_KEY environment variable for search functionality

Example (PowerShell):
    $env:OPENAI_API_KEY="your-openai-api-key"
    $env:SERP_API_KEY="your-serpapi-api-key"
"""

import os
from chemcrow.agents import ChemCrow

def main():
    # Check if API key is set
    if not os.getenv("OPENAI_API_KEY"):
        print("Please set your OPENAI_API_KEY environment variable before using Chem.")
        print("You can do this by running:")
        print('$env:OPENAI_API_KEY="your-openai-api-key"')
        return
    
    # Initialize Chem with correct parameters
    print("Initializing Chem...")
    # Note: Use gpt-3.5-turbo as default since gpt-4 requires higher API access
    chem_model = ChemCrow(model="gpt-3.5-turbo", temp=0.1, verbose=True)
    
    # Example query
    question = "What is the molecular weight of tylenol?"
    print(f"Question: {question}")
    
    # Run the query
    try:
        result = chem_model.run(question)
        print(f"Answer: {result}")
    except Exception as e:
        print(f"Error: {e}")
        if "quota" in str(e).lower() or "billing" in str(e).lower():
            print("\n🚨 API Quota Issue:")
            print("Your OpenAI API key has exceeded its quota or has billing issues.")
            print("Please check your OpenAI account at https://platform.openai.com/")
            print("You may need to:")
            print("1. Add billing information")
            print("2. Purchase credits")
            print("3. Wait for quota reset")
        elif "does not exist" in str(e) or "access" in str(e):
            print("\n🚨 Model Access Issue:")
            print("Your API key doesn't have access to the requested model.")
            print("Try using 'gpt-3.5-turbo' instead of 'gpt-4' models.")
        else:
            print("Make sure your OpenAI API key is valid and you have sufficient credits.")

if __name__ == "__main__":
    main()
