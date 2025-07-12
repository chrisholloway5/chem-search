# Chem-search

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![PHP](https://img.shields.io/badge/PHP-7.4%2B-blue.svg)](https://www.php.net/)
[![Python](https://img.shields.io/badge/Python-3.8%2B-green.svg)](https://www.python.org/)
[![ChemCrow](https://img.shields.io/badge/ChemCrow-Powered-orange.svg)](https://github.com/ur-whitelab/chemcrow-public)

> 🧪 An AI-powered chemistry platform for molecular analysis, property calculations, and chemical research

Chem-search is a comprehensive web-based chemistry platform that combines the power of ChemCrow AI agents with an intuitive web interface. Built with PHP and Python, it provides tools for molecular analysis, format conversion, database searches, and AI-assisted chemical research.

![Chem-search Interface](screenshot.png)

## ✨ Features

- **🔬 Molecular Analysis**: Calculate molecular properties, descriptors, and chemical characteristics
- **🤖 AI Chemistry Assistant**: Ask complex chemistry questions powered by ChemCrow's AI agent system
- **🔄 Format Converter**: Convert between SMILES, InChI, molecular formulas, and other chemical formats
- **🔍 Database Search**: Search chemical databases like PubChem for compound information
- **📊 Structure Viewer**: Visualize molecular structures in 2D and 3D formats
- **📚 Research Papers**: Find and analyze chemistry research papers using AI-powered search
- **📱 Responsive Design**: Mobile-friendly interface with Bootstrap 5

## 🚀 Quick Start

### Prerequisites

- **PHP 7.4+** with extensions: `curl`, `json`, `mbstring`
- **Python 3.8+** with pip
- **XAMPP** (recommended for Windows) or other web server
- **OpenAI API Key** for AI functionality

### Installation

1. **Set up the project**
   
   **Option A: Clone directly to htdocs (traditional)**
   ```bash
   cd C:\xampp\htdocs
   git clone https://github.com/chrisholloway5/Chem-search.git chem-search
   cd chem-search
   ```
   
   **Option B: Use shortcut (recommended)**
   ```bash
   # Clone to your preferred location
   git clone https://github.com/chrisholloway5/Chem-search.git
   cd Chem-search
   

2. **Set up Python environment and install chemistry packages**
   ```bash
   cd project
   python -m venv .venv
   
   # Activate virtual environment
   # Windows
   .venv\Scripts\activate
   
   # Linux/Mac
   source .venv/bin/activate
   
   # Install required chemistry packages
   pip install --upgrade pip
   pip install -r requirements.txt
   
   # Or install manually:
   pip install rdkit chemcrow openai langchain python-dotenv requests
   ```
   
   **Required Chemistry Packages:**
   - **RDKit**: Chemistry toolkit for molecular analysis
   - **ChemCrow**: AI agent framework for chemistry
   - **OpenAI**: API for GPT models
   - **LangChain**: Framework for AI agent workflows
   
   **Installation troubleshooting:**
   ```bash
   # If RDKit installation fails, try conda:
   conda install -c conda-forge rdkit
   
   # For Windows users with compilation issues:
   pip install rdkit-pypi  # Alternative RDKit package
   ```

3. **Configure environment variables**
   ```bash
   cp .env.example .env
   ```
   
   Edit `.env` with your settings:
   ```env
   # API Keys
   OPENAI_API_KEY=your-openai-api-key-here
   SERP_API_KEY=your-serpapi-key-here  # Optional
   
   # Application Settings
   APP_NAME=Chem Search
   APP_VERSION=1.0.0
   APP_DEBUG=false
   APP_URL=http://localhost/chem-search
   
   # Python Configuration
   PYTHON_PATH=C:/xampp/htdocs/chem-search/project/.venv/Scripts/python.exe
   
   # AI Model Settings
   DEFAULT_AI_MODEL=gpt-3.5-turbo
   DEFAULT_TEMPERATURE=0.1
   ```

4. **Set up database (optional)**
   ```bash
   cd database
   # Edit schema.sql if needed, then import
   mysql -u username -p database_name < schema.sql
   ```

5. **Start XAMPP and access the application**
   ```bash
   # Start XAMPP Control Panel
   # Enable Apache and MySQL services
   
   # Access your application at:
   # http://localhost/chem-search
   ```

6. **Verify installation**
   - Open `http://localhost/chem-search` in your browser
   - Check setup status at: `http://localhost/chem-search/test_setup.php`
   - Try the molecule analyzer with caffeine: `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`

## 🧪 Python Requirements & Verification

### Required Chemistry Packages

| Package | Purpose | Status Check |
|---------|---------|--------------|
| **RDKit** | Chemistry toolkit for molecular analysis | `python -c "import rdkit; print('✅ RDKit:', rdkit.__version__)"` |
| **ChemCrow** | AI agent framework for chemistry | `python -c "import chemcrow; print('✅ ChemCrow:', chemcrow.__version__)"` |
| **OpenAI** | API for GPT models | `python -c "import openai; print('✅ OpenAI:', openai.__version__)"` |
| **LangChain** | Framework for AI agent workflows | `python -c "import langchain; print('✅ LangChain:', langchain.__version__)"` |

### Installation Commands

```bash
# Navigate to project directory
cd project

# Activate virtual environment
.venv\Scripts\activate  # Windows
source .venv/bin/activate  # Linux/Mac

# Install all required packages
pip install -r requirements.txt

# Or install core packages individually
pip install rdkit chemcrow openai langchain python-dotenv requests
```

### Verify Installation

```bash
# Run verification script
python -c "
import sys
packages = ['rdkit', 'chemcrow', 'openai', 'langchain', 'requests', 'dotenv']
for pkg in packages:
    try:
        __import__(pkg)
        print(f'✅ {pkg} - OK')
    except ImportError:
        print(f'❌ {pkg} - MISSING')
"

# Or check via web interface
# Visit: http://localhost/chem-search/test_setup.php
```

### Common Installation Issues

#### RDKit Installation Problems
```bash
# If standard installation fails:
pip install rdkit-pypi

# Or use conda (recommended for complex environments):
conda install -c conda-forge rdkit
```

#### ChemCrow Dependencies
```bash
# If ChemCrow installation is slow or fails:
pip install --no-cache-dir chemcrow

# Install dependencies first:
pip install numpy pandas langchain openai
pip install chemcrow
```

#### Windows-specific Issues
```bash
# For Visual C++ compilation errors:
# Install Microsoft C++ Build Tools
# Or use pre-compiled wheels:
pip install --only-binary=all rdkit chemcrow
```

## 🏗️ Project Structure

```
Chem-search/
├── assets/                 # Frontend assets
│   ├── css/
│   └── js/
├── database/              # Database schema and setup
│   ├── schema.sql
│   └── setup.py
├── includes/              # PHP includes and configuration
│   ├── config.php
│   ├── database.php
│   └── header.php
├── logs/                  # Application logs
├── project/               # Python/ChemCrow integration
│   └── chemcrow/         # ChemCrow library
├── python/               # Python scripts
│   ├── analyze_molecule.py
│   ├── chem_assistant.py
│   └── format_converter.py
├── tools/                # Chemistry tools (PHP)
│   ├── ai_assistant.php
│   ├── molecule_analyzer.php
│   ├── format_converter.php
│   ├── database_search.php
│   ├── structure_viewer.php
│   └── paper_search.php
├── index.php             # Main landing page
├── about.php             # About page
└── README.md
```

## 🛠️ Technology Stack

### Backend
- **PHP 7.4+**: Web application framework
- **Python 3.8+**: ChemCrow integration and molecular analysis
- **ChemCrow**: AI chemistry agent system
- **RDKit**: Cheminformatics toolkit
- **OpenAI API**: Large language models for AI assistance

### Frontend
- **Bootstrap 5**: Responsive UI framework
- **JavaScript (ES6+)**: Interactive functionality
- **Font Awesome 6.4.0**: Icons and visual elements (CDN)
- **CSS3**: Custom styling and animations

### Integrations
- **LangChain**: Framework for AI agent workflows
- **PubChem**: Chemical database integration
- **Paper-QA**: Research paper analysis
- **Chemical databases**: Various chemistry data sources

## 📖 Usage Examples

### Molecular Analysis
```php
// Via web interface: tools/molecule_analyzer.php
// Input: SMILES notation
// Output: Molecular properties, descriptors, and visualizations
```

### AI Chemistry Assistant
```python
# Python integration
from chemcrow.agents import ChemCrow

chem_model = ChemCrow(model="gpt-4-0613", temp=0.1)
result = chem_model.run("What is the molecular weight of aspirin?")
```

### Format Conversion
```php
// Convert between chemical formats
// SMILES → InChI → Molecular Formula
// Via tools/format_converter.php
```

## 🔧 Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `OPENAI_API_KEY` | OpenAI API key for AI functionality | Required |
| `SERP_API_KEY` | SerpAPI key for web searches | Optional |
| `PYTHON_PATH` | Path to Python executable | Auto-detect |
| `DEFAULT_AI_MODEL` | OpenAI model to use | `gpt-3.5-turbo` |
| `DEFAULT_TEMPERATURE` | AI model temperature | `0.1` |
| `APP_DEBUG` | Enable debug mode | `false` |

### Database Configuration

```php
// database/schema.sql contains the database structure
// Optional: Used for storing search history and results
```

## 🧪 Available Tools

1. **Molecule Analyzer** (`tools/molecule_analyzer.php`)
   - Calculate molecular properties
   - Generate molecular descriptors
   - Visualize structures

2. **AI Assistant** (`tools/ai_assistant.php`)
   - Ask chemistry questions
   - Get synthesis suggestions
   - Analyze chemical reactions

3. **Format Converter** (`tools/format_converter.php`)
   - Convert between chemical formats
   - Validate chemical structures
   - Generate canonical representations

4. **Database Search** (`tools/database_search.php`)
   - Search PubChem database
   - Find chemical information
   - Retrieve compound data

5. **Structure Viewer** (`tools/structure_viewer.php`)
   - 2D/3D molecular visualization
   - Interactive structure display
   - Chemical drawing tools

6. **Paper Search** (`tools/paper_search.php`)
   - AI-powered literature search
   - Research paper analysis
   - Chemical knowledge extraction

## 🤝 Contributing

We welcome contributions to Chem-search! Here's how you can help:

1. **Fork the repository**
2. **Create a feature branch** (`git checkout -b feature/amazing-feature`)
3. **Commit your changes** (`git commit -m 'Add amazing feature'`)
4. **Push to the branch** (`git push origin feature/amazing-feature`)
5. **Open a Pull Request**

### Development Guidelines

- Follow PSR-12 coding standards for PHP
- Use PEP 8 for Python code
- Add tests for new features
- Update documentation as needed
- Ensure mobile responsiveness

## 🐛 Troubleshooting

### Common Issues

1. **Python module not found**
   ```bash
   # Ensure ChemCrow is installed
   pip install chemcrow rdkit
   ```

2. **OpenAI API errors**
   ```bash
   # Check your API key in .env
   export OPENAI_API_KEY=your-key-here
   ```

3. **Permission errors**
   ```bash
   # Ensure proper file permissions
   chmod 755 python/
   chmod +x python/*.py
   ```

4. **Database connection issues**
   ```bash
   # Check database configuration in includes/config.php
   # Ensure database user has proper permissions
   ```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **ChemCrow Team**: For the amazing AI chemistry agent framework
- **RDKit Community**: For the comprehensive cheminformatics toolkit
- **OpenAI**: For providing powerful language models
- **Bootstrap Team**: For the excellent UI framework

## 📚 References

### ChemCrow Paper
```bibtex
@article{bran2023chemcrow,
    title={ChemCrow: Augmenting large-language models with chemistry tools},
    author={Andres M Bran and Sam Cox and Oliver Schilter and Carlo Baldassari and Andrew D White and Philippe Schwaller},
    year={2023},
    eprint={2304.05376},
    archivePrefix={arXiv},
    primaryClass={physics.chem-ph},
    publisher={arXiv}
}
```

## 📞 Support

- **Issues**: [GitHub Issues](https://github.com/chrisholloway5/Chem-search/issues)
- **Discussions**: [GitHub Discussions](https://github.com/chrisholloway5/Chem-search/discussions)
- **Email**: [chris@example.com](mailto:chris@example.com)

---

<div align="center">
  <p><strong>Built with ❤️ for the chemistry community</strong></p>
  <p>
    <a href="https://github.com/chrisholloway5/Chem-search">⭐ Star this repository</a> |
    <a href="https://github.com/chrisholloway5/Chem-search/fork">🍴 Fork it</a> |
    <a href="https://github.com/chrisholloway5/Chem-search/issues">🐛 Report bugs</a>
  </p>
</div>
