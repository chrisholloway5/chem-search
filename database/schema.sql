-- Database schema for Subchems Search Chemical Profiles
-- MySQL/MariaDB SQL script

CREATE DATABASE IF NOT EXISTS subchems_search 
CHARACTER SET utf8mb4 
COLLATE utf8mb4_unicode_ci;

USE subchems_search;

-- Main chemicals table with identifiers and basic properties
CREATE TABLE chemicals (
    id INT AUTO_INCREMENT PRIMARY KEY,
    smiles VARCHAR(1000) NOT NULL,
    smiles_canonical VARCHAR(1000),
    iupac_name TEXT,
    common_name VARCHAR(500),
    inchi TEXT,
    inchi_key VARCHAR(255),
    molecular_formula VARCHAR(255),
    cas_number VARCHAR(50),
    pubchem_cid BIGINT,
    chemspider_id BIGINT,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    -- Indexes for fast lookups
    UNIQUE KEY idx_smiles_canonical (smiles_canonical),
    INDEX idx_inchi_key (inchi_key),
    INDEX idx_cas_number (cas_number),
    INDEX idx_pubchem_cid (pubchem_cid),
    INDEX idx_molecular_formula (molecular_formula)
);

-- Computed molecular properties
CREATE TABLE molecular_properties (
    id INT AUTO_INCREMENT PRIMARY KEY,
    chemical_id INT NOT NULL,
    
    -- Basic molecular properties
    molecular_weight DECIMAL(10,4),
    exact_mass DECIMAL(15,8),
    monoisotopic_mass DECIMAL(15,8),
    
    -- Lipophilicity and solubility
    xlogp3 DECIMAL(8,4),
    clogp DECIMAL(8,4),
    
    -- Hydrogen bonding
    hbd_count INT, -- Hydrogen Bond Donor Count
    hba_count INT, -- Hydrogen Bond Acceptor Count
    
    -- Structural properties
    rotatable_bond_count INT,
    heavy_atom_count INT,
    aromatic_ring_count INT,
    aliphatic_ring_count INT,
    
    -- Surface area and volume
    tpsa DECIMAL(8,3), -- Topological Polar Surface Area
    molecular_volume DECIMAL(10,3),
    
    -- Formal properties
    formal_charge INT,
    total_charge INT,
    
    -- Complexity measures
    complexity_score DECIMAL(10,3),
    bertz_complexity DECIMAL(10,3),
    
    -- Atom counts
    isotope_atom_count INT,
    radical_electron_count INT,
    
    -- Stereochemistry
    defined_atom_stereocenter_count INT,
    undefined_atom_stereocenter_count INT,
    defined_bond_stereocenter_count INT,
    undefined_bond_stereocenter_count INT,
    
    -- Unit properties
    covalently_bonded_unit_count INT,
    component_count INT,
    
    -- Validation flags
    compound_is_canonicalized BOOLEAN DEFAULT FALSE,
    has_stereochemistry BOOLEAN DEFAULT FALSE,
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    FOREIGN KEY (chemical_id) REFERENCES chemicals(id) ON DELETE CASCADE,
    INDEX idx_molecular_weight (molecular_weight),
    INDEX idx_xlogp3 (xlogp3),
    INDEX idx_tpsa (tpsa)
);

-- Drug-like properties and ADMET predictions
CREATE TABLE drug_properties (
    id INT AUTO_INCREMENT PRIMARY KEY,
    chemical_id INT NOT NULL,
    
    -- Lipinski's Rule of Five
    lipinski_violations INT DEFAULT 0,
    lipinski_mw_violation BOOLEAN DEFAULT FALSE,
    lipinski_logp_violation BOOLEAN DEFAULT FALSE,
    lipinski_hbd_violation BOOLEAN DEFAULT FALSE,
    lipinski_hba_violation BOOLEAN DEFAULT FALSE,
    
    -- Veber's Rule
    veber_violations INT DEFAULT 0,
    veber_rotatable_bonds_violation BOOLEAN DEFAULT FALSE,
    veber_tpsa_violation BOOLEAN DEFAULT FALSE,
    
    -- Ghose filter
    ghose_violations INT DEFAULT 0,
    
    -- PAINS (Pan Assay Interference Compounds)
    pains_alerts INT DEFAULT 0,
    pains_description TEXT,
    
    -- Synthetic accessibility
    sa_score DECIMAL(4,2), -- 1 (easy) to 10 (difficult)
    
    -- Blood-brain barrier penetration
    bbb_penetration_probability DECIMAL(4,3),
    
    -- Human intestinal absorption
    hia_probability DECIMAL(4,3),
    
    -- Cytochrome P450 interactions
    cyp3a4_substrate_probability DECIMAL(4,3),
    cyp3a4_inhibitor_probability DECIMAL(4,3),
    cyp2d6_substrate_probability DECIMAL(4,3),
    cyp2d6_inhibitor_probability DECIMAL(4,3),
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    FOREIGN KEY (chemical_id) REFERENCES chemicals(id) ON DELETE CASCADE
);

-- Spectroscopic and physical properties
CREATE TABLE physical_properties (
    id INT AUTO_INCREMENT PRIMARY KEY,
    chemical_id INT NOT NULL,
    
    -- Physical state and appearance
    physical_state VARCHAR(50), -- solid, liquid, gas
    color VARCHAR(100),
    odor VARCHAR(200),
    
    -- Thermodynamic properties
    melting_point DECIMAL(8,2), -- Celsius
    boiling_point DECIMAL(8,2), -- Celsius
    flash_point DECIMAL(8,2), -- Celsius
    autoignition_temperature DECIMAL(8,2), -- Celsius
    
    -- Density and related
    density DECIMAL(8,4), -- g/cm³
    vapor_pressure DECIMAL(15,8), -- mmHg at 25°C
    vapor_density DECIMAL(8,3), -- relative to air
    
    -- Solubility
    water_solubility DECIMAL(15,8), -- mg/L
    water_solubility_unit VARCHAR(20) DEFAULT 'mg/L',
    logS DECIMAL(6,3), -- Log aqueous solubility
    
    -- Refractive index
    refractive_index DECIMAL(8,5),
    
    -- Stability
    stability VARCHAR(500),
    shelf_life VARCHAR(200),
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    FOREIGN KEY (chemical_id) REFERENCES chemicals(id) ON DELETE CASCADE
);

-- Safety and toxicology data
CREATE TABLE safety_properties (
    id INT AUTO_INCREMENT PRIMARY KEY,
    chemical_id INT NOT NULL,
    
    -- GHS Classification
    ghs_pictograms JSON,
    ghs_signal_word VARCHAR(50), -- Danger, Warning
    ghs_hazard_statements JSON,
    ghs_precautionary_statements JSON,
    
    -- Toxicity data
    oral_ld50_rat DECIMAL(12,3), -- mg/kg
    dermal_ld50_rabbit DECIMAL(12,3), -- mg/kg
    inhalation_lc50_rat DECIMAL(12,3), -- mg/L/4h
    
    -- Mutagenicity and carcinogenicity
    ames_test_result VARCHAR(20), -- positive, negative, inconclusive
    carcinogenicity_classification VARCHAR(100),
    
    -- Environmental data
    bioconcentration_factor DECIMAL(12,3),
    biodegradation_probability DECIMAL(4,3),
    
    -- Reactivity
    incompatible_materials TEXT,
    hazardous_decomposition TEXT,
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    
    FOREIGN KEY (chemical_id) REFERENCES chemicals(id) ON DELETE CASCADE
);

-- Synonyms and alternative names
CREATE TABLE chemical_synonyms (
    id INT AUTO_INCREMENT PRIMARY KEY,
    chemical_id INT NOT NULL,
    synonym VARCHAR(500) NOT NULL,
    synonym_type VARCHAR(50), -- IUPAC, common, trade, systematic
    source VARCHAR(100), -- PubChem, ChemSpider, etc.
    
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    FOREIGN KEY (chemical_id) REFERENCES chemicals(id) ON DELETE CASCADE,
    INDEX idx_synonym (synonym),
    INDEX idx_synonym_type (synonym_type)
);

-- Analysis history and caching
CREATE TABLE analysis_cache (
    id INT AUTO_INCREMENT PRIMARY KEY,
    chemical_id INT NOT NULL,
    analysis_type VARCHAR(100) NOT NULL, -- rdkit_properties, pubchem_lookup, etc.
    raw_data JSON,
    analysis_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    source VARCHAR(100), -- rdkit, pubchem, chemspider
    version VARCHAR(50), -- version of analysis tool/database
    
    FOREIGN KEY (chemical_id) REFERENCES chemicals(id) ON DELETE CASCADE,
    INDEX idx_analysis_type (analysis_type),
    INDEX idx_analysis_date (analysis_date)
);

-- User search and analysis history
CREATE TABLE search_history (
    id INT AUTO_INCREMENT PRIMARY KEY,
    session_id VARCHAR(255),
    chemical_id INT,
    search_query VARCHAR(1000),
    search_type VARCHAR(50), -- smiles, name, cas, etc.
    ip_address VARCHAR(45),
    user_agent TEXT,
    search_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    
    FOREIGN KEY (chemical_id) REFERENCES chemicals(id) ON DELETE SET NULL,
    INDEX idx_session_id (session_id),
    INDEX idx_search_date (search_date),
    INDEX idx_search_type (search_type)
);

-- Views for common queries
CREATE VIEW chemical_full_profile AS
SELECT 
    c.*,
    mp.molecular_weight,
    mp.exact_mass,
    mp.xlogp3,
    mp.hbd_count,
    mp.hba_count,
    mp.rotatable_bond_count,
    mp.tpsa,
    mp.heavy_atom_count,
    mp.formal_charge,
    mp.complexity_score,
    dp.lipinski_violations,
    dp.sa_score,
    pp.water_solubility,
    pp.melting_point,
    pp.boiling_point
FROM chemicals c
LEFT JOIN molecular_properties mp ON c.id = mp.chemical_id
LEFT JOIN drug_properties dp ON c.id = dp.chemical_id
LEFT JOIN physical_properties pp ON c.id = pp.chemical_id;
