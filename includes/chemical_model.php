<?php
/**
 * Chemical Model - Handles all chemical data operations
 */

require_once __DIR__ . '/database.php';

class Chemical {
    private $db;
    
    public function __construct() {
        $this->db = Database::getInstance();
    }
    
    /**
     * Find chemical by SMILES or create new entry
     */
    public function findOrCreateBySmiles($smiles) {
        // First try to find by canonical SMILES
        $chemical = $this->findBySmiles($smiles);
        
        if ($chemical) {
            return $chemical;
        }
        
        // Create new chemical entry
        return $this->createChemical(['smiles' => $smiles]);
    }
    
    /**
     * Find chemical by SMILES notation
     */
    public function findBySmiles($smiles) {
        $sql = "SELECT * FROM chemicals WHERE smiles = :smiles1 OR smiles_canonical = :smiles2 LIMIT 1";
        return $this->db->fetchOne($sql, ['smiles1' => $smiles, 'smiles2' => $smiles]);
    }
    
    /**
     * Find chemical by InChI Key
     */
    public function findByInChIKey($inchiKey) {
        $sql = "SELECT * FROM chemicals WHERE inchi_key = :inchi_key LIMIT 1";
        return $this->db->fetchOne($sql, ['inchi_key' => $inchiKey]);
    }
    
    /**
     * Find chemical by CAS number
     */
    public function findByCAS($casNumber) {
        $sql = "SELECT * FROM chemicals WHERE cas_number = :cas_number LIMIT 1";
        return $this->db->fetchOne($sql, ['cas_number' => $casNumber]);
    }
    
    /**
     * Search chemicals by name or synonym
     */
    public function searchByName($name) {
        $searchTerm = '%' . $name . '%';
        $sql = "
            SELECT DISTINCT c.* FROM chemicals c
            LEFT JOIN chemical_synonyms cs ON c.id = cs.chemical_id
            WHERE c.iupac_name LIKE :name 
               OR c.common_name LIKE :name
               OR cs.synonym LIKE :name
            ORDER BY 
                CASE 
                    WHEN c.common_name LIKE :name THEN 1
                    WHEN c.iupac_name LIKE :name THEN 2
                    ELSE 3
                END
            LIMIT 20
        ";
        return $this->db->fetchAll($sql, ['name' => $searchTerm]);
    }
    
    /**
     * Create new chemical entry
     */
    public function createChemical($data) {
        $chemicalId = $this->db->insert('chemicals', $data);
        return $this->getById($chemicalId);
    }
    
    /**
     * Get chemical by ID
     */
    public function getById($id) {
        $sql = "SELECT * FROM chemicals WHERE id = :id";
        return $this->db->fetchOne($sql, ['id' => $id]);
    }
    
    /**
     * Get full chemical profile with all properties
     */
    public function getFullProfile($id) {
        $sql = "SELECT * FROM chemical_full_profile WHERE id = :id";
        return $this->db->fetchOne($sql, ['id' => $id]);
    }
    
    /**
     * Update chemical basic information
     */
    public function updateChemical($id, $data) {
        return $this->db->update('chemicals', $data, 'id = :id', ['id' => $id]);
    }
    
    /**
     * Store molecular properties
     */
    public function storeMolecularProperties($chemicalId, $properties) {
        // Check if properties already exist
        $existing = $this->db->fetchOne(
            "SELECT id FROM molecular_properties WHERE chemical_id = :id",
            ['id' => $chemicalId]
        );
        
        $properties['chemical_id'] = $chemicalId;
        
        if ($existing) {
            return $this->db->update(
                'molecular_properties', 
                $properties, 
                'chemical_id = :id', 
                ['id' => $chemicalId]
            );
        } else {
            return $this->db->insert('molecular_properties', $properties);
        }
    }
    
    /**
     * Store drug properties
     */
    public function storeDrugProperties($chemicalId, $properties) {
        $existing = $this->db->fetchOne(
            "SELECT id FROM drug_properties WHERE chemical_id = :id",
            ['id' => $chemicalId]
        );
        
        $properties['chemical_id'] = $chemicalId;
        
        if ($existing) {
            return $this->db->update(
                'drug_properties', 
                $properties, 
                'chemical_id = :id', 
                ['id' => $chemicalId]
            );
        } else {
            return $this->db->insert('drug_properties', $properties);
        }
    }
    
    /**
     * Store physical properties
     */
    public function storePhysicalProperties($chemicalId, $properties) {
        $existing = $this->db->fetchOne(
            "SELECT id FROM physical_properties WHERE chemical_id = :id",
            ['id' => $chemicalId]
        );
        
        $properties['chemical_id'] = $chemicalId;
        
        if ($existing) {
            return $this->db->update(
                'physical_properties', 
                $properties, 
                'chemical_id = :id', 
                ['id' => $chemicalId]
            );
        } else {
            return $this->db->insert('physical_properties', $properties);
        }
    }
    
    /**
     * Store safety properties
     */
    public function storeSafetyProperties($chemicalId, $properties) {
        $existing = $this->db->fetchOne(
            "SELECT id FROM safety_properties WHERE chemical_id = :id",
            ['id' => $chemicalId]
        );
        
        $properties['chemical_id'] = $chemicalId;
        
        if ($existing) {
            return $this->db->update(
                'safety_properties', 
                $properties, 
                'chemical_id = :id', 
                ['id' => $chemicalId]
            );
        } else {
            return $this->db->insert('safety_properties', $properties);
        }
    }
    
    /**
     * Add synonym for chemical
     */
    public function addSynonym($chemicalId, $synonym, $type = 'common', $source = null) {
        // Check if synonym already exists
        $existing = $this->db->fetchOne(
            "SELECT id FROM chemical_synonyms WHERE chemical_id = :id AND synonym = :synonym",
            ['id' => $chemicalId, 'synonym' => $synonym]
        );
        
        if (!$existing) {
            return $this->db->insert('chemical_synonyms', [
                'chemical_id' => $chemicalId,
                'synonym' => $synonym,
                'synonym_type' => $type,
                'source' => $source
            ]);
        }
        
        return $existing['id'];
    }
    
    /**
     * Cache analysis data
     */
    public function cacheAnalysis($chemicalId, $analysisType, $data, $source = 'rdkit', $version = '1.0') {
        return $this->db->insert('analysis_cache', [
            'chemical_id' => $chemicalId,
            'analysis_type' => $analysisType,
            'raw_data' => json_encode($data),
            'source' => $source,
            'version' => $version
        ]);
    }
    
    /**
     * Get cached analysis data
     */
    public function getCachedAnalysis($chemicalId, $analysisType, $maxAge = 86400) {
        $sql = "
            SELECT * FROM analysis_cache 
            WHERE chemical_id = :id 
              AND analysis_type = :type 
              AND analysis_date > DATE_SUB(NOW(), INTERVAL :max_age SECOND)
            ORDER BY analysis_date DESC 
            LIMIT 1
        ";
        
        $result = $this->db->fetchOne($sql, [
            'id' => $chemicalId,
            'type' => $analysisType,
            'max_age' => $maxAge
        ]);
        
        if ($result) {
            $result['raw_data'] = json_decode($result['raw_data'], true);
        }
        
        return $result;
    }
    
    /**
     * Log search activity
     */
    public function logSearch($sessionId, $chemicalId, $query, $searchType, $ipAddress, $userAgent) {
        return $this->db->insert('search_history', [
            'session_id' => $sessionId,
            'chemical_id' => $chemicalId,
            'search_query' => $query,
            'search_type' => $searchType,
            'ip_address' => $ipAddress,
            'user_agent' => $userAgent
        ]);
    }
    
    /**
     * Get popular searches
     */
    public function getPopularSearches($limit = 10, $days = 30) {
        $sql = "
            SELECT 
                search_query,
                search_type,
                COUNT(*) as search_count,
                MAX(search_date) as last_search
            FROM search_history 
            WHERE search_date > DATE_SUB(NOW(), INTERVAL :days DAY)
              AND chemical_id IS NOT NULL
            GROUP BY search_query, search_type
            ORDER BY search_count DESC, last_search DESC
            LIMIT :limit
        ";
        
        return $this->db->fetchAll($sql, ['days' => $days, 'limit' => $limit]);
    }
    
    /**
     * Get recent analyses
     */
    public function getRecentAnalyses($limit = 20) {
        $sql = "
            SELECT c.*, sh.search_date 
            FROM chemicals c
            JOIN search_history sh ON c.id = sh.chemical_id
            WHERE sh.search_date > DATE_SUB(NOW(), INTERVAL 7 DAY)
            ORDER BY sh.search_date DESC
            LIMIT :limit
        ";
        
        return $this->db->fetchAll($sql, ['limit' => $limit]);
    }
}
?>
