<?php
/**
 * Database configuration and connection handler for Chem Search
 */

class Database {
    private static $instance = null;
    private $connection;
    
    // Database configuration
    private $host;
    private $database;
    private $username;
    private $password;
    private $charset = 'utf8mb4';
    
    private function __construct() {
        // Load database configuration from environment or defaults
        $this->host = $_ENV['DB_HOST'] ?? 'localhost';
        $this->database = $_ENV['DB_NAME'] ?? 'chem_search';
        $this->username = $_ENV['DB_USER'] ?? 'root';
        $this->password = $_ENV['DB_PASS'] ?? '';
        
        $this->connect();
    }
    
    public static function getInstance() {
        if (self::$instance === null) {
            self::$instance = new self();
        }
        return self::$instance;
    }
    
    private function connect() {
        try {
            $dsn = "mysql:host={$this->host};dbname={$this->database};charset={$this->charset}";
            $options = [
                PDO::ATTR_ERRMODE => PDO::ERRMODE_EXCEPTION,
                PDO::ATTR_DEFAULT_FETCH_MODE => PDO::FETCH_ASSOC,
                PDO::ATTR_EMULATE_PREPARES => false
            ];
            
            $this->connection = new PDO($dsn, $this->username, $this->password, $options);
            
            // Set charset after connection
            $this->connection->exec("SET NAMES {$this->charset}");
        } catch (PDOException $e) {
            ErrorHandler::logError("Database connection failed: " . $e->getMessage());
            throw new Exception("Database connection failed");
        }
    }
    
    public function getConnection() {
        return $this->connection;
    }
    
    public function query($sql, $params = []) {
        try {
            $stmt = $this->connection->prepare($sql);
            $stmt->execute($params);
            return $stmt;
        } catch (PDOException $e) {
            ErrorHandler::logError("Database query failed: " . $e->getMessage() . " SQL: " . $sql);
            throw new Exception("Database query failed");
        }
    }
    
    public function fetchOne($sql, $params = []) {
        $stmt = $this->query($sql, $params);
        return $stmt->fetch();
    }
    
    public function fetchAll($sql, $params = []) {
        $stmt = $this->query($sql, $params);
        return $stmt->fetchAll();
    }
    
    public function insert($table, $data) {
        $columns = implode(',', array_keys($data));
        $placeholders = ':' . implode(', :', array_keys($data));
        
        $sql = "INSERT INTO {$table} ({$columns}) VALUES ({$placeholders})";
        $this->query($sql, $data);
        
        return $this->connection->lastInsertId();
    }
    
    public function update($table, $data, $where, $whereParams = []) {
        $set = [];
        foreach (array_keys($data) as $key) {
            $set[] = "{$key} = :{$key}";
        }
        $setClause = implode(', ', $set);
        
        $sql = "UPDATE {$table} SET {$setClause} WHERE {$where}";
        $params = array_merge($data, $whereParams);
        
        return $this->query($sql, $params);
    }
    
    public function delete($table, $where, $params = []) {
        $sql = "DELETE FROM {$table} WHERE {$where}";
        return $this->query($sql, $params);
    }
    
    public function beginTransaction() {
        return $this->connection->beginTransaction();
    }
    
    public function commit() {
        return $this->connection->commit();
    }
    
    public function rollback() {
        return $this->connection->rollback();
    }
    
    public function tableExists($tableName) {
        $sql = "SELECT COUNT(*) as count FROM information_schema.tables WHERE table_schema = DATABASE() AND table_name = ?";
        $result = $this->fetchOne($sql, [$tableName]);
        return $result && $result['count'] > 0;
    }
    
    public function initializeDatabase() {
        // Check if tables exist, if not, create them
        $requiredTables = [
            'chemicals', 'molecular_properties', 'drug_properties', 
            'physical_properties', 'safety_properties', 'chemical_synonyms',
            'analysis_cache', 'search_history'
        ];
        
        $missingTables = [];
        foreach ($requiredTables as $table) {
            if (!$this->tableExists($table)) {
                $missingTables[] = $table;
            }
        }
        
        if (!empty($missingTables)) {
            throw new Exception("Missing database tables: " . implode(', ', $missingTables) . 
                              ". Please run the database schema script.");
        }
        
        return true;
    }
}

// Prevent cloning and unserialization
/**
 * Additional Database class methods for security
 */
class DatabaseSingleton extends Database {
    public function __clone() {
        throw new Exception("Cannot clone Database instance");
    }
    
    public function __wakeup() {
        throw new Exception("Cannot unserialize Database instance");
    }
}
?>
