<?php
// Configuration file for Chem Search

// Load environment variables if .env file exists
$envFile = __DIR__ . '/../.env';
if (file_exists($envFile)) {
    $env = parse_ini_file($envFile, false, INI_SCANNER_TYPED);
    if ($env !== false && is_array($env)) {
        foreach ($env as $key => $value) {
            if (!isset($_ENV[$key])) {
                $_ENV[$key] = $value;
            }
        }
    }
}

// Python environment paths
define('PYTHON_PATH', $_ENV['PYTHON_PATH'] ?? 'C:/xampp/htdocs/chemcrow/project/.venv/Scripts/python.exe');
define('PYTHON_SCRIPTS_DIR', realpath(__DIR__ . '/../python'));

// API configuration
define('OPENAI_API_KEY', $_ENV['OPENAI_API_KEY'] ?? getenv('OPENAI_API_KEY') ?: '');

// Application settings
define('APP_NAME', $_ENV['APP_NAME'] ?? 'Chem Search');
define('APP_VERSION', $_ENV['APP_VERSION'] ?? '1.0.0');
define('APP_DEBUG', filter_var($_ENV['APP_DEBUG'] ?? false, FILTER_VALIDATE_BOOLEAN));
define('APP_URL', $_ENV['APP_URL'] ?? 'http://localhost/chemcrow');

// Default models and settings
define('DEFAULT_AI_MODEL', $_ENV['DEFAULT_AI_MODEL'] ?? 'gpt-3.5-turbo');
define('DEFAULT_TEMPERATURE', (float)($_ENV['DEFAULT_TEMPERATURE'] ?? 0.1));

// Security settings
define('SESSION_LIFETIME', (int)($_ENV['SESSION_LIFETIME'] ?? 3600));
define('RATE_LIMIT_REQUESTS', (int)($_ENV['RATE_LIMIT_REQUESTS'] ?? 10));
define('RATE_LIMIT_WINDOW', (int)($_ENV['RATE_LIMIT_WINDOW'] ?? 60));

// PHP settings
ini_set('display_errors', APP_DEBUG ? 1 : 0);
error_reporting(APP_DEBUG ? E_ALL : E_ERROR);

// Include error handler AFTER constants are defined
require_once __DIR__ . '/error_handler.php';

// Include database if enabled
$useDatabase = filter_var($_ENV['USE_DATABASE'] ?? true, FILTER_VALIDATE_BOOLEAN);
if ($useDatabase) {
    require_once __DIR__ . '/database.php';
    
    // Initialize database connection
    try {
        $db = Database::getInstance();
        $db->initializeDatabase();
        define('DATABASE_ENABLED', true);
    } catch (Exception $e) {
        ErrorHandler::logError("Database initialization failed: " . $e->getMessage());
        define('DATABASE_ENABLED', false);
        
        if (APP_DEBUG) {
            echo "<div class='alert alert-warning'>Database not available: " . $e->getMessage() . "</div>";
        }
    }
} else {
    define('DATABASE_ENABLED', false);
}

// Security and validation functions
function sanitizeInput($input) {
    if (is_array($input)) {
        return array_map('sanitizeInput', $input);
    }
    return htmlspecialchars(trim($input), ENT_QUOTES, 'UTF-8');
}

function validateSMILES($smiles) {
    // Basic SMILES validation - check for common chemical notation patterns
    if (empty($smiles) || strlen($smiles) > 500) {
        return false;
    }
    
    // Allow only valid SMILES characters
    if (!preg_match('/^[A-Za-z0-9\[\]()=#@+\-\/\\\\\.]+$/', $smiles)) {
        return false;
    }
    
    return true;
}

function validateFilename($filename) {
    // Validate filename for upload/processing
    if (empty($filename) || strlen($filename) > 255) {
        return false;
    }
    
    // Check for path traversal attempts
    if (strpos($filename, '..') !== false || strpos($filename, '/') !== false || strpos($filename, '\\') !== false) {
        return false;
    }
    
    return true;
}

function generateCSRFToken() {
    if (session_status() == PHP_SESSION_NONE) {
        session_start();
    }
    
    if (!isset($_SESSION['csrf_token'])) {
        $_SESSION['csrf_token'] = bin2hex(random_bytes(32));
    }
    
    return $_SESSION['csrf_token'];
}

function validateCSRFToken($token) {
    if (session_status() == PHP_SESSION_NONE) {
        session_start();
    }
    
    return isset($_SESSION['csrf_token']) && hash_equals($_SESSION['csrf_token'], $token);
}

function rateLimitCheck($action = 'default', $limit = null, $window = null) {
    if (session_status() == PHP_SESSION_NONE) {
        session_start();
    }
    
    $limit = $limit ?? RATE_LIMIT_REQUESTS;
    $window = $window ?? RATE_LIMIT_WINDOW;
    
    $key = 'rate_limit_' . $action;
    $now = time();
    
    if (!isset($_SESSION[$key])) {
        $_SESSION[$key] = ['count' => 1, 'start' => $now];
        return true;
    }
    
    $data = $_SESSION[$key];
    
    if (($now - $data['start']) > $window) {
        $_SESSION[$key] = ['count' => 1, 'start' => $now];
        return true;
    }
    
    if ($data['count'] >= $limit) {
        return false;
    }
    
    $_SESSION[$key]['count']++;
    return true;
}

function formatError($message) {
    return [
        'error' => true,
        'message' => $message,
        'timestamp' => date('Y-m-d H:i:s')
    ];
}

function formatSuccess($data) {
    return [
        'success' => true,
        'data' => $data,
        'timestamp' => date('Y-m-d H:i:s')
    ];
}

// Check if Python environment is available
function checkPythonEnvironment() {
    $command = PYTHON_PATH . ' --version 2>&1';
    $output = shell_exec($command);
    return !empty($output) && strpos($output, 'Python') !== false;
}

// Check if required Python packages are installed
function checkRequiredPackages() {
    $required_packages = ['rdkit', 'chemcrow'];
    $missing = [];
    
    foreach ($required_packages as $package) {
        $command = PYTHON_PATH . ' -c "import ' . $package . '" 2>&1';
        $output = shell_exec($command);
        if (!empty($output) && strpos($output, 'ModuleNotFoundError') !== false) {
            $missing[] = $package;
        }
    }
    
    return $missing;
}
?>
