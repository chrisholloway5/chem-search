<?php
/**
 * Centralized error handling and logging for Chem Search
 */

class ErrorHandler {
    private static $logFile = __DIR__ . '/../logs/error.log';
    
    public static function init() {
        // Ensure logs directory exists
        $logDir = dirname(self::$logFile);
        if (!is_dir($logDir)) {
            mkdir($logDir, 0755, true);
        }
        
        // Set custom error handler
        set_error_handler([self::class, 'handleError']);
        set_exception_handler([self::class, 'handleException']);
        register_shutdown_function([self::class, 'handleFatalError']);
    }
    
    public static function handleError($severity, $message, $file, $line) {
        if (!(error_reporting() & $severity)) {
            return false;
        }
        
        $errorMsg = sprintf(
            "[%s] Error: %s in %s on line %d",
            date('Y-m-d H:i:s'),
            $message,
            $file,
            $line
        );
        
        self::logError($errorMsg);
        
        if (defined('APP_DEBUG') && APP_DEBUG) {
            echo "<div class='alert alert-danger'>$errorMsg</div>";
        }
        
        return true;
    }
    
    public static function handleException($exception) {
        $errorMsg = sprintf(
            "[%s] Uncaught Exception: %s in %s on line %d",
            date('Y-m-d H:i:s'),
            $exception->getMessage(),
            $exception->getFile(),
            $exception->getLine()
        );
        
        self::logError($errorMsg);
        
        if (defined('APP_DEBUG') && APP_DEBUG) {
            echo "<div class='alert alert-danger'>$errorMsg</div>";
        } else {
            echo "<div class='alert alert-danger'>An unexpected error occurred. Please try again later.</div>";
        }
    }
    
    public static function handleFatalError() {
        $error = error_get_last();
        if ($error && $error['type'] === E_ERROR) {
            $errorMsg = sprintf(
                "[%s] Fatal Error: %s in %s on line %d",
                date('Y-m-d H:i:s'),
                $error['message'],
                $error['file'],
                $error['line']
            );
            
            self::logError($errorMsg);
        }
    }
    
    public static function logError($message) {
        $timestamp = date('Y-m-d H:i:s');
        $logEntry = "[$timestamp] $message" . PHP_EOL;
        file_put_contents(self::$logFile, $logEntry, FILE_APPEND | LOCK_EX);
    }
    
    public static function logActivity($action, $details = '') {
        $logFile = dirname(self::$logFile) . '/activity.log';
        $timestamp = date('Y-m-d H:i:s');
        $ip = $_SERVER['REMOTE_ADDR'] ?? 'unknown';
        $userAgent = $_SERVER['HTTP_USER_AGENT'] ?? 'unknown';
        
        $logEntry = sprintf(
            "[%s] IP: %s | Action: %s | Details: %s | User-Agent: %s",
            $timestamp,
            $ip,
            $action,
            $details,
            $userAgent
        ) . PHP_EOL;
        
        file_put_contents($logFile, $logEntry, FILE_APPEND | LOCK_EX);
    }
}

// Initialize error handling
ErrorHandler::init();
?>
