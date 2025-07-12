<?php
// Test PHP configuration for Chem Search
echo "<h1>Chem Search - PHP Configuration Test</h1>\n";
echo "<hr>\n";

// Basic PHP info
echo "<h2>PHP Version</h2>\n";
echo "<p>PHP Version: " . phpversion() . "</p>\n";
echo "<p>Server Software: " . $_SERVER['SERVER_SOFTWARE'] . "</p>\n";
echo "<p>Document Root: " . $_SERVER['DOCUMENT_ROOT'] . "</p>\n";
echo "<p>Current Script: " . $_SERVER['SCRIPT_NAME'] . "</p>\n";

// Check required extensions
echo "<h2>Required Extensions</h2>\n";
$required_extensions = ['curl', 'json', 'mbstring', 'openssl'];
echo "<ul>\n";
foreach ($required_extensions as $ext) {
    $status = extension_loaded($ext) ? '✅ Loaded' : '❌ Missing';
    echo "<li><strong>$ext:</strong> $status</li>\n";
}
echo "</ul>\n";

// Check if .env file exists and is readable
echo "<h2>Environment Configuration</h2>\n";
$env_file = __DIR__ . '/.env';
if (file_exists($env_file)) {
    echo "<p>✅ .env file found: $env_file</p>\n";
    if (is_readable($env_file)) {
        echo "<p>✅ .env file is readable</p>\n";
    } else {
        echo "<p>❌ .env file is not readable</p>\n";
    }
} else {
    echo "<p>❌ .env file not found: $env_file</p>\n";
}

// Check Python path
echo "<h2>Python Configuration</h2>\n";
$python_paths = [
    'C:/xampp/htdocs/chem-search/project/.venv/Scripts/python.exe',
    'C:/Python311/python.exe',
    'C:/Python310/python.exe',
    'C:/Python39/python.exe'
];

foreach ($python_paths as $path) {
    if (file_exists($path)) {
        echo "<p>✅ Python found: $path</p>\n";
        
        // Test Python packages
        echo "<h3>Python Package Verification</h3>\n";
        $packages = ['rdkit', 'chemcrow', 'openai', 'langchain', 'requests', 'dotenv'];
        
        $cmd = '"' . $path . '" -c "';
        $cmd .= 'import sys; ';
        $cmd .= 'packages = ' . json_encode($packages) . '; ';
        $cmd .= 'for pkg in packages: ';
        $cmd .= '    try: ';
        $cmd .= '        __import__(pkg); ';
        $cmd .= '        print(f\\"✅ {pkg} - OK\\"); ';
        $cmd .= '    except ImportError: ';
        $cmd .= '        print(f\\"❌ {pkg} - MISSING\\")';
        $cmd .= '"';
        
        echo "<div style='background: #f8f9fa; padding: 10px; border: 1px solid #dee2e6; border-radius: 5px; font-family: monospace;'>\n";
        echo "<strong>Package Status:</strong><br>\n";
        
        $output = shell_exec($cmd . ' 2>&1');
        if ($output) {
            echo nl2br(htmlspecialchars($output));
        } else {
            echo "❌ Unable to check Python packages<br>\n";
        }
        
        echo "</div>\n";
        echo "<p><small><a href='project/verify_packages.py' target='_blank'>📋 Run detailed verification script</a></small></p>\n";
        
        break;
    }
}

// Check project directories
echo "<h2>Project Structure</h2>\n";
$important_dirs = ['tools', 'includes', 'assets', 'project', 'python'];
echo "<ul>\n";
foreach ($important_dirs as $dir) {
    $status = is_dir(__DIR__ . '/' . $dir) ? '✅ Exists' : '❌ Missing';
    echo "<li><strong>$dir/:</strong> $status</li>\n";
}
echo "</ul>\n";

// Check important files
echo "<h2>Important Files</h2>\n";
$important_files = [
    'index.php',
    'includes/config.php',
    'tools/molecule_analyzer.php',
    'tools/ai_assistant.php'
];
echo "<ul>\n";
foreach ($important_files as $file) {
    $status = file_exists(__DIR__ . '/' . $file) ? '✅ Exists' : '❌ Missing';
    echo "<li><strong>$file:</strong> $status</li>\n";
}
echo "</ul>\n";

echo "<hr>\n";
echo "<h2>Quick Links</h2>\n";
echo "<ul>\n";
echo "<li><a href='index.php'>🏠 Home Page</a></li>\n";
echo "<li><a href='tools/molecule_analyzer.php'>🧪 Molecule Analyzer</a></li>\n";
echo "<li><a href='tools/ai_assistant.php'>🤖 AI Assistant</a></li>\n";
echo "<li><a href='about.php'>ℹ️ About</a></li>\n";
echo "</ul>\n";

echo "<hr>\n";
echo "<p><em>Generated on: " . date('Y-m-d H:i:s') . "</em></p>\n";
?>
