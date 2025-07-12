<?php
/**
 * Link Checker Script - Tests all links in the website
 */

echo "<h1>🔗 Chem Search - Link Verification</h1>";
echo "<hr>";

// Base URLs to test
$base_url = "http://localhost/chem-search";
$pages_to_test = [
    '' => 'Home Page',
    '/about.php' => 'About Page',
    '/api_status.php' => 'API Status Page',
    '/test_setup.php' => 'Setup Test Page',
    '/tools/molecule_analyzer.php' => 'Molecule Analyzer',
    '/tools/format_converter.php' => 'Format Converter',
    '/tools/ai_assistant.php' => 'AI Assistant',
    '/tools/database_search.php' => 'Database Search',
    '/tools/structure_viewer.php' => 'Structure Viewer',
    '/tools/paper_search.php' => 'Paper Search',
    '/tools/database_manager.php' => 'Database Manager'
];

echo "<h2>📋 Testing Page Accessibility</h2>";
echo "<table border='1' style='border-collapse: collapse; width: 100%;'>";
echo "<tr><th>Page</th><th>URL</th><th>Status</th><th>Response</th></tr>";

$working_count = 0;
$broken_count = 0;

foreach ($pages_to_test as $path => $name) {
    $full_url = $base_url . $path;
    
    // Test if page exists and loads
    $ch = curl_init();
    curl_setopt($ch, CURLOPT_URL, $full_url);
    curl_setopt($ch, CURLOPT_RETURNTRANSFER, true);
    curl_setopt($ch, CURLOPT_FOLLOWLOCATION, true);
    curl_setopt($ch, CURLOPT_TIMEOUT, 10);
    curl_setopt($ch, CURLOPT_NOBODY, true); // Head request only
    
    $response = curl_exec($ch);
    $http_code = curl_getinfo($ch, CURLINFO_HTTP_CODE);
    $error = curl_error($ch);
    curl_close($ch);
    
    if ($http_code == 200) {
        echo "<tr style='background-color: #d4edda;'>";
        echo "<td>$name</td>";
        echo "<td><a href='$full_url' target='_blank'>$full_url</a></td>";
        echo "<td>✅ Working</td>";
        echo "<td>HTTP $http_code</td>";
        echo "</tr>";
        $working_count++;
    } else {
        echo "<tr style='background-color: #f8d7da;'>";
        echo "<td>$name</td>";
        echo "<td>$full_url</td>";
        echo "<td>❌ Broken</td>";
        echo "<td>HTTP $http_code" . ($error ? " - $error" : "") . "</td>";
        echo "</tr>";
        $broken_count++;
    }
}

echo "</table>";

echo "<h2>📊 Summary</h2>";
echo "<p><strong>✅ Working Links:</strong> $working_count</p>";
echo "<p><strong>❌ Broken Links:</strong> $broken_count</p>";

// Test file existence for assets
echo "<h2>📁 Testing Asset Files</h2>";
echo "<table border='1' style='border-collapse: collapse; width: 100%;'>";
echo "<tr><th>Asset Type</th><th>File Path</th><th>Status</th></tr>";

$assets_to_test = [
    'CSS' => [
        'assets/css/style.css'
    ],
    'JavaScript' => [
        'assets/js/main.js'
    ],
    'Images' => [
        // Add any images that should exist
    ]
];

foreach ($assets_to_test as $type => $files) {
    foreach ($files as $file) {
        $full_path = __DIR__ . '/' . $file;
        $url_path = $base_url . '/' . $file;
        
        if (file_exists($full_path)) {
            echo "<tr style='background-color: #d4edda;'>";
            echo "<td>$type</td>";
            echo "<td><a href='$url_path' target='_blank'>$file</a></td>";
            echo "<td>✅ Exists</td>";
            echo "</tr>";
        } else {
            echo "<tr style='background-color: #f8d7da;'>";
            echo "<td>$type</td>";
            echo "<td>$file</td>";
            echo "<td>❌ Missing</td>";
            echo "</tr>";
        }
    }
}

echo "</table>";

// Test navigation consistency
echo "<h2>🧭 Navigation Link Analysis</h2>";

function checkNavigationLinks($file_path, $page_name) {
    if (!file_exists($file_path)) {
        echo "<p>❌ <strong>$page_name</strong>: File not found</p>";
        return;
    }
    
    $content = file_get_contents($file_path);
    
    // Extract href links
    preg_match_all('/href=["\']([^"\']*)["\']/', $content, $matches);
    $links = $matches[1];
    
    // Filter relative links (exclude external URLs and anchors)
    $relative_links = array_filter($links, function($link) {
        return !preg_match('/^(https?:\/\/|#|mailto:)/', $link) && !empty($link);
    });
    
    echo "<details>";
    echo "<summary><strong>$page_name</strong> (" . count($relative_links) . " internal links)</summary>";
    echo "<ul>";
    
    foreach ($relative_links as $link) {
        // Determine the expected file path
        $base_dir = dirname($file_path);
        
        if (strpos($link, '/') === 0) {
            // Absolute path from web root
            $expected_path = __DIR__ . $link;
        } else {
            // Relative path
            $expected_path = $base_dir . '/' . $link;
        }
        
        // Clean up path
        $expected_path = realpath($expected_path);
        
        if ($expected_path && file_exists($expected_path)) {
            echo "<li>✅ <code>$link</code></li>";
        } else {
            echo "<li>❌ <code>$link</code> → Expected: " . ($expected_path ?: 'Path resolution failed') . "</li>";
        }
    }
    
    echo "</ul>";
    echo "</details>";
}

// Check key pages
$pages_to_analyze = [
    __DIR__ . '/index.php' => 'Home Page',
    __DIR__ . '/includes/header.php' => 'Header Navigation',
    __DIR__ . '/includes/footer.php' => 'Footer Navigation',
    __DIR__ . '/tools/molecule_analyzer.php' => 'Molecule Analyzer',
    __DIR__ . '/tools/ai_assistant.php' => 'AI Assistant'
];

foreach ($pages_to_analyze as $file => $name) {
    checkNavigationLinks($file, $name);
}

echo "<hr>";
echo "<p><em>Link check completed on " . date('Y-m-d H:i:s') . "</em></p>";
?>
