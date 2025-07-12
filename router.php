<?php
// Router script for PHP development server
$request_uri = parse_url($_SERVER['REQUEST_URI'], PHP_URL_PATH);

// Remove query parameters
$path = $request_uri;

// If requesting root, serve index.php
if ($path === '/' || $path === '') {
    include 'index.php';
    return true;
}

// If file exists, serve it
$file = __DIR__ . $path;
if (is_file($file)) {
    return false; // Let PHP serve the file
}

// For other routes, try to find the file
if (file_exists(__DIR__ . $path . '.php')) {
    include __DIR__ . $path . '.php';
    return true;
}

// If nothing found, return 404
http_response_code(404);
echo "404 - File not found";
return true;
?>
