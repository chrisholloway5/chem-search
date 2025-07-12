<?php
session_start();
require_once 'includes/config.php';

// Check system status
$python_available = checkPythonEnvironment();
$missing_packages = checkRequiredPackages();
$api_key_set = !empty(OPENAI_API_KEY);
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>API Status - ChemCrow Web</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <link href="assets/css/style.css" rel="stylesheet">
</head>
<body>
    <?php include 'includes/header.php'; ?>
    
    <main class="container my-5">
        <div class="row">
            <div class="col-lg-12">
                <h1 class="mb-4">
                    <i class="fas fa-server text-info"></i> System Status
                </h1>
                <p class="lead">Check the status of ChemCrow Web components and dependencies.</p>
            </div>
        </div>

        <!-- System Status Cards -->
        <div class="row g-4">
            <!-- Python Environment -->
            <div class="col-md-6 col-lg-4">
                <div class="card h-100">
                    <div class="card-body text-center">
                        <?php if ($python_available): ?>
                        <i class="fas fa-check-circle fa-3x text-success mb-3"></i>
                        <h5 class="card-title text-success">Python Environment</h5>
                        <p class="card-text">Python is properly configured and accessible.</p>
                        <span class="badge bg-success">Available</span>
                        <?php else: ?>
                        <i class="fas fa-times-circle fa-3x text-danger mb-3"></i>
                        <h5 class="card-title text-danger">Python Environment</h5>
                        <p class="card-text">Python environment is not accessible.</p>
                        <span class="badge bg-danger">Unavailable</span>
                        <?php endif; ?>
                    </div>
                </div>
            </div>

            <!-- Required Packages -->
            <div class="col-md-6 col-lg-4">
                <div class="card h-100">
                    <div class="card-body text-center">
                        <?php if (empty($missing_packages)): ?>
                        <i class="fas fa-check-circle fa-3x text-success mb-3"></i>
                        <h5 class="card-title text-success">Python Packages</h5>
                        <p class="card-text">All required packages are installed.</p>
                        <span class="badge bg-success">Complete</span>
                        <?php else: ?>
                        <i class="fas fa-exclamation-triangle fa-3x text-warning mb-3"></i>
                        <h5 class="card-title text-warning">Python Packages</h5>
                        <p class="card-text">Some packages are missing.</p>
                        <span class="badge bg-warning">Incomplete</span>
                        <?php endif; ?>
                    </div>
                </div>
            </div>

            <!-- OpenAI API -->
            <div class="col-md-6 col-lg-4">
                <div class="card h-100">
                    <div class="card-body text-center">
                        <?php if ($api_key_set): ?>
                        <i class="fas fa-key fa-3x text-info mb-3"></i>
                        <h5 class="card-title text-info">OpenAI API Key</h5>
                        <p class="card-text">API key is configured.</p>
                        <span class="badge bg-info">Configured</span>
                        <?php else: ?>
                        <i class="fas fa-times-circle fa-3x text-danger mb-3"></i>
                        <h5 class="card-title text-danger">OpenAI API Key</h5>
                        <p class="card-text">API key is not set.</p>
                        <span class="badge bg-danger">Missing</span>
                        <?php endif; ?>
                    </div>
                </div>
            </div>
        </div>

        <!-- Detailed Status -->
        <div class="row mt-5">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-list-ul"></i> Detailed Status Information</h5>
                    </div>
                    <div class="card-body">
                        <!-- Python Status -->
                        <div class="mb-4">
                            <h6><i class="fab fa-python text-primary"></i> Python Environment</h6>
                            <ul class="list-group">
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    Python Executable
                                    <?php if ($python_available): ?>
                                    <span class="badge bg-success rounded-pill">Available</span>
                                    <?php else: ?>
                                    <span class="badge bg-danger rounded-pill">Not Found</span>
                                    <?php endif; ?>
                                </li>
                                <li class="list-group-item">
                                    <small class="text-muted">Path: <?php echo htmlspecialchars(PYTHON_PATH); ?></small>
                                </li>
                            </ul>
                        </div>

                        <!-- Package Status -->
                        <div class="mb-4">
                            <h6><i class="fas fa-box text-success"></i> Required Packages</h6>
                            <ul class="list-group">
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    RDKit (Chemistry Toolkit)
                                    <?php if (!in_array('rdkit', $missing_packages)): ?>
                                    <span class="badge bg-success rounded-pill">Installed</span>
                                    <?php else: ?>
                                    <span class="badge bg-danger rounded-pill">Missing</span>
                                    <?php endif; ?>
                                </li>
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    ChemCrow (AI Agent)
                                    <?php if (!in_array('chemcrow', $missing_packages)): ?>
                                    <span class="badge bg-success rounded-pill">Installed</span>
                                    <?php else: ?>
                                    <span class="badge bg-danger rounded-pill">Missing</span>
                                    <?php endif; ?>
                                </li>
                            </ul>
                        </div>

                        <!-- API Status -->
                        <div class="mb-4">
                            <h6><i class="fas fa-cloud text-info"></i> API Configuration</h6>
                            <ul class="list-group">
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    OpenAI API Key
                                    <?php if ($api_key_set): ?>
                                    <span class="badge bg-success rounded-pill">Set</span>
                                    <?php else: ?>
                                    <span class="badge bg-danger rounded-pill">Not Set</span>
                                    <?php endif; ?>
                                </li>
                                <?php if ($api_key_set): ?>
                                <li class="list-group-item">
                                    <small class="text-muted">Key: <?php echo substr(OPENAI_API_KEY, 0, 10) . '...' . substr(OPENAI_API_KEY, -4); ?></small>
                                </li>
                                <?php endif; ?>
                            </ul>
                        </div>

                        <!-- Features Availability -->
                        <div class="mb-4">
                            <h6><i class="fas fa-tools text-warning"></i> Feature Availability</h6>
                            <ul class="list-group">
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    Molecule Analyzer
                                    <?php if ($python_available && !in_array('rdkit', $missing_packages)): ?>
                                    <span class="badge bg-success rounded-pill">Available</span>
                                    <?php else: ?>
                                    <span class="badge bg-danger rounded-pill">Unavailable</span>
                                    <?php endif; ?>
                                </li>
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    Format Converter
                                    <?php if ($python_available && !in_array('rdkit', $missing_packages)): ?>
                                    <span class="badge bg-success rounded-pill">Available</span>
                                    <?php else: ?>
                                    <span class="badge bg-danger rounded-pill">Unavailable</span>
                                    <?php endif; ?>
                                </li>
                                <li class="list-group-item d-flex justify-content-between align-items-center">
                                    AI Assistant
                                    <?php if ($python_available && !in_array('chemcrow', $missing_packages) && $api_key_set): ?>
                                    <span class="badge bg-success rounded-pill">Available</span>
                                    <?php else: ?>
                                    <span class="badge bg-warning rounded-pill">Limited</span>
                                    <?php endif; ?>
                                </li>
                            </ul>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Troubleshooting -->
        <?php if (!$python_available || !empty($missing_packages) || !$api_key_set): ?>
        <div class="row mt-5">
            <div class="col-lg-12">
                <div class="card border-warning">
                    <div class="card-header bg-warning text-dark">
                        <h5><i class="fas fa-wrench"></i> Troubleshooting</h5>
                    </div>
                    <div class="card-body">
                        <?php if (!$python_available): ?>
                        <div class="alert alert-danger">
                            <h6><i class="fas fa-exclamation-triangle"></i> Python Environment Issue</h6>
                            <p>Python is not accessible at the configured path. Please:</p>
                            <ul>
                                <li>Verify Python is installed</li>
                                <li>Check the path in <code>includes/config.php</code></li>
                                <li>Ensure the virtual environment is activated</li>
                            </ul>
                        </div>
                        <?php endif; ?>

                        <?php if (!empty($missing_packages)): ?>
                        <div class="alert alert-warning">
                            <h6><i class="fas fa-box-open"></i> Missing Packages</h6>
                            <p>The following packages are missing:</p>
                            <ul>
                                <?php foreach ($missing_packages as $package): ?>
                                <li><code><?php echo htmlspecialchars($package); ?></code></li>
                                <?php endforeach; ?>
                            </ul>
                            <p>Install them using: <code>pip install <?php echo implode(' ', $missing_packages); ?></code></p>
                        </div>
                        <?php endif; ?>

                        <?php if (!$api_key_set): ?>
                        <div class="alert alert-info">
                            <h6><i class="fas fa-key"></i> OpenAI API Key</h6>
                            <p>Set your OpenAI API key as an environment variable:</p>
                            <code>$env:OPENAI_API_KEY="your-api-key-here"</code>
                        </div>
                        <?php endif; ?>
                    </div>
                </div>
            </div>
        </div>
        <?php endif; ?>
    </main>

    <?php include 'includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
