<?php
session_start();
require_once dirname(__DIR__) . '/includes/config.php';

// Initialize CSRF token
$csrfToken = generateCSRFToken();

// Handle form submission
$result = null;
$error = null;

if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    // Validate CSRF token
    if (!validateCSRFToken($_POST['csrf_token'] ?? '')) {
        $error = "Invalid security token. Please refresh the page and try again.";
    } else {
        // Rate limiting
        if (!rateLimitCheck('molecule_analysis', 5, 60)) {
            $error = "Too many requests. Please wait a moment before trying again.";
        } else {
            $smiles = sanitizeInput($_POST['smiles'] ?? '');
            
            if (!validateSMILES($smiles)) {
                $error = "Invalid SMILES notation. Please check your input.";
            } else {
                try {
                    // Log the analysis request
                    ErrorHandler::logActivity('molecule_analysis', "SMILES: $smiles");
                    
                    // Call Python script to analyze molecule
                    $result = analyzeMolecule($smiles);
                } catch (Exception $e) {
                    $error = "Analysis failed: " . $e->getMessage();
                    ErrorHandler::logError("Molecule analysis error: " . $e->getMessage());
                }
            }
        }
    }
} elseif (!empty($_GET['smiles'])) {
    // Handle GET request (for backwards compatibility)
    $smiles = sanitizeInput($_GET['smiles']);
    if (validateSMILES($smiles)) {
        try {
            $result = analyzeMolecule($smiles);
        } catch (Exception $e) {
            $error = "Analysis failed: " . $e->getMessage();
        }
    } else {
        $error = "Invalid SMILES notation.";
    }
}

function analyzeMolecule($smiles) {
    // Validate inputs
    if (!validateSMILES($smiles)) {
        throw new Exception("Invalid SMILES notation");
    }
    
    // Initialize chemical model
    require_once dirname(__DIR__) . '/includes/chemical_model.php';
    $chemical = new Chemical();
    
    try {
        // Start database transaction
        $db = Database::getInstance();
        $db->beginTransaction();
        
        // Check if chemical already exists in database
        $existingChemical = $chemical->findBySmiles($smiles);
        
        if ($existingChemical) {
            // Check for cached analysis (24 hours)
            $cachedAnalysis = $chemical->getCachedAnalysis($existingChemical['id'], 'comprehensive_analysis', 86400);
            
            if ($cachedAnalysis) {
                $db->commit();
                return array_merge($existingChemical, [
                    'properties' => $cachedAnalysis['raw_data'],
                    'from_cache' => true,
                    'cache_date' => $cachedAnalysis['analysis_date']
                ]);
            }
        }
        
        // Execute Python analysis
        $analysisResult = executePythonAnalysis($smiles);
        
        if (!$analysisResult['success']) {
            $db->rollback();
            throw new Exception($analysisResult['error']);
        }
        
        // Store or update chemical in database
        if ($existingChemical) {
            $chemicalId = $existingChemical['id'];
            
            // Update basic properties if we have new data
            if (isset($analysisResult['basic_properties'])) {
                $chemical->updateChemical($chemicalId, $analysisResult['basic_properties']);
            }
        } else {
            // Create new chemical entry
            $basicData = array_merge(['smiles' => $smiles], $analysisResult['basic_properties'] ?? []);
            $newChemical = $chemical->createChemical($basicData);
            $chemicalId = $newChemical['id'];
        }
        
        // Store molecular properties
        if (isset($analysisResult['molecular_properties'])) {
            $chemical->storeMolecularProperties($chemicalId, $analysisResult['molecular_properties']);
        }
        
        // Store drug properties
        if (isset($analysisResult['drug_properties'])) {
            $chemical->storeDrugProperties($chemicalId, $analysisResult['drug_properties']);
        }
        
        // Store physical properties
        if (isset($analysisResult['physical_properties'])) {
            $chemical->storePhysicalProperties($chemicalId, $analysisResult['physical_properties']);
        }
        
        // Cache the full analysis
        $chemical->cacheAnalysis(
            $chemicalId, 
            'comprehensive_analysis', 
            $analysisResult,
            'rdkit',
            $analysisResult['analysis_metadata']['rdkit_version'] ?? '1.0'
        );
        
        // Log the search
        $sessionId = session_id();
        $ipAddress = $_SERVER['REMOTE_ADDR'] ?? 'unknown';
        $userAgent = $_SERVER['HTTP_USER_AGENT'] ?? 'unknown';
        $chemical->logSearch($sessionId, $chemicalId, $smiles, 'smiles', $ipAddress, $userAgent);
        
        $db->commit();
        
        // Get the full profile
        $fullProfile = $chemical->getFullProfile($chemicalId);
        
        return array_merge($fullProfile, [
            'detailed_analysis' => $analysisResult,
            'from_cache' => false,
            'analysis_date' => date('Y-m-d H:i:s')
        ]);
        
    } catch (Exception $e) {
        $db->rollback();
        throw $e;
    }
}

function executePythonAnalysis($smiles) {
    // Validate Python path
    if (!file_exists(PYTHON_PATH)) {
        throw new Exception("Python interpreter not found");
    }
    
    // Path to Python script
    $pythonScript = realpath(dirname(__DIR__) . '/python/analyze_molecule.py');
    if (!$pythonScript || !file_exists($pythonScript)) {
        throw new Exception("Python analysis script not found");
    }
    
    // Escape arguments properly
    $escapedSmiles = escapeshellarg($smiles);
    $escapedPythonPath = escapeshellarg(PYTHON_PATH);
    $escapedScriptPath = escapeshellarg($pythonScript);
    
    // Execute Python script with timeout
    $command = "$escapedPythonPath $escapedScriptPath $escapedSmiles 2>&1";
    
    // Set timeout to prevent hanging
    $descriptorspec = [
        0 => ["pipe", "r"],
        1 => ["pipe", "w"],
        2 => ["pipe", "w"]
    ];
    
    $process = proc_open($command, $descriptorspec, $pipes);
    
    if (!is_resource($process)) {
        throw new Exception("Failed to execute Python script");
    }
    
    // Close stdin and read output with timeout
    fclose($pipes[0]);
    
    stream_set_timeout($pipes[1], 30); // 30 second timeout
    $output = stream_get_contents($pipes[1]);
    $error_output = stream_get_contents($pipes[2]);
    
    fclose($pipes[1]);
    fclose($pipes[2]);
    
    $return_value = proc_close($process);
    
    if ($return_value !== 0) {
        throw new Exception("Python script execution failed: " . $error_output);
    }
    
    if (empty($output)) {
        throw new Exception("No output from Python script");
    }
    
    $data = json_decode($output, true);
    
    if (json_last_error() !== JSON_ERROR_NONE) {
        throw new Exception("Invalid JSON response: " . json_last_error_msg());
    }
    
    return $data;
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Molecule Analyzer - Chem Search</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <link href="../assets/css/style.css" rel="stylesheet">
</head>
<body>
    <?php include '../includes/header.php'; ?>
    
    <main class="container my-5">
        <div class="row">
            <div class="col-lg-12">
                <h1 class="mb-4">
                    <i class="fas fa-atom text-primary"></i> Molecule Analyzer
                </h1>
                <p class="lead">Enter a SMILES notation to analyze molecular properties using RDKit.</p>
            </div>
        </div>

        <!-- Input Form -->
        <div class="row mb-4">
            <div class="col-lg-8 mx-auto">
                <div class="card shadow">
                    <div class="card-body">
                        <form method="POST" class="needs-validation" novalidate>
                            <input type="hidden" name="csrf_token" value="<?php echo $csrfToken; ?>">
                            <div class="mb-3">
                                <label for="smiles" class="form-label">SMILES Notation</label>
                                <input type="text" class="form-control" id="smiles" name="smiles" 
                                       value="<?php echo htmlspecialchars($_POST['smiles'] ?? $_GET['smiles'] ?? ''); ?>" 
                                       placeholder="e.g., CN1C=NC2=C1C(=O)N(C(=O)N2C)C (caffeine)"
                                       required>
                                <div class="invalid-feedback">
                                    Please enter a valid SMILES notation.
                                </div>
                            </div>
                            <div class="d-grid">
                                <button type="submit" class="btn btn-primary btn-lg">
                                    <i class="fas fa-search"></i> Analyze Molecule
                                </button>
                            </div>
                        </form>
                    </div>
                </div>
            </div>
        </div>

        <!-- Example Molecules -->
        <div class="row mb-4">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-lightbulb text-warning"></i> Example Molecules</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-4 mb-2">
                                <a href="?smiles=O" class="btn btn-outline-primary btn-sm w-100">
                                    Water (H₂O)
                                </a>
                            </div>
                            <div class="col-md-4 mb-2">
                                <a href="?smiles=CCO" class="btn btn-outline-primary btn-sm w-100">
                                    Ethanol
                                </a>
                            </div>
                            <div class="col-md-4 mb-2">
                                <a href="?smiles=CN1C=NC2=C1C(=O)N(C(=O)N2C)C" class="btn btn-outline-primary btn-sm w-100">
                                    Caffeine
                                </a>
                            </div>
                            <div class="col-md-4 mb-2">
                                <a href="?smiles=CC(=O)OC1=CC=CC=C1C(=O)O" class="btn btn-outline-primary btn-sm w-100">
                                    Aspirin
                                </a>
                            </div>
                            <div class="col-md-4 mb-2">
                                <a href="?smiles=CC(=O)NC1=CC=C(C=C1)O" class="btn btn-outline-primary btn-sm w-100">
                                    Tylenol
                                </a>
                            </div>
                            <div class="col-md-4 mb-2">
                                <a href="?smiles=CC(C)CC1=CC=C(C=C1)C(C)C(=O)O" class="btn btn-outline-primary btn-sm w-100">
                                    Ibuprofen
                                </a>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <?php if ($error): ?>
        <div class="alert alert-danger" role="alert">
            <i class="fas fa-exclamation-triangle"></i> <strong>Error:</strong> <?php echo htmlspecialchars($error); ?>
        </div>
        <?php endif; ?>

        <?php if ($result): ?>
        <!-- Results -->
        <div class="row">
            <div class="col-lg-12">
                <div class="molecule-card fade-in">
                    <h3 class="text-center mb-4">
                        <i class="fas fa-flask text-success"></i> Analysis Results
                    </h3>
                    
                    <!-- Basic Information -->
                    <div class="card mb-4">
                        <div class="card-header bg-primary text-white">
                            <h5><i class="fas fa-info-circle"></i> Basic Information</h5>
                        </div>
                        <div class="card-body">
                            <div class="row">
                                <div class="col-md-6">
                                    <strong>SMILES:</strong> <code><?php echo htmlspecialchars($result['smiles']); ?></code>
                                </div>
                                <div class="col-md-6">
                                    <strong>Molecular Formula:</strong> <code><?php echo htmlspecialchars($result['molecular_formula']); ?></code>
                                </div>
                            </div>
                            <?php if (isset($result['inchi'])): ?>
                            <div class="row mt-2">
                                <div class="col-12">
                                    <strong>InChI:</strong> <small><code><?php echo htmlspecialchars($result['inchi']); ?></code></small>
                                </div>
                            </div>
                            <?php endif; ?>
                            <?php if (isset($result['inchi_key'])): ?>
                            <div class="row mt-2">
                                <div class="col-12">
                                    <strong>InChI Key:</strong> <code><?php echo htmlspecialchars($result['inchi_key']); ?></code>
                                </div>
                            </div>
                            <?php endif; ?>
                        </div>
                    </div>

                    <!-- Molecular Properties -->
                    <div class="card">
                        <div class="card-header bg-success text-white">
                            <h5><i class="fas fa-calculator"></i> Molecular Properties</h5>
                        </div>
                        <div class="card-body">
                            <div class="property-grid">
                                <div class="property-item">
                                    <h6>Molecular Weight</h6>
                                    <div class="value"><?php echo number_format($result['molecular_weight'], 2); ?></div>
                                    <div class="unit">g/mol</div>
                                </div>
                                
                                <div class="property-item">
                                    <h6>LogP</h6>
                                    <div class="value"><?php echo number_format($result['logp'], 2); ?></div>
                                    <div class="unit">lipophilicity</div>
                                </div>
                                
                                <div class="property-item">
                                    <h6>H-Bond Donors</h6>
                                    <div class="value"><?php echo $result['hbd']; ?></div>
                                    <div class="unit">count</div>
                                </div>
                                
                                <div class="property-item">
                                    <h6>H-Bond Acceptors</h6>
                                    <div class="value"><?php echo $result['hba']; ?></div>
                                    <div class="unit">count</div>
                                </div>
                                
                                <?php if (isset($result['tpsa'])): ?>
                                <div class="property-item">
                                    <h6>TPSA</h6>
                                    <div class="value"><?php echo number_format($result['tpsa'], 2); ?></div>
                                    <div class="unit">Ų</div>
                                </div>
                                <?php endif; ?>
                                
                                <?php if (isset($result['rotatable_bonds'])): ?>
                                <div class="property-item">
                                    <h6>Rotatable Bonds</h6>
                                    <div class="value"><?php echo $result['rotatable_bonds']; ?></div>
                                    <div class="unit">count</div>
                                </div>
                                <?php endif; ?>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        <?php endif; ?>
    </main>

    <?php include '../includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script>
        // Form validation
        (function() {
            'use strict';
            window.addEventListener('load', function() {
                var forms = document.getElementsByClassName('needs-validation');
                var validation = Array.prototype.filter.call(forms, function(form) {
                    form.addEventListener('submit', function(event) {
                        if (form.checkValidity() === false) {
                            event.preventDefault();
                            event.stopPropagation();
                        }
                        form.classList.add('was-validated');
                    }, false);
                });
            }, false);
        })();
    </script>
</body>
</html>
