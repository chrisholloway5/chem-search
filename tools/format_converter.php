<?php
session_start();
require_once '../includes/config.php';

// Handle form submission
$result = null;
$error = null;

if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    $input = $_POST['input'] ?? '';
    $input_type = $_POST['input_type'] ?? 'smiles';
    $output_type = $_POST['output_type'] ?? 'smiles';
    
    if (!empty($input)) {
        try {
            $result = convertFormat($input, $input_type, $output_type);
        } catch (Exception $e) {
            $error = $e->getMessage();
        }
    } else {
        $error = "Please enter a chemical structure.";
    }
}

function convertFormat($input, $input_type, $output_type) {
    $pythonScript = realpath('../python/format_converter.py');
    $pythonPath = 'C:/xampp/htdocs/chemcrow/project/.venv/Scripts/python.exe';
    
    $escapedInput = escapeshellarg($input);
    $escapedInputType = escapeshellarg($input_type);
    $escapedOutputType = escapeshellarg($output_type);
    
    $command = "$pythonPath \"$pythonScript\" $escapedInput $escapedInputType $escapedOutputType";
    $output = shell_exec($command);
    
    if ($output === null) {
        throw new Exception("Failed to execute conversion script");
    }
    
    $data = json_decode($output, true);
    
    if (json_last_error() !== JSON_ERROR_NONE) {
        throw new Exception("Invalid JSON response from conversion script");
    }
    
    if (isset($data['error'])) {
        throw new Exception($data['error']);
    }
    
    return $data;
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Format Converter - Chem Search</title>
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
                    <i class="fas fa-exchange-alt text-success"></i> Chemical Format Converter
                </h1>
                <p class="lead">Convert between different chemical structure formats (SMILES, InChI, molecular formula, etc.).</p>
            </div>
        </div>

        <!-- Conversion Form -->
        <div class="row mb-4">
            <div class="col-lg-10 mx-auto">
                <div class="card shadow">
                    <div class="card-body">
                        <form method="POST" class="needs-validation" novalidate>
                            <div class="row mb-3">
                                <div class="col-md-6">
                                    <label for="input_type" class="form-label">Input Format</label>
                                    <select class="form-select" id="input_type" name="input_type" required>
                                        <option value="smiles" <?php echo ($_POST['input_type'] ?? '') === 'smiles' ? 'selected' : ''; ?>>SMILES</option>
                                        <option value="inchi" <?php echo ($_POST['input_type'] ?? '') === 'inchi' ? 'selected' : ''; ?>>InChI</option>
                                        <option value="mol" <?php echo ($_POST['input_type'] ?? '') === 'mol' ? 'selected' : ''; ?>>MOL Block</option>
                                    </select>
                                </div>
                                <div class="col-md-6">
                                    <label for="output_type" class="form-label">Output Format</label>
                                    <select class="form-select" id="output_type" name="output_type" required>
                                        <option value="smiles" <?php echo ($_POST['output_type'] ?? '') === 'smiles' ? 'selected' : ''; ?>>SMILES</option>
                                        <option value="inchi" <?php echo ($_POST['output_type'] ?? '') === 'inchi' ? 'selected' : ''; ?>>InChI</option>
                                        <option value="inchi_key" <?php echo ($_POST['output_type'] ?? '') === 'inchi_key' ? 'selected' : ''; ?>>InChI Key</option>
                                        <option value="formula" <?php echo ($_POST['output_type'] ?? '') === 'formula' ? 'selected' : ''; ?>>Molecular Formula</option>
                                        <option value="mol" <?php echo ($_POST['output_type'] ?? '') === 'mol' ? 'selected' : ''; ?>>MOL Block</option>
                                    </select>
                                </div>
                            </div>
                            
                            <div class="mb-3">
                                <label for="input" class="form-label">Chemical Structure</label>
                                <textarea class="form-control" id="input" name="input" rows="4" 
                                          placeholder="Enter your chemical structure here..."
                                          required><?php echo htmlspecialchars($_POST['input'] ?? ''); ?></textarea>
                                <div class="invalid-feedback">
                                    Please enter a chemical structure.
                                </div>
                            </div>
                            
                            <div class="d-grid">
                                <button type="submit" class="btn btn-success btn-lg">
                                    <i class="fas fa-exchange-alt"></i> Convert Format
                                </button>
                            </div>
                        </form>
                    </div>
                </div>
            </div>
        </div>

        <!-- Examples -->
        <div class="row mb-4">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-lightbulb text-warning"></i> Example Conversions</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-4">
                                <h6>SMILES Examples:</h6>
                                <ul class="list-unstyled">
                                    <li><code>CCO</code> - Ethanol</li>
                                    <li><code>CN1C=NC2=C1C(=O)N(C(=O)N2C)C</code> - Caffeine</li>
                                    <li><code>CC(=O)OC1=CC=CC=C1C(=O)O</code> - Aspirin</li>
                                </ul>
                            </div>
                            <div class="col-md-4">
                                <h6>InChI Examples:</h6>
                                <ul class="list-unstyled">
                                    <li><small><code>InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3</code></small> - Ethanol</li>
                                </ul>
                            </div>
                            <div class="col-md-4">
                                <h6>Supported Formats:</h6>
                                <ul class="list-unstyled">
                                    <li>✅ SMILES notation</li>
                                    <li>✅ InChI strings</li>
                                    <li>✅ InChI Keys</li>
                                    <li>✅ Molecular formulas</li>
                                    <li>✅ MOL file format</li>
                                </ul>
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
        <!-- Conversion Results -->
        <div class="row">
            <div class="col-lg-12">
                <div class="card border-success shadow fade-in">
                    <div class="card-header bg-success text-white">
                        <h5><i class="fas fa-check-circle"></i> Conversion Results</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h6>Input (<?php echo strtoupper(htmlspecialchars($_POST['input_type'])); ?>):</h6>
                                <div class="bg-light p-3 rounded mb-3">
                                    <code><?php echo htmlspecialchars($_POST['input']); ?></code>
                                </div>
                            </div>
                            <div class="col-md-6">
                                <h6>Output (<?php echo strtoupper(htmlspecialchars($_POST['output_type'])); ?>):</h6>
                                <div class="bg-success bg-opacity-10 p-3 rounded mb-3">
                                    <code class="text-success"><?php echo htmlspecialchars($result['output']); ?></code>
                                </div>
                                <button class="btn btn-outline-success btn-sm" onclick="copyToClipboard('<?php echo htmlspecialchars($result['output']); ?>')">
                                    <i class="fas fa-copy"></i> Copy Result
                                </button>
                            </div>
                        </div>
                        
                        <?php if (isset($result['additional_info'])): ?>
                        <hr>
                        <h6>Additional Information:</h6>
                        <div class="row">
                            <?php foreach ($result['additional_info'] as $key => $value): ?>
                            <div class="col-md-6 col-lg-4 mb-2">
                                <strong><?php echo ucfirst(str_replace('_', ' ', htmlspecialchars($key))); ?>:</strong>
                                <span class="text-muted"><?php echo htmlspecialchars($value); ?></span>
                            </div>
                            <?php endforeach; ?>
                        </div>
                        <?php endif; ?>
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

        // Copy to clipboard function
        function copyToClipboard(text) {
            navigator.clipboard.writeText(text).then(function() {
                // Show success message
                const btn = event.target.closest('button');
                const originalText = btn.innerHTML;
                btn.innerHTML = '<i class="fas fa-check"></i> Copied!';
                btn.classList.remove('btn-outline-success');
                btn.classList.add('btn-success');
                
                setTimeout(function() {
                    btn.innerHTML = originalText;
                    btn.classList.remove('btn-success');
                    btn.classList.add('btn-outline-success');
                }, 2000);
            });
        }
    </script>
</body>
</html>
