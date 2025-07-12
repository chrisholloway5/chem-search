<?php
session_start();
require_once '../includes/config.php';
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Structure Viewer - Chem Search</title>
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
                    <i class="fas fa-project-diagram text-danger"></i> Molecular Structure Viewer
                </h1>
                <p class="lead">Visualize molecular structures in 2D and 3D formats.</p>
            </div>
        </div>

        <!-- Input Form -->
        <div class="row mb-4">
            <div class="col-lg-8 mx-auto">
                <div class="card shadow">
                    <div class="card-body">
                        <form id="structureForm" class="needs-validation" novalidate>
                            <div class="mb-3">
                                <label for="smiles" class="form-label">SMILES Notation</label>
                                <input type="text" class="form-control" id="smiles" name="smiles" 
                                       value="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
                                       placeholder="e.g., CN1C=NC2=C1C(=O)N(C(=O)N2C)C (caffeine)"
                                       required>
                                <div class="invalid-feedback">
                                    Please enter a valid SMILES notation.
                                </div>
                            </div>
                            <div class="d-grid">
                                <button type="submit" class="btn btn-danger btn-lg">
                                    <i class="fas fa-eye"></i> Visualize Structure
                                </button>
                            </div>
                        </form>
                    </div>
                </div>
            </div>
        </div>

        <!-- Structure Display -->
        <div class="row">
            <div class="col-lg-12">
                <div class="card shadow">
                    <div class="card-header bg-danger text-white">
                        <h5><i class="fas fa-project-diagram"></i> Molecular Structure</h5>
                    </div>
                    <div class="card-body">
                        <!-- 2D Structure -->
                        <div class="row">
                            <div class="col-md-6">
                                <h6><i class="fas fa-draw-polygon"></i> 2D Structure</h6>
                                <div id="structure2d" class="structure-viewer">
                                    <div class="text-center py-5">
                                        <i class="fas fa-molecule fa-4x text-muted mb-3"></i>
                                        <p class="text-muted">Enter a SMILES notation to view the 2D structure</p>
                                    </div>
                                </div>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-cube"></i> 3D Structure (Coming Soon)</h6>
                                <div id="structure3d" class="structure-viewer">
                                    <div class="text-center py-5">
                                        <i class="fas fa-cubes fa-4x text-muted mb-3"></i>
                                        <p class="text-muted">3D visualization will be available soon</p>
                                        <small class="text-muted">Will support interactive 3D molecular models</small>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Visualization Options -->
        <div class="row mt-4">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-cog"></i> Visualization Options</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-4">
                                <h6><i class="fas fa-palette text-primary"></i> Display Options</h6>
                                <div class="form-check">
                                    <input class="form-check-input" type="checkbox" id="showAtomLabels" checked>
                                    <label class="form-check-label" for="showAtomLabels">
                                        Show atom labels
                                    </label>
                                </div>
                                <div class="form-check">
                                    <input class="form-check-input" type="checkbox" id="showBondLabels">
                                    <label class="form-check-label" for="showBondLabels">
                                        Show bond labels
                                    </label>
                                </div>
                                <div class="form-check">
                                    <input class="form-check-input" type="checkbox" id="showHydrogens">
                                    <label class="form-check-label" for="showHydrogens">
                                        Show explicit hydrogens
                                    </label>
                                </div>
                            </div>
                            <div class="col-md-4">
                                <h6><i class="fas fa-download text-success"></i> Export Options</h6>
                                <button class="btn btn-outline-success btn-sm mb-2 w-100" onclick="exportStructure('png')">
                                    <i class="fas fa-image"></i> Export as PNG
                                </button>
                                <button class="btn btn-outline-success btn-sm mb-2 w-100" onclick="exportStructure('svg')">
                                    <i class="fas fa-vector-square"></i> Export as SVG
                                </button>
                                <button class="btn btn-outline-success btn-sm mb-2 w-100" onclick="exportStructure('mol')">
                                    <i class="fas fa-file-alt"></i> Export as MOL
                                </button>
                            </div>
                            <div class="col-md-4">
                                <h6><i class="fas fa-lightbulb text-warning"></i> Quick Examples</h6>
                                <button class="btn btn-outline-warning btn-sm mb-1 w-100" onclick="loadExample('caffeine')">
                                    Caffeine
                                </button>
                                <button class="btn btn-outline-warning btn-sm mb-1 w-100" onclick="loadExample('aspirin')">
                                    Aspirin
                                </button>
                                <button class="btn btn-outline-warning btn-sm mb-1 w-100" onclick="loadExample('tylenol')">
                                    Tylenol
                                </button>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Feature Notice -->
        <div class="row mt-4">
            <div class="col-lg-12">
                <div class="alert alert-info" role="alert">
                    <i class="fas fa-info-circle"></i> 
                    <strong>Note:</strong> The structure viewer is currently using a placeholder implementation. 
                    Full 2D/3D visualization will be implemented using libraries like RDKit.js or ChemDoodle.
                </div>
            </div>
        </div>
    </main>

    <?php include '../includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script src="../assets/js/main.js"></script>
    <script>
        // Structure viewer functionality
        document.getElementById('structureForm').addEventListener('submit', function(e) {
            e.preventDefault();
            
            const smiles = document.getElementById('smiles').value;
            if (!smiles) {
                ChemWeb.showToast('Please enter a SMILES notation', 'warning');
                return;
            }
            
            // Placeholder for structure visualization
            displayStructure(smiles);
        });

        function displayStructure(smiles) {
            const structure2d = document.getElementById('structure2d');
            
            // Placeholder implementation - would integrate with RDKit.js or similar
            structure2d.innerHTML = `
                <div class="text-center py-4">
                    <i class="fas fa-molecule fa-5x text-primary mb-3"></i>
                    <h5>Structure for: ${smiles}</h5>
                    <p class="text-muted">2D structure visualization would appear here</p>
                    <small class="text-muted">Integration with RDKit.js or ChemDoodle pending</small>
                </div>
            `;
            
            ChemWeb.showToast('Structure loaded (placeholder)', 'info');
        }

        function loadExample(molecule) {
            const examples = {
                'caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
                'aspirin': 'CC(=O)OC1=CC=CC=C1C(=O)O',
                'tylenol': 'CC(=O)NC1=CC=C(C=C1)O'
            };
            
            if (examples[molecule]) {
                document.getElementById('smiles').value = examples[molecule];
                displayStructure(examples[molecule]);
            }
        }

        function exportStructure(format) {
            ChemWeb.showToast(`Export as ${format.toUpperCase()} (feature coming soon)`, 'info');
        }

        // Initialize with caffeine example
        document.addEventListener('DOMContentLoaded', function() {
            displayStructure('CN1C=NC2=C1C(=O)N(C(=O)N2C)C');
        });
    </script>
</body>
</html>
