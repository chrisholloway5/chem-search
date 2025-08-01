<?php
session_start();
require_once '../includes/config.php';

// Handle SMILES to 3D coordinate conversion request
$moleculeData = null;
$errorMessage = null;
$smilesInput = '';

if ($_POST && isset($_POST['smiles']) && !empty(trim($_POST['smiles']))) {
    $smilesInput = trim($_POST['smiles']);
    
    // Generate 3D coordinates and molecular data
    $moleculeData = generateMoleculeData($smilesInput);
    if (!$moleculeData && !isset($moleculeData['error'])) {
        $errorMessage = "Unable to generate 3D coordinates for the provided SMILES notation.";
    }
}

function generateMoleculeData($smiles) {
    try {
        // In a real implementation, this would use RDKit Python to generate 3D coordinates
        // For now, we'll create placeholder data and use external services
        
        // Get basic molecular information from PubChem
        $pubchemData = getPubChemData($smiles);
        
        // Generate placeholder 3D coordinates (in real implementation, use RDKit)
        $atoms = generatePlaceholder3DCoordinates($smiles);
        
        return [
            'smiles' => $smiles,
            'atoms' => $atoms,
            'pubchem' => $pubchemData,
            'formula' => $pubchemData['MolecularFormula'] ?? 'Unknown',
            'weight' => $pubchemData['MolecularWeight'] ?? 0,
            'name' => $pubchemData['Title'] ?? $pubchemData['IUPACName'] ?? 'Unknown Compound'
        ];
        
    } catch (Exception $e) {
        return ['error' => $e->getMessage()];
    }
}

function getPubChemData($smiles) {
    $url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/" . urlencode($smiles) . "/property/Title,MolecularFormula,MolecularWeight,IUPACName/JSON";
    
    $context = stream_context_create([
        'http' => [
            'timeout' => 5,
            'user_agent' => 'ChemSearch/1.0'
        ]
    ]);
    
    $response = @file_get_contents($url, false, $context);
    
    if ($response !== false) {
        $data = json_decode($response, true);
        if (isset($data['PropertyTable']['Properties'][0])) {
            return $data['PropertyTable']['Properties'][0];
        }
    }
    
    return [];
}

function generatePlaceholder3DCoordinates($smiles) {
    // This is a placeholder - in real implementation, would use RDKit to generate actual 3D coordinates
    // For demonstration, we'll create some basic molecular geometries
    
    $commonMolecules = [
        'CCO' => [ // Ethanol
            ['C', 1.0, 0.0, 0.0],
            ['C', 0.0, 0.0, 0.0],
            ['O', -1.0, 0.0, 0.0],
            ['H', 1.5, 0.5, 0.5],
            ['H', 1.5, -0.5, -0.5],
            ['H', 1.5, 0.0, -1.0],
            ['H', 0.0, 1.0, 0.0],
            ['H', 0.0, -1.0, 0.0],
            ['H', -1.5, 0.0, 0.0]
        ],
        'O' => [ // Water
            ['O', 0.0, 0.0, 0.0],
            ['H', 0.76, 0.59, 0.0],
            ['H', -0.76, 0.59, 0.0]
        ],
        'CN1C=NC2=C1C(=O)N(C(=O)N2C)C' => [ // Caffeine (simplified)
            ['C', 0.0, 0.0, 0.0],
            ['N', 1.4, 0.0, 0.0],
            ['C', 2.1, 1.2, 0.0],
            ['N', 1.4, 2.4, 0.0],
            ['C', 0.0, 2.4, 0.0],
            ['C', -0.7, 1.2, 0.0],
            ['C', -2.1, 1.2, 0.0],
            ['O', -2.8, 2.4, 0.0],
            ['N', -2.8, 0.0, 0.0],
            ['C', -2.1, -1.2, 0.0],
            ['O', -2.8, -2.4, 0.0],
            ['N', -0.7, -1.2, 0.0],
            ['C', 2.1, 3.6, 0.0],
            ['C', -4.2, 0.0, 0.0],
            ['C', 0.0, -2.4, 0.0]
        ]
    ];
    
    if (isset($commonMolecules[$smiles])) {
        return $commonMolecules[$smiles];
    }
    
    // Generate random coordinates for unknown molecules
    $atomCount = min(20, max(5, strlen($smiles) / 3));
    $atoms = [];
    
    for ($i = 0; $i < $atomCount; $i++) {
        $element = ['C', 'N', 'O', 'H'][rand(0, 3)];
        $x = (rand(-100, 100) / 50.0);
        $y = (rand(-100, 100) / 50.0);
        $z = (rand(-100, 100) / 50.0);
        $atoms[] = [$element, $x, $y, $z];
    }
    
    return $atoms;
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>3D Structure Viewer - ChemSearch</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
    <link href="../assets/css/style.css" rel="stylesheet">
    
    <!-- 3D Molecular Visualization Libraries -->
    <script src="../assets/js/3Dmol.js"></script>
    
    <style>
        .molecule-viewer {
            border: 2px solid #dee2e6;
            border-radius: 8px;
            background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
            position: relative;
            overflow: hidden;
            height: 600px;
        }
        
        .viewer-controls {
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            border-bottom: 1px solid #dee2e6;
            padding: 10px;
            position: relative;
            z-index: 1000;
        }
        
        .control-group {
            margin-bottom: 10px;
        }
        
        .control-group:last-child {
            margin-bottom: 0;
        }
        
        .molecular-info {
            background: linear-gradient(135deg, #e3f2fd 0%, #bbdefb 100%);
            border-radius: 8px;
            padding: 20px;
            margin-top: 20px;
        }
        
        .info-item {
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 8px 0;
            border-bottom: 1px solid rgba(0,0,0,0.1);
        }
        
        .info-item:last-child {
            border-bottom: none;
        }
        
        .loading-spinner {
            position: absolute;
            top: 50%;
            left: 50%;
            transform: translate(-50%, -50%);
            z-index: 10;
        }
        
        #viewer-3d {
            width: 100%;
            height: 540px;
            min-height: 540px;
            background: radial-gradient(circle at center, #ffffff 0%, #f8f9fa 100%);
            border: 1px solid #dee2e6;
            border-radius: 8px;
            position: relative;
            overflow: hidden;
        }
        
        .feature-badge {
            display: inline-block;
            padding: 4px 8px;
            background: #28a745;
            color: white;
            border-radius: 12px;
            font-size: 0.75rem;
            margin: 2px;
        }
        
        .example-molecules {
            background: #f8f9fa;
            border-radius: 8px;
            padding: 15px;
            margin-top: 20px;
        }
        
        .example-btn {
            margin: 2px;
            font-size: 0.85rem;
        }
        
        .viewer-mode-tabs .nav-link {
            border-radius: 0;
            border: 1px solid #dee2e6;
            border-bottom: none;
            background: #f8f9fa;
        }
        
        .viewer-mode-tabs .nav-link.active {
            background: #007bff;
            color: white;
            border-color: #007bff;
        }
    </style>
</head>
<body>
    <?php include '../includes/header.php'; ?>
    
    <main class="container my-5">
        <div class="row">
            <div class="col-lg-12">
                <h1 class="mb-4">
                    <i class="fas fa-cube text-primary"></i> 3D Molecular Structure Viewer
                </h1>
                <p class="lead">Interactive 3D visualization of molecular structures with advanced controls.</p>
                
                <!-- Feature badges -->
                <div class="mb-4">
                    <span class="feature-badge"><i class="fas fa-cube"></i> 3D Visualization</span>
                    <span class="feature-badge"><i class="fas fa-palette"></i> Multiple Styles</span>
                    <span class="feature-badge"><i class="fas fa-expand-arrows-alt"></i> Interactive Controls</span>
                    <span class="feature-badge"><i class="fas fa-info-circle"></i> Molecular Properties</span>
                    <span class="feature-badge"><i class="fas fa-download"></i> Export Options</span>
                </div>
            </div>
        </div>

        <!-- Input Form -->
        <div class="row mb-4">
            <div class="col-lg-12">
                <div class="card shadow">
                    <div class="card-body">
                        <form method="POST" class="needs-validation" novalidate>
                            <div class="row">
                                <div class="col-md-8">
                                    <label for="smiles" class="form-label">
                                        <i class="fas fa-dna"></i> SMILES Notation
                                    </label>
                                    <input type="text" class="form-control" id="smiles" name="smiles" 
                                           value="<?php echo htmlspecialchars($smilesInput); ?>" 
                                           placeholder="e.g., CN1C=NC2=C1C(=O)N(C(=O)N2C)C (caffeine)"
                                           required>
                                    <div class="invalid-feedback">
                                        Please provide a valid SMILES notation.
                                    </div>
                                </div>
                                <div class="col-md-4 d-flex align-items-end">
                                    <button type="submit" class="btn btn-primary w-100">
                                        <i class="fas fa-search"></i> Generate 3D Structure
                                    </button>
                                </div>
                            </div>
                        </form>
                        
                        <!-- Example molecules -->
                        <div class="example-molecules">
                            <h6><i class="fas fa-flask"></i> Quick Examples:</h6>
                            <button type="button" class="btn btn-outline-secondary btn-sm example-btn" onclick="loadExample('O')">
                                Water (H₂O)
                            </button>
                            <button type="button" class="btn btn-outline-secondary btn-sm example-btn" onclick="loadExample('CCO')">
                                Ethanol
                            </button>
                            <button type="button" class="btn btn-outline-secondary btn-sm example-btn" onclick="loadExample('CC(=O)O')">
                                Acetic Acid
                            </button>
                            <button type="button" class="btn btn-outline-secondary btn-sm example-btn" onclick="loadExample('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')">
                                Caffeine
                            </button>
                            <button type="button" class="btn btn-outline-secondary btn-sm example-btn" onclick="loadExample('CC(C)CC1=CC=C(C=C1)C(C)C(=O)O')">
                                Ibuprofen
                            </button>
                            <button type="button" class="btn btn-outline-secondary btn-sm example-btn" onclick="loadExample('CC1=CC=C(C=C1)C(=O)O')">
                                p-Cresol
                            </button>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <?php if ($moleculeData && !isset($moleculeData['error'])): ?>
        <!-- 3D Viewer Section -->
        <div class="row">
            <div class="col-lg-8">
                <div class="card shadow">
                    <div class="card-header">
                        <ul class="nav nav-tabs viewer-mode-tabs" role="tablist">
                            <li class="nav-item">
                                <a class="nav-link active" id="3d-tab" data-bs-toggle="tab" href="#3d-viewer" role="tab">
                                    <i class="fas fa-cube"></i> 3D View
                                </a>
                            </li>
                            <li class="nav-item">
                                <a class="nav-link" id="info-tab" data-bs-toggle="tab" href="#molecule-info" role="tab">
                                    <i class="fas fa-info-circle"></i> Properties
                                </a>
                            </li>
                        </ul>
                    </div>
                    <div class="card-body p-0">
                        <div class="tab-content">
                            <div class="tab-pane fade show active" id="3d-viewer" role="tabpanel">
                                <div class="molecule-viewer">
                                    <!-- Viewer Controls -->
                                    <div class="viewer-controls">
                                        <div class="row">
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <label class="form-label"><i class="fas fa-palette"></i> Style:</label>
                                                    <select id="styleSelect" class="form-select form-select-sm">
                                                        <option value="sphere">Ball & Stick</option>
                                                        <option value="stick">Stick</option>
                                                        <option value="line">Wireframe</option>
                                                        <option value="cartoon">Cartoon</option>
                                                    </select>
                                                </div>
                                            </div>
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <label class="form-label"><i class="fas fa-eye"></i> Color:</label>
                                                    <select id="colorSelect" class="form-select form-select-sm">
                                                        <option value="element">By Element</option>
                                                        <option value="spectrum">Spectrum</option>
                                                        <option value="chain">By Chain</option>
                                                        <option value="residue">By Residue</option>
                                                    </select>
                                                </div>
                                            </div>
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <label class="form-label"><i class="fas fa-tag"></i> Labels:</label>
                                                    <div class="form-check form-switch">
                                                        <input class="form-check-input" type="checkbox" id="showLabels">
                                                        <label class="form-check-label" for="showLabels">Show</label>
                                                    </div>
                                                </div>
                                            </div>
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <button type="button" class="btn btn-outline-secondary btn-sm w-100" onclick="resetView()">
                                                        <i class="fas fa-redo"></i> Reset View
                                                    </button>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                    
                                    <!-- 3D Viewer Container -->
                                    <div id="viewer-3d">
                                        <div class="loading-spinner">
                                            <div class="spinner-border text-primary" role="status">
                                                <span class="visually-hidden">Loading 3D structure...</span>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                            <div class="tab-pane fade" id="molecule-info" role="tabpanel">
                                <div class="molecular-info">
                                    <h5><i class="fas fa-flask"></i> Molecular Properties</h5>
                                    <div class="info-item">
                                        <strong><i class="fas fa-tag"></i> Name:</strong>
                                        <span><?php echo htmlspecialchars($moleculeData['name']); ?></span>
                                    </div>
                                    <div class="info-item">
                                        <strong><i class="fas fa-dna"></i> SMILES:</strong>
                                        <span class="font-monospace"><?php echo htmlspecialchars($moleculeData['smiles']); ?></span>
                                    </div>
                                    <div class="info-item">
                                        <strong><i class="fas fa-atom"></i> Formula:</strong>
                                        <span><?php echo htmlspecialchars($moleculeData['formula']); ?></span>
                                    </div>
                                    <div class="info-item">
                                        <strong><i class="fas fa-weight"></i> Molecular Weight:</strong>
                                        <span><?php echo number_format($moleculeData['weight'], 2); ?> g/mol</span>
                                    </div>
                                    <div class="info-item">
                                        <strong><i class="fas fa-cube"></i> Atoms:</strong>
                                        <span><?php echo count($moleculeData['atoms']); ?> atoms</span>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            
            <!-- Control Panel -->
            <div class="col-lg-4">
                <div class="card shadow">
                    <div class="card-header">
                        <h5><i class="fas fa-sliders-h"></i> Visualization Controls</h5>
                    </div>
                    <div class="card-body">
                        <!-- Export Options -->
                        <div class="mb-3">
                            <h6><i class="fas fa-download"></i> Export Options:</h6>
                            <div class="d-grid gap-2">
                                <button type="button" class="btn btn-outline-primary btn-sm" onclick="exportPNG()">
                                    <i class="fas fa-image"></i> Export PNG
                                </button>
                                <button type="button" class="btn btn-outline-secondary btn-sm" onclick="exportPDB()">
                                    <i class="fas fa-file-code"></i> Export PDB
                                </button>
                            </div>
                        </div>
                        
                        <!-- Animation Controls -->
                        <div class="mb-3">
                            <h6><i class="fas fa-play"></i> Animation:</h6>
                            <div class="d-grid gap-2">
                                <button type="button" class="btn btn-outline-success btn-sm" onclick="startRotation()">
                                    <i class="fas fa-sync-alt"></i> Auto Rotate
                                </button>
                                <button type="button" class="btn btn-outline-danger btn-sm" onclick="stopRotation()">
                                    <i class="fas fa-stop"></i> Stop Rotation
                                </button>
                            </div>
                        </div>
                        
                        <!-- View Presets -->
                        <div class="mb-3">
                            <h6><i class="fas fa-camera"></i> View Presets:</h6>
                            <div class="d-grid gap-2">
                                <button type="button" class="btn btn-outline-info btn-sm" onclick="setView('front')">
                                    <i class="fas fa-eye"></i> Front View
                                </button>
                                <button type="button" class="btn btn-outline-info btn-sm" onclick="setView('side')">
                                    <i class="fas fa-eye"></i> Side View
                                </button>
                                <button type="button" class="btn btn-outline-info btn-sm" onclick="setView('top')">
                                    <i class="fas fa-eye"></i> Top View
                                </button>
                            </div>
                        </div>
                        
                        <!-- Instructions -->
                        <div class="alert alert-info">
                            <h6><i class="fas fa-info-circle"></i> Controls:</h6>
                            <small>
                                • <strong>Mouse:</strong> Drag to rotate<br>
                                • <strong>Scroll:</strong> Zoom in/out<br>
                                • <strong>Right-click:</strong> Pan view<br>
                                • <strong>Double-click:</strong> Center view
                            </small>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        
        <?php elseif (isset($moleculeData['error'])): ?>
        <div class="row">
            <div class="col-lg-12">
                <div class="alert alert-danger">
                    <h5><i class="fas fa-exclamation-triangle"></i> Error</h5>
                    <?php echo htmlspecialchars($moleculeData['error']); ?>
                </div>
            </div>
        </div>
        
        <?php elseif ($errorMessage): ?>
        <div class="row">
            <div class="col-lg-12">
                <div class="alert alert-warning">
                    <h5><i class="fas fa-exclamation-circle"></i> Warning</h5>
                    <?php echo htmlspecialchars($errorMessage); ?>
                </div>
            </div>
        </div>
        
        <?php else: ?>
        <!-- Welcome Section with Default 3D Viewer -->
        <div class="row">
            <div class="col-lg-8">
                <div class="card shadow">
                    <div class="card-header">
                        <ul class="nav nav-tabs viewer-mode-tabs" role="tablist">
                            <li class="nav-item">
                                <a class="nav-link active" id="3d-tab" data-bs-toggle="tab" href="#3d-viewer" role="tab">
                                    <i class="fas fa-cube"></i> 3D Demo
                                </a>
                            </li>
                            <li class="nav-item">
                                <a class="nav-link" id="welcome-tab" data-bs-toggle="tab" href="#welcome-info" role="tab">
                                    <i class="fas fa-info-circle"></i> Welcome
                                </a>
                            </li>
                        </ul>
                    </div>
                    <div class="card-body p-0">
                        <div class="tab-content">
                            <div class="tab-pane fade show active" id="3d-viewer" role="tabpanel">
                                <div class="molecule-viewer">
                                    <!-- Viewer Controls -->
                                    <div class="viewer-controls">
                                        <div class="row">
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <label class="form-label"><i class="fas fa-palette"></i> Style:</label>
                                                    <select id="styleSelect" class="form-select form-select-sm">
                                                        <option value="sphere">Ball & Stick</option>
                                                        <option value="stick">Stick</option>
                                                        <option value="line">Wireframe</option>
                                                        <option value="cartoon">Cartoon</option>
                                                    </select>
                                                </div>
                                            </div>
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <label class="form-label"><i class="fas fa-eye"></i> Color:</label>
                                                    <select id="colorSelect" class="form-select form-select-sm">
                                                        <option value="element">By Element</option>
                                                        <option value="spectrum">Spectrum</option>
                                                        <option value="chain">By Chain</option>
                                                        <option value="residue">By Residue</option>
                                                    </select>
                                                </div>
                                            </div>
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <label class="form-label"><i class="fas fa-tag"></i> Labels:</label>
                                                    <div class="form-check form-switch">
                                                        <input class="form-check-input" type="checkbox" id="showLabels">
                                                        <label class="form-check-label" for="showLabels">Show</label>
                                                    </div>
                                                </div>
                                            </div>
                                            <div class="col-md-3">
                                                <div class="control-group">
                                                    <button type="button" class="btn btn-outline-secondary btn-sm w-100" onclick="resetView()">
                                                        <i class="fas fa-redo"></i> Reset View
                                                    </button>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                    
                                    <!-- 3D Viewer Container -->
                                    <div id="viewer-3d">
                                        <div class="loading-spinner">
                                            <div class="spinner-border text-primary" role="status">
                                                <span class="visually-hidden">Loading 3D structure...</span>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                            <div class="tab-pane fade" id="welcome-info" role="tabpanel">
                                <div class="molecular-info text-center">
                                    <i class="fas fa-cube fa-4x text-primary mb-3"></i>
                                    <h3>Interactive 3D Molecular Visualization</h3>
                                    <p class="lead">Enter a SMILES notation above to generate and visualize molecular structures in 3D.</p>
                                    
                                    <div class="row mt-4">
                                        <div class="col-md-4">
                                            <div class="feature-card p-3">
                                                <i class="fas fa-cube fa-2x text-primary mb-2"></i>
                                                <h5>3D Visualization</h5>
                                                <p>Interactive 3D molecular models with multiple rendering styles</p>
                                            </div>
                                        </div>
                                        <div class="col-md-4">
                                            <div class="feature-card p-3">
                                                <i class="fas fa-palette fa-2x text-success mb-2"></i>
                                                <h5>Multiple Styles</h5>
                                                <p>Ball & stick, wireframe, cartoon, and more visualization modes</p>
                                            </div>
                                        </div>
                                        <div class="col-md-4">
                                            <div class="feature-card p-3">
                                                <i class="fas fa-download fa-2x text-info mb-2"></i>
                                                <h5>Export Options</h5>
                                                <p>Export structures as PNG images or PDB files</p>
                                            </div>
                                        </div>
                                    </div>
                                    
                                    <div class="alert alert-info mt-4">
                                        <i class="fas fa-lightbulb"></i>
                                        <strong>Demo:</strong> A sample caffeine molecule is loaded in the 3D viewer. Try the controls above to explore different visualization styles!
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            
            <!-- Control Panel -->
            <div class="col-lg-4">
                <div class="card shadow">
                    <div class="card-header">
                        <h5><i class="fas fa-sliders-h"></i> Visualization Controls</h5>
                    </div>
                    <div class="card-body">
                        <!-- Export Options -->
                        <div class="mb-3">
                            <h6><i class="fas fa-download"></i> Export Options:</h6>
                            <div class="d-grid gap-2">
                                <button type="button" class="btn btn-outline-primary btn-sm" onclick="exportPNG()">
                                    <i class="fas fa-image"></i> Export PNG
                                </button>
                                <button type="button" class="btn btn-outline-secondary btn-sm" onclick="exportPDB()">
                                    <i class="fas fa-file-code"></i> Export PDB
                                </button>
                            </div>
                        </div>
                        
                        <!-- Animation Controls -->
                        <div class="mb-3">
                            <h6><i class="fas fa-play"></i> Animation:</h6>
                            <div class="d-grid gap-2">
                                <button type="button" class="btn btn-outline-success btn-sm" onclick="startRotation()">
                                    <i class="fas fa-sync-alt"></i> Auto Rotate
                                </button>
                                <button type="button" class="btn btn-outline-danger btn-sm" onclick="stopRotation()">
                                    <i class="fas fa-stop"></i> Stop Rotation
                                </button>
                            </div>
                        </div>
                        
                        <!-- View Presets -->
                        <div class="mb-3">
                            <h6><i class="fas fa-camera"></i> View Presets:</h6>
                            <div class="d-grid gap-2">
                                <button type="button" class="btn btn-outline-info btn-sm" onclick="setView('front')">
                                    <i class="fas fa-eye"></i> Front View
                                </button>
                                <button type="button" class="btn btn-outline-info btn-sm" onclick="setView('side')">
                                    <i class="fas fa-eye"></i> Side View
                                </button>
                                <button type="button" class="btn btn-outline-info btn-sm" onclick="setView('top')">
                                    <i class="fas fa-eye"></i> Top View
                                </button>
                            </div>
                        </div>
                        
                        <!-- Instructions -->
                        <div class="alert alert-info">
                            <h6><i class="fas fa-info-circle"></i> Controls:</h6>
                            <small>
                                • <strong>Mouse:</strong> Drag to rotate<br>
                                • <strong>Scroll:</strong> Zoom in/out<br>
                                • <strong>Right-click:</strong> Pan view<br>
                                • <strong>Double-click:</strong> Center view
                            </small>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        <?php endif; ?>
    </main>

    <?php include '../includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
    <script src="../assets/js/main.js"></script>
    
    <script>
        // Global 3DMol viewer instance
        let viewer;
        let isRotating = false;
        let rotationInterval;
        
        // Initialize 3D viewer when document is ready
        document.addEventListener('DOMContentLoaded', function() {
            // Setup form validation first
            setupFormValidation();
            
            // Wait for all resources to load before initializing viewer
            if (document.readyState === 'complete') {
                initViewer();
            } else {
                window.addEventListener('load', initViewer);
            }
        });

        function initViewer() {
            // Check if 3DMol is available
            if (typeof $3Dmol === 'undefined') {
                console.error('3DMol.js library not loaded');
                showAlert('3D visualization library failed to load. Please refresh the page.', 'danger');
                return;
            }
            
            // Wait a moment for DOM to settle, then initialize
            setTimeout(() => {
                const viewerContainer = document.getElementById('viewer-3d');
                if (viewerContainer) {
                    // Ensure container has dimensions
                    const containerRect = viewerContainer.getBoundingClientRect();
                    if (containerRect.width === 0 || containerRect.height === 0) {
                        console.warn('Viewer container has no dimensions, retrying in 100ms');
                        setTimeout(initViewer, 100);
                        return;
                    }
                    
                    initializeViewer();
                    
                    <?php if ($moleculeData && !isset($moleculeData['error'])): ?>
                    // Load molecule data if available
                    loadMoleculeData();
                    <?php else: ?>
                    // Load a default molecule for demonstration
                    loadDefaultMolecule();
                    <?php endif; ?>
                }
            }, 100);
        }
        
        function initializeViewer() {
            try {
                console.log('Initializing 3D viewer...');
                
                // Check for WebGL support
                if (!window.WebGLRenderingContext) {
                    throw new Error('WebGL is not supported by this browser');
                }
                
                // Test WebGL context
                const canvas = document.createElement('canvas');
                const gl = canvas.getContext('webgl') || canvas.getContext('experimental-webgl');
                if (!gl) {
                    throw new Error('WebGL context could not be created');
                }
                
                // Initialize 3DMol.js viewer
                const element = document.getElementById('viewer-3d');
                if (!element) {
                    console.error('viewer-3d element not found');
                    return;
                }
                
                // Get container dimensions
                const containerRect = element.getBoundingClientRect();
                const width = Math.max(containerRect.width, 400);
                const height = 540;
                
                // Ensure element has explicit dimensions
                element.style.width = width + 'px';
                element.style.height = height + 'px';
                
                const config = { 
                    backgroundColor: 'white',
                    width: width,
                    height: height,
                    antialias: true,
                    camerax: 0,
                    cameray: 0,
                    cameraz: 50
                };
                
                console.log('Creating viewer with config:', config);
                viewer = $3Dmol.createViewer(element, config);
                
                // Wait for viewer to be ready
                if (viewer) {
                    console.log('3D viewer created successfully');
                    
                    // Force a resize to ensure proper canvas dimensions
                    setTimeout(() => {
                        if (viewer.resize) {
                            viewer.resize();
                        }
                    }, 50);
                } else {
                    throw new Error('Failed to create 3DMol viewer');
                }
                
                // Hide loading spinner
                const spinner = element.querySelector('.loading-spinner');
                if (spinner) {
                    spinner.style.display = 'none';
                }
                
                // Add window resize handler
                window.addEventListener('resize', function() {
                    if (viewer && viewer.resize) {
                        setTimeout(() => {
                            const newWidth = Math.max(element.getBoundingClientRect().width, 400);
                            element.style.width = newWidth + 'px';
                            viewer.resize();
                        }, 100);
                    }
                });
                
                return true;
            } catch (error) {
                console.error('Error initializing 3D viewer:', error);
                showAlert('Failed to initialize 3D viewer: ' + error.message, 'danger');
                return false;
            }
        }
        
        function loadMoleculeData() {
            if (!viewer) {
                console.error('Viewer not initialized');
                return;
            }
            
            try {
                console.log('Loading molecule data...');
                
                // Get molecular data from PHP
                const moleculeData = <?php echo json_encode($moleculeData); ?>;
                
                // Clear existing models
                viewer.clear();
                
                // Create molecule from atom coordinates
                if (moleculeData.atoms && moleculeData.atoms.length > 0) {
                    // Convert atom data to XYZ format
                    let xyzData = moleculeData.atoms.length + "\n";
                    xyzData += "Generated from " + (moleculeData.smiles || "unknown") + "\n";
                    
                    moleculeData.atoms.forEach(atom => {
                        const [element, x, y, z] = atom;
                        xyzData += element + " " + x + " " + y + " " + z + "\n";
                    });
                    
                    viewer.addModel(xyzData, "xyz");
                } else {
                    showAlert('No molecular coordinate data found.', 'danger');
                    return;
                }
                
                setStyle('sphere');
                viewer.zoomTo();
                viewer.render();
                
                console.log('Molecule data loaded successfully');
                
            } catch (error) {
                console.error('Error loading molecule data:', error);
                showAlert('Failed to load molecule data: ' + error.message, 'danger');
            }
        }
        
        function loadDefaultMolecule() {
            if (!viewer) {
                console.error('Viewer not initialized');
                return;
            }
            
            try {
                console.log('Loading default molecule (caffeine)...');
                
                // Clear existing models
                viewer.clear();
                
                // Create caffeine molecule using XYZ format
                const caffeineXYZ = `15
Caffeine molecule
C 0.0 0.0 0.0
N 1.4 0.0 0.0
C 2.1 1.2 0.0
N 1.4 2.4 0.0
C 0.0 2.4 0.0
C -0.7 1.2 0.0
C -2.1 1.2 0.0
O -2.8 2.4 0.0
N -2.8 0.0 0.0
C -2.1 -1.2 0.0
O -2.8 -2.4 0.0
N -0.7 -1.2 0.0
C 2.1 3.6 0.0
C -4.2 0.0 0.0
C 0.0 -2.4 0.0`;
                
                viewer.addModel(caffeineXYZ, "xyz");
                
                setStyle('sphere');
                viewer.zoomTo();
                viewer.render();
                
                console.log('Default molecule loaded successfully');
                
            } catch (error) {
                console.error('Error loading default molecule:', error);
                showAlert('Failed to load default molecule: ' + error.message, 'danger');
            }
        }
        
        function loadExampleMolecule(smiles) {
            if (!viewer) {
                console.error('Viewer not initialized');
                return;
            }
            
            try {
                console.log('Loading example molecule:', smiles);
                
                // Clear existing models
                viewer.clear();
                
                // Define example molecules in XYZ format
                const molecules = {
                    'O': `3
Water molecule
O 0.0 0.0 0.0
H 0.76 0.59 0.0
H -0.76 0.59 0.0`,
                    'CCO': `9
Ethanol molecule
C 1.0 0.0 0.0
C 0.0 0.0 0.0
O -1.0 0.0 0.0
H 1.5 0.5 0.5
H 1.5 -0.5 -0.5
H 1.5 0.0 -1.0
H 0.0 1.0 0.0
H 0.0 -1.0 0.0
H -1.5 0.0 0.0`,
                    'CC(=O)O': `8
Acetic acid molecule
C 1.0 0.0 0.0
C 0.0 0.0 0.0
O -0.5 1.0 0.0
O -0.5 -1.0 0.0
H 1.5 0.5 0.5
H 1.5 -0.5 -0.5
H 1.5 0.0 -1.0
H -1.0 -1.5 0.0`,
                    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': `15
Caffeine molecule
C 0.0 0.0 0.0
N 1.4 0.0 0.0
C 2.1 1.2 0.0
N 1.4 2.4 0.0
C 0.0 2.4 0.0
C -0.7 1.2 0.0
C -2.1 1.2 0.0
O -2.8 2.4 0.0
N -2.8 0.0 0.0
C -2.1 -1.2 0.0
O -2.8 -2.4 0.0
N -0.7 -1.2 0.0
C 2.1 3.6 0.0
C -4.2 0.0 0.0
C 0.0 -2.4 0.0`
                };
                
                const moleculeXYZ = molecules[smiles];
                if (!moleculeXYZ) {
                    console.warn('Molecule not found in predefined set, loading default');
                    loadDefaultMolecule();
                    return;
                }
                
                viewer.addModel(moleculeXYZ, "xyz");
                
                setStyle('sphere');
                viewer.zoomTo();
                viewer.render();
                
                console.log('Example molecule loaded successfully');
                showAlert('Molecule loaded: ' + smiles, 'success');
                
            } catch (error) {
                console.error('Error loading example molecule:', error);
                showAlert('Failed to load example molecule: ' + error.message, 'danger');
            }
        }
        
        function calculateDistance(atom1, atom2) {
            const dx = atom1[1] - atom2[1];
            const dy = atom1[2] - atom2[2];
            const dz = atom1[3] - atom2[3];
            return Math.sqrt(dx*dx + dy*dy + dz*dz);
        }
        
        function setStyle(style) {
            if (!viewer) {
                console.error('Viewer not initialized');
                return;
            }
            
            try {
                viewer.setStyle({}, {}); // Clear existing styles
                
                switch(style) {
                    case 'sphere':
                        viewer.setStyle({}, {
                            sphere: {radius: 0.3},
                            stick: {radius: 0.1}
                        });
                        break;
                    case 'stick':
                        viewer.setStyle({}, {
                            stick: {radius: 0.15}
                        });
                        break;
                    case 'line':
                        viewer.setStyle({}, {
                            line: {linewidth: 2}
                        });
                        break;
                    case 'cartoon':
                        viewer.setStyle({}, {
                            cartoon: {color: 'spectrum'}
                        });
                        break;
                    default:
                        viewer.setStyle({}, {
                            sphere: {radius: 0.3},
                            stick: {radius: 0.1}
                        });
                }
                
                viewer.render();
                console.log('Style changed to:', style);
                
            } catch (error) {
                console.error('Error setting style:', error);
                showAlert('Failed to change visualization style: ' + error.message, 'warning');
            }
        }
        
        function setColorScheme(scheme) {
            if (!viewer) return;
            
            const currentStyle = document.getElementById('styleSelect').value;
            
            switch(scheme) {
                case 'element':
                    setStyle(currentStyle); // Default coloring by element
                    break;
                case 'spectrum':
                    viewer.setStyle({}, {
                        [currentStyle]: {colorscheme: 'spectrum'}
                    });
                    break;
                case 'chain':
                    viewer.setStyle({}, {
                        [currentStyle]: {colorscheme: 'chain'}
                    });
                    break;
                case 'residue':
                    viewer.setStyle({}, {
                        [currentStyle]: {colorscheme: 'residue'}
                    });
                    break;
            }
            
            viewer.render();
        }
        
        function toggleLabels(show) {
            if (!viewer) return;
            
            try {
                if (show) {
                    // Add atom labels using 3DMol.js API
                    const atoms = viewer.getModel().selectedAtoms({});
                    atoms.forEach((atom, index) => {
                        viewer.addLabel(atom.elem, {
                            position: atom,
                            fontSize: 12,
                            fontColor: 'black',
                            backgroundColor: 'white',
                            backgroundOpacity: 0.8
                        });
                    });
                } else {
                    // Remove all labels
                    viewer.removeAllLabels();
                }
                
                viewer.render();
            } catch (error) {
                console.warn('Label functionality not fully supported:', error);
                // Silently fail - labels are optional
            }
        }
        
        function resetView() {
            if (!viewer) return;
            
            viewer.zoomTo();
            viewer.render();
            stopRotation();
        }
        
        function startRotation() {
            if (!viewer || isRotating) return;
            
            isRotating = true;
            rotationInterval = setInterval(() => {
                viewer.rotate(1, 'y');
                viewer.render();
            }, 50);
        }
        
        function stopRotation() {
            if (rotationInterval) {
                clearInterval(rotationInterval);
                rotationInterval = null;
            }
            isRotating = false;
        }
        
        function setView(viewType) {
            if (!viewer) return;
            
            switch(viewType) {
                case 'front':
                    viewer.rotate(0, 'x');
                    viewer.rotate(0, 'y');
                    break;
                case 'side':
                    viewer.rotate(90, 'y');
                    break;
                case 'top':
                    viewer.rotate(90, 'x');
                    break;
            }
            
            viewer.render();
        }
        
        function exportPNG() {
            if (!viewer) {
                showAlert('3D viewer not initialized', 'warning');
                return;
            }
            
            try {
                viewer.render(); // Ensure rendering is complete
                viewer.pngURI(null, null, null, null, function(uri) {
                    const link = document.createElement('a');
                    link.download = 'molecule_3d.png';
                    link.href = uri;
                    link.click();
                    showAlert('PNG export successful', 'success');
                });
            } catch (error) {
                showAlert('PNG export failed: ' + error.message, 'danger');
            }
        }
        
        function exportPDB() {
            showAlert('PDB export feature coming soon', 'info');
        }
        
        function setupFormValidation() {
            // Form validation
            const forms = document.querySelectorAll('.needs-validation');
            Array.prototype.slice.call(forms).forEach(function(form) {
                form.addEventListener('submit', function(event) {
                    if (!form.checkValidity()) {
                        event.preventDefault();
                        event.stopPropagation();
                    }
                    form.classList.add('was-validated');
                }, false);
            });
        }
        
        function loadExample(smiles) {
            document.getElementById('smiles').value = smiles;
            // Immediately load the molecule in the 3D viewer
            loadExampleMolecule(smiles);
        }
        
        function loadExampleMolecule(smiles) {
            if (!viewer) {
                console.error('Viewer not initialized');
                return;
            }
            
            try {
                console.log('Loading example molecule:', smiles);
                
                // Clear existing models
                viewer.clear();
                
                // Define example molecules in XYZ format
                const molecules = {
                    'O': `3
Water molecule
O 0.0 0.0 0.0
H 0.76 0.59 0.0
H -0.76 0.59 0.0`,
                    'CCO': `9
Ethanol molecule
C 1.0 0.0 0.0
C 0.0 0.0 0.0
O -1.0 0.0 0.0
H 1.5 0.5 0.5
H 1.5 -0.5 -0.5
H 1.5 0.0 -1.0
H 0.0 1.0 0.0
H 0.0 -1.0 0.0
H -1.5 0.0 0.0`,
                    'CC(=O)O': `8
Acetic acid molecule
C 1.0 0.0 0.0
C 0.0 0.0 0.0
O -0.5 1.0 0.0
O -0.5 -1.0 0.0
H 1.5 0.5 0.5
H 1.5 -0.5 -0.5
H 1.5 0.0 -1.0
H -1.0 -1.5 0.0`,
                    'CN1C=NC2=C1C(=O)N(C(=O)N2C)C': `15
Caffeine molecule
C 0.0 0.0 0.0
N 1.4 0.0 0.0
C 2.1 1.2 0.0
N 1.4 2.4 0.0
C 0.0 2.4 0.0
C -0.7 1.2 0.0
C -2.1 1.2 0.0
O -2.8 2.4 0.0
N -2.8 0.0 0.0
C -2.1 -1.2 0.0
O -2.8 -2.4 0.0
N -0.7 -1.2 0.0
C 2.1 3.6 0.0
C -4.2 0.0 0.0
C 0.0 -2.4 0.0`
                };
                
                const moleculeXYZ = molecules[smiles];
                if (!moleculeXYZ) {
                    console.warn('Molecule not found in predefined set, loading default');
                    loadDefaultMolecule();
                    return;
                }
                
                viewer.addModel(moleculeXYZ, "xyz");
                
                setStyle('sphere');
                viewer.zoomTo();
                viewer.render();
                
                console.log('Example molecule loaded successfully');
                showAlert('Molecule loaded: ' + smiles, 'success');
                
            } catch (error) {
                console.error('Error loading example molecule:', error);
                showAlert('Failed to load example molecule: ' + error.message, 'danger');
            }
        }
        
        function showAlert(message, type) {
            // Create and show bootstrap alert
            const alertDiv = document.createElement('div');
            alertDiv.className = `alert alert-${type} alert-dismissible fade show`;
            alertDiv.innerHTML = `
                ${message}
                <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
            `;
            
            // Insert at top of main container
            const main = document.querySelector('main');
            main.insertBefore(alertDiv, main.firstChild);
            
            // Auto-dismiss after 3 seconds
            setTimeout(() => {
                alertDiv.remove();
            }, 3000);
        }
        
        // Event listeners for controls
        document.addEventListener('DOMContentLoaded', function() {
            // Style selector
            const styleSelect = document.getElementById('styleSelect');
            if (styleSelect) {
                styleSelect.addEventListener('change', function() {
                    setStyle(this.value);
                });
            }
            
            // Color selector
            const colorSelect = document.getElementById('colorSelect');
            if (colorSelect) {
                colorSelect.addEventListener('change', function() {
                    setColorScheme(this.value);
                });
            }
            
            // Labels toggle
            const showLabels = document.getElementById('showLabels');
            if (showLabels) {
                showLabels.addEventListener('change', function() {
                    toggleLabels(this.checked);
                });
            }
        });
    </script>
</body>
</html>
