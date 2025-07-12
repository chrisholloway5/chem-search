<?php
session_start();
require_once '../includes/config.php';

// Handle search request
$searchResults = null;
$errorMessage = null;
$searchQuery = '';
$searchType = 'name';

if ($_POST && isset($_POST['search_query']) && !empty(trim($_POST['search_query']))) {
    $searchQuery = trim($_POST['search_query']);
    $searchType = $_POST['search_type'] ?? 'name';
    
    // Perform database search
    $searchResults = performDatabaseSearch($searchQuery, $searchType);
}

function performDatabaseSearch($query, $type) {
    $results = [];
    
    try {
        // PubChem API search
        $pubchemData = searchPubChem($query, $type);
        if ($pubchemData) {
            $results['pubchem'] = $pubchemData;
        }
        
        // If no results from PubChem, try alternative search methods
        if (empty($results)) {
            $results['message'] = "No results found in PubChem. Try a different search term or method.";
        }
        
    } catch (Exception $e) {
        $results['error'] = "Search error: " . $e->getMessage();
    }
    
    return $results;
}

function searchPubChem($query, $type) {
    $baseUrl = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";
    
    // Determine search endpoint based on type
    switch ($type) {
        case 'name':
            $searchUrl = "$baseUrl/compound/name/" . urlencode($query) . "/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON";
            break;
        case 'smiles':
            $searchUrl = "$baseUrl/compound/smiles/" . urlencode($query) . "/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON";
            break;
        case 'inchi':
            $searchUrl = "$baseUrl/compound/inchi/" . urlencode($query) . "/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON";
            break;
        case 'formula':
            $searchUrl = "$baseUrl/compound/formula/" . urlencode($query) . "/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON";
            break;
        case 'cas':
            // CAS numbers are searched as names in PubChem
            $searchUrl = "$baseUrl/compound/name/" . urlencode($query) . "/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON";
            break;
        default:
            $searchUrl = "$baseUrl/compound/name/" . urlencode($query) . "/property/Title,MolecularFormula,MolecularWeight,CanonicalSMILES,InChI,InChIKey,IUPACName/JSON";
    }
    
    // Make API request
    $context = stream_context_create([
        'http' => [
            'timeout' => 10,
            'user_agent' => 'ChemSearch/1.0'
        ]
    ]);
    
    $response = @file_get_contents($searchUrl, false, $context);
    
    if ($response === false) {
        return null;
    }
    
    $data = json_decode($response, true);
    
    if (isset($data['PropertyTable']['Properties'])) {
        return $data['PropertyTable']['Properties'];
    }
    
    return null;
}

function getSafetyInfo($cid) {
    // This would typically fetch safety information from PubChem
    // For now, return placeholder data
    return [
        'hazards' => ['May cause skin irritation', 'Avoid inhalation'],
        'precautions' => ['Wear protective gloves', 'Use in well-ventilated area'],
        'ghs_classification' => 'Not available'
    ];
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Database Search - Chem Search</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <link href="../assets/css/style.css" rel="stylesheet">
</head>
<body>
    <?php include '../includes/header.php'; ?>
    
    <main class="container my-5">
        <!-- Header -->
        <div class="row mb-4">
            <div class="col-lg-12">
                <h1 class="mb-3">
                    <i class="fas fa-search text-info"></i> Chemical Database Search
                </h1>
                <p class="lead">Search chemical databases like PubChem for comprehensive compound information.</p>
            </div>
        </div>

        <!-- Search Form -->
        <div class="row mb-4">
            <div class="col-lg-12">
                <div class="card shadow">
                    <div class="card-header bg-info text-white">
                        <h5><i class="fas fa-database"></i> Search Chemical Databases</h5>
                    </div>
                    <div class="card-body">
                        <form method="POST" class="mb-3">
                            <div class="row g-3">
                                <div class="col-md-6">
                                    <label for="search_query" class="form-label">Search Query</label>
                                    <input type="text" class="form-control" id="search_query" name="search_query" 
                                           value="<?php echo htmlspecialchars($searchQuery); ?>" 
                                           placeholder="Enter compound name, SMILES, formula, etc." required>
                                </div>
                                <div class="col-md-4">
                                    <label for="search_type" class="form-label">Search Type</label>
                                    <select class="form-select" id="search_type" name="search_type">
                                        <option value="name" <?php echo $searchType === 'name' ? 'selected' : ''; ?>>Compound Name</option>
                                        <option value="smiles" <?php echo $searchType === 'smiles' ? 'selected' : ''; ?>>SMILES Notation</option>
                                        <option value="inchi" <?php echo $searchType === 'inchi' ? 'selected' : ''; ?>>InChI/InChI Key</option>
                                        <option value="formula" <?php echo $searchType === 'formula' ? 'selected' : ''; ?>>Molecular Formula</option>
                                        <option value="cas" <?php echo $searchType === 'cas' ? 'selected' : ''; ?>>CAS Number</option>
                                    </select>
                                </div>
                                <div class="col-md-2 d-flex align-items-end">
                                    <button type="submit" class="btn btn-info w-100">
                                        <i class="fas fa-search"></i> Search
                                    </button>
                                </div>
                            </div>
                        </form>
                        
                        <!-- Search Examples -->
                        <div class="row">
                            <div class="col-12">
                                <small class="text-muted">
                                    <strong>Examples:</strong> 
                                    <span class="badge bg-light text-dark me-1">aspirin</span>
                                    <span class="badge bg-light text-dark me-1">CCO</span>
                                    <span class="badge bg-light text-dark me-1">C8H9NO2</span>
                                    <span class="badge bg-light text-dark me-1">50-78-2</span>
                                </small>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Search Results -->
        <?php if ($searchResults !== null): ?>
        <div class="row">
            <div class="col-lg-12">
                <?php if (isset($searchResults['error'])): ?>
                    <div class="alert alert-danger">
                        <i class="fas fa-exclamation-triangle"></i> <?php echo htmlspecialchars($searchResults['error']); ?>
                    </div>
                <?php elseif (isset($searchResults['message'])): ?>
                    <div class="alert alert-warning">
                        <i class="fas fa-info-circle"></i> <?php echo htmlspecialchars($searchResults['message']); ?>
                    </div>
                <?php elseif (isset($searchResults['pubchem']) && is_array($searchResults['pubchem'])): ?>
                    <div class="card shadow">
                        <div class="card-header bg-success text-white">
                            <h5><i class="fas fa-check-circle"></i> Search Results from PubChem</h5>
                        </div>
                        <div class="card-body">
                            <?php 
                            $compounds = $searchResults['pubchem'];
                            $compoundCount = count($compounds);
                            ?>
                            <p class="mb-3"><strong>Found <?php echo $compoundCount; ?> compound(s)</strong></p>
                            
                            <?php foreach ($compounds as $index => $compound): ?>
                            <div class="compound-result border rounded p-3 mb-3 <?php echo $index % 2 === 0 ? 'bg-light' : ''; ?>">
                                <div class="row">
                                    <div class="col-md-8">
                                        <h6 class="text-primary">
                                            <i class="fas fa-flask"></i> 
                                            <?php echo htmlspecialchars($compound['Title'] ?? $compound['IUPACName'] ?? 'Unknown Compound'); ?>
                                        </h6>
                                        
                                        <div class="row g-2 mt-2">
                                            <div class="col-sm-6">
                                                <small class="text-muted">Molecular Formula:</small><br>
                                                <code class="bg-primary text-white px-2 py-1 rounded"><?php echo htmlspecialchars($compound['MolecularFormula'] ?? 'N/A'); ?></code>
                                            </div>
                                            <div class="col-sm-6">
                                                <small class="text-muted">Molecular Weight:</small><br>
                                                <strong><?php echo number_format($compound['MolecularWeight'] ?? 0, 2); ?> g/mol</strong>
                                            </div>
                                        </div>
                                        
                                        <?php if (!empty($compound['CanonicalSMILES'])): ?>
                                        <div class="mt-2">
                                            <small class="text-muted">SMILES:</small><br>
                                            <code class="bg-success text-white px-2 py-1 rounded small"><?php echo htmlspecialchars($compound['CanonicalSMILES']); ?></code>
                                        </div>
                                        <?php endif; ?>
                                        
                                        <?php if (!empty($compound['InChIKey'])): ?>
                                        <div class="mt-2">
                                            <small class="text-muted">InChI Key:</small><br>
                                            <code class="bg-info text-white px-2 py-1 rounded small"><?php echo htmlspecialchars($compound['InChIKey']); ?></code>
                                        </div>
                                        <?php endif; ?>
                                    </div>
                                    <div class="col-md-4">
                                        <div class="text-center">
                                            <div class="bg-white border rounded p-2 mb-2">
                                                <i class="fas fa-atom fa-3x text-info"></i>
                                                <br><small class="text-muted">Structure Placeholder</small>
                                            </div>
                                            <div class="d-grid gap-1">
                                                <a href="molecule_analyzer.php?smiles=<?php echo urlencode($compound['CanonicalSMILES'] ?? ''); ?>" 
                                                   class="btn btn-sm btn-primary">
                                                    <i class="fas fa-calculator"></i> Analyze
                                                </a>
                                                <button class="btn btn-sm btn-outline-secondary" onclick="copyToClipboard('<?php echo addslashes($compound['CanonicalSMILES'] ?? ''); ?>')">
                                                    <i class="fas fa-copy"></i> Copy SMILES
                                                </button>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                                
                                <!-- Additional Information -->
                                <div class="row mt-3">
                                    <div class="col-12">
                                        <div class="accordion" id="accordion<?php echo $index; ?>">
                                            <div class="accordion-item">
                                                <h2 class="accordion-header">
                                                    <button class="accordion-button collapsed" type="button" data-bs-toggle="collapse" 
                                                            data-bs-target="#collapse<?php echo $index; ?>">
                                                        <i class="fas fa-info-circle me-2"></i> Additional Information
                                                    </button>
                                                </h2>
                                                <div id="collapse<?php echo $index; ?>" class="accordion-collapse collapse" 
                                                     data-bs-parent="#accordion<?php echo $index; ?>">
                                                    <div class="accordion-body">
                                                        <div class="row">
                                                            <div class="col-md-6">
                                                                <h6><i class="fas fa-tag text-primary"></i> Identifiers</h6>
                                                                <table class="table table-sm">
                                                                    <tr><td>IUPAC Name:</td><td><?php echo htmlspecialchars($compound['IUPACName'] ?? 'N/A'); ?></td></tr>
                                                                    <tr><td>InChI:</td><td class="small"><?php echo htmlspecialchars(substr($compound['InChI'] ?? 'N/A', 0, 50)) . (strlen($compound['InChI'] ?? '') > 50 ? '...' : ''); ?></td></tr>
                                                                </table>
                                                            </div>
                                                            <div class="col-md-6">
                                                                <h6><i class="fas fa-shield-alt text-warning"></i> Safety Information</h6>
                                                                <div class="alert alert-info small">
                                                                    <i class="fas fa-info-circle"></i> 
                                                                    For detailed safety information, consult the official PubChem page or safety data sheets.
                                                                </div>
                                                            </div>
                                                        </div>
                                                    </div>
                                                </div>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                            </div>
                            <?php endforeach; ?>
                        </div>
                    </div>
                <?php endif; ?>
            </div>
        </div>
        <?php endif; ?>

        <!-- Database Information -->
        <div class="row mt-4">
            <div class="col-lg-12">
                <div class="card border-secondary">
                    <div class="card-header bg-secondary text-white">
                        <h5><i class="fas fa-database"></i> Database Information</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h6><i class="fas fa-check text-success"></i> Available Features:</h6>
                                <ul class="list-unstyled">
                                    <li><i class="fas fa-check-circle text-success"></i> PubChem compound search</li>
                                    <li><i class="fas fa-check-circle text-success"></i> Chemical property lookup</li>
                                    <li><i class="fas fa-check-circle text-success"></i> Multiple search methods</li>
                                    <li><i class="fas fa-check-circle text-success"></i> Molecular information</li>
                                    <li><i class="fas fa-check-circle text-success"></i> Structure identifiers</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-road text-info"></i> Coming Soon:</h6>
                                <ul class="list-unstyled">
                                    <li><i class="fas fa-clock text-muted"></i> ChEMBL database integration</li>
                                    <li><i class="fas fa-clock text-muted"></i> Biological activity data</li>
                                    <li><i class="fas fa-clock text-muted"></i> Toxicity information</li>
                                    <li><i class="fas fa-clock text-muted"></i> 3D structure visualization</li>
                                    <li><i class="fas fa-clock text-muted"></i> Batch search capabilities</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Related Tools -->
        <div class="row mt-4">
            <div class="col-lg-12">
                <div class="card border-primary">
                    <div class="card-header bg-primary text-white">
                        <h5><i class="fas fa-tools"></i> Related Tools</h5>
                    </div>
                    <div class="card-body">
                        <div class="row g-3">
                            <div class="col-md-4">
                                <div class="d-grid">
                                    <a href="molecule_analyzer.php" class="btn btn-outline-primary">
                                        <i class="fas fa-atom"></i> Molecule Analyzer
                                        <br><small>Analyze molecular properties</small>
                                    </a>
                                </div>
                            </div>
                            <div class="col-md-4">
                                <div class="d-grid">
                                    <a href="format_converter.php" class="btn btn-outline-success">
                                        <i class="fas fa-exchange-alt"></i> Format Converter
                                        <br><small>Convert between formats</small>
                                    </a>
                                </div>
                            </div>
                            <div class="col-md-4">
                                <div class="d-grid">
                                    <a href="ai_assistant.php" class="btn btn-outline-warning">
                                        <i class="fas fa-robot"></i> AI Assistant
                                        <br><small>Get chemistry help</small>
                                    </a>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </main>

    <?php include '../includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script src="../assets/js/main.js"></script>
    
    <script>
    function copyToClipboard(text) {
        if (navigator.clipboard && window.isSecureContext) {
            navigator.clipboard.writeText(text).then(function() {
                showAlert('SMILES copied to clipboard!', 'success');
            }).catch(function(err) {
                showAlert('Failed to copy to clipboard', 'danger');
            });
        } else {
            // Fallback for older browsers
            const textArea = document.createElement('textarea');
            textArea.value = text;
            document.body.appendChild(textArea);
            textArea.select();
            try {
                document.execCommand('copy');
                showAlert('SMILES copied to clipboard!', 'success');
            } catch (err) {
                showAlert('Failed to copy to clipboard', 'danger');
            }
            document.body.removeChild(textArea);
        }
    }

    function showAlert(message, type) {
        const alertDiv = document.createElement('div');
        alertDiv.className = `alert alert-${type} alert-dismissible fade show position-fixed`;
        alertDiv.style.cssText = 'top: 20px; right: 20px; z-index: 9999; max-width: 300px;';
        alertDiv.innerHTML = `
            ${message}
            <button type="button" class="btn-close" data-bs-dismiss="alert"></button>
        `;
        document.body.appendChild(alertDiv);
        
        setTimeout(() => {
            alertDiv.remove();
        }, 3000);
    }

    // Add loading state to search button
    document.querySelector('form').addEventListener('submit', function() {
        const button = this.querySelector('button[type="submit"]');
        const originalText = button.innerHTML;
        button.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Searching...';
        button.disabled = true;
        
        setTimeout(() => {
            button.innerHTML = originalText;
            button.disabled = false;
        }, 10000); // Reset after 10 seconds if no response
    });

    // Add example click handlers
    document.addEventListener('DOMContentLoaded', function() {
        const examples = document.querySelectorAll('.badge');
        examples.forEach(function(badge) {
            badge.style.cursor = 'pointer';
            badge.addEventListener('click', function() {
                const searchInput = document.getElementById('search_query');
                const searchType = document.getElementById('search_type');
                
                const text = this.textContent.trim();
                searchInput.value = text;
                
                // Auto-detect search type
                if (text.match(/^[A-Z][a-z]*(\d+[A-Z][a-z]*\d*)*$/)) {
                    searchType.value = 'formula';
                } else if (text.match(/^\d{2,7}-\d{2}-\d$/)) {
                    searchType.value = 'cas';
                } else if (text.match(/^[A-Z\[\]@\+\-\(\)=\\\/\#\d]+$/i)) {
                    searchType.value = 'smiles';
                } else {
                    searchType.value = 'name';
                }
                
                showAlert('Example filled in! Click Search to continue.', 'info');
            });
        });
    });
    </script>
</body>
</html>
