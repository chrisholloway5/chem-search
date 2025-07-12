<?php
session_start();
require_once '../includes/config.php';
require_once '../includes/chemical_model.php';

// Check if database is enabled
if (!DATABASE_ENABLED) {
    die('<div class="alert alert-danger">Database is not enabled or not available.</div>');
}

$chemical = new Chemical();
$db = Database::getInstance();

// Handle actions
$action = $_GET['action'] ?? 'dashboard';
$message = '';

switch ($action) {
    case 'stats':
        $stats = getDatabaseStats();
        break;
    case 'recent':
        $recentAnalyses = $chemical->getRecentAnalyses(20);
        break;
    case 'popular':
        $popularSearches = $chemical->getPopularSearches(15, 30);
        break;
    case 'export':
        handleExport();
        break;
}

function getDatabaseStats() {
    global $db;
    
    $stats = [];
    
    // Chemical counts
    $result = $db->fetchOne("SELECT COUNT(*) as count FROM chemicals");
    $stats['total_chemicals'] = $result['count'];
    
    // Analysis counts
    $result = $db->fetchOne("SELECT COUNT(*) as count FROM analysis_cache");
    $stats['total_analyses'] = $result['count'];
    
    // Search counts
    $result = $db->fetchOne("SELECT COUNT(*) as count FROM search_history WHERE search_date > DATE_SUB(NOW(), INTERVAL 30 DAY)");
    $stats['searches_last_30_days'] = $result['count'];
    
    // Properties distribution
    $result = $db->fetchOne("
        SELECT 
            AVG(molecular_weight) as avg_mw,
            MIN(molecular_weight) as min_mw,
            MAX(molecular_weight) as max_mw,
            AVG(xlogp3) as avg_logp,
            COUNT(*) as count_with_props
        FROM molecular_properties 
        WHERE molecular_weight IS NOT NULL
    ");
    $stats['properties'] = $result;
    
    return $stats;
}

function handleExport() {
    global $chemical;
    
    if (!rateLimitCheck('export', 2, 3600)) {
        die(json_encode(['error' => 'Export rate limit exceeded. Please wait before exporting again.']));
    }
    
    $format = $_GET['format'] ?? 'json';
    $limit = min((int)($_GET['limit'] ?? 100), 1000); // Max 1000 records
    
    $sql = "
        SELECT 
            c.smiles, c.smiles_canonical, c.iupac_name, c.common_name, 
            c.molecular_formula, c.cas_number, c.inchi_key,
            mp.molecular_weight, mp.xlogp3, mp.hbd_count, mp.hba_count,
            mp.rotatable_bond_count, mp.tpsa, mp.heavy_atom_count,
            dp.lipinski_violations, dp.drug_like
        FROM chemicals c
        LEFT JOIN molecular_properties mp ON c.id = mp.chemical_id
        LEFT JOIN drug_properties dp ON c.id = dp.chemical_id
        ORDER BY c.created_at DESC
        LIMIT :limit
    ";
    
    $db = Database::getInstance();
    $data = $db->fetchAll($sql, ['limit' => $limit]);
    
    switch ($format) {
        case 'csv':
            header('Content-Type: text/csv');
            header('Content-Disposition: attachment; filename="chemicals_export.csv"');
            
            if (!empty($data)) {
                $output = fopen('php://output', 'w');
                fputcsv($output, array_keys($data[0]));
                foreach ($data as $row) {
                    fputcsv($output, $row);
                }
                fclose($output);
            }
            exit;
            
        case 'json':
        default:
            header('Content-Type: application/json');
            header('Content-Disposition: attachment; filename="chemicals_export.json"');
            echo json_encode($data, JSON_PRETTY_PRINT);
            exit;
    }
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Database Management - Chem Search</title>
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
                    <i class="fas fa-database text-primary"></i> Database Management
                </h1>
                
                <!-- Navigation Tabs -->
                <ul class="nav nav-tabs mb-4" id="dbTabs" role="tablist">
                    <li class="nav-item" role="presentation">
                        <button class="nav-link <?php echo $action === 'dashboard' ? 'active' : ''; ?>" 
                                onclick="location.href='?action=dashboard'">
                            <i class="fas fa-chart-bar"></i> Dashboard
                        </button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link <?php echo $action === 'recent' ? 'active' : ''; ?>" 
                                onclick="location.href='?action=recent'">
                            <i class="fas fa-clock"></i> Recent Analyses
                        </button>
                    </li>
                    <li class="nav-item" role="presentation">
                        <button class="nav-link <?php echo $action === 'popular' ? 'active' : ''; ?>" 
                                onclick="location.href='?action=popular'">
                            <i class="fas fa-fire"></i> Popular Searches
                        </button>
                    </li>
                </ul>
                
                <?php if ($action === 'dashboard' || $action === 'stats'): 
                    $stats = getDatabaseStats();
                ?>
                <!-- Dashboard -->
                <div class="row mb-4">
                    <div class="col-md-3">
                        <div class="card bg-primary text-white">
                            <div class="card-body">
                                <div class="d-flex justify-content-between">
                                    <div>
                                        <h4 class="mb-0"><?php echo number_format($stats['total_chemicals']); ?></h4>
                                        <p class="mb-0">Total Chemicals</p>
                                    </div>
                                    <i class="fas fa-flask fa-2x"></i>
                                </div>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card bg-success text-white">
                            <div class="card-body">
                                <div class="d-flex justify-content-between">
                                    <div>
                                        <h4 class="mb-0"><?php echo number_format($stats['total_analyses']); ?></h4>
                                        <p class="mb-0">Total Analyses</p>
                                    </div>
                                    <i class="fas fa-chart-line fa-2x"></i>
                                </div>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card bg-info text-white">
                            <div class="card-body">
                                <div class="d-flex justify-content-between">
                                    <div>
                                        <h4 class="mb-0"><?php echo number_format($stats['searches_last_30_days']); ?></h4>
                                        <p class="mb-0">Searches (30d)</p>
                                    </div>
                                    <i class="fas fa-search fa-2x"></i>
                                </div>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-3">
                        <div class="card bg-warning text-dark">
                            <div class="card-body">
                                <div class="d-flex justify-content-between">
                                    <div>
                                        <h4 class="mb-0"><?php echo number_format($stats['properties']['count_with_props']); ?></h4>
                                        <p class="mb-0">With Properties</p>
                                    </div>
                                    <i class="fas fa-cogs fa-2x"></i>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
                
                <!-- Property Statistics -->
                <?php if ($stats['properties']['count_with_props'] > 0): ?>
                <div class="row mb-4">
                    <div class="col-lg-6">
                        <div class="card">
                            <div class="card-header">
                                <h5><i class="fas fa-weight"></i> Molecular Weight Distribution</h5>
                            </div>
                            <div class="card-body">
                                <p><strong>Average:</strong> <?php echo round($stats['properties']['avg_mw'], 2); ?> Da</p>
                                <p><strong>Range:</strong> <?php echo round($stats['properties']['min_mw'], 2); ?> - <?php echo round($stats['properties']['max_mw'], 2); ?> Da</p>
                            </div>
                        </div>
                    </div>
                    <div class="col-lg-6">
                        <div class="card">
                            <div class="card-header">
                                <h5><i class="fas fa-tint"></i> Lipophilicity (LogP)</h5>
                            </div>
                            <div class="card-body">
                                <p><strong>Average LogP:</strong> <?php echo round($stats['properties']['avg_logp'], 2); ?></p>
                                <p><em>Optimal drug-like range: -0.4 to 5.6</em></p>
                            </div>
                        </div>
                    </div>
                </div>
                <?php endif; ?>
                
                <!-- Export Options -->
                <div class="card mb-4">
                    <div class="card-header">
                        <h5><i class="fas fa-download"></i> Export Data</h5>
                    </div>
                    <div class="card-body">
                        <p>Export chemical data in various formats for analysis or backup.</p>
                        <div class="btn-group" role="group">
                            <a href="?action=export&format=json&limit=100" class="btn btn-outline-primary">
                                <i class="fas fa-file-code"></i> JSON (100 records)
                            </a>
                            <a href="?action=export&format=csv&limit=100" class="btn btn-outline-success">
                                <i class="fas fa-file-csv"></i> CSV (100 records)
                            </a>
                            <a href="?action=export&format=json&limit=1000" class="btn btn-outline-warning">
                                <i class="fas fa-file-code"></i> JSON (1000 records)
                            </a>
                        </div>
                    </div>
                </div>
                
                <?php elseif ($action === 'recent'): ?>
                <!-- Recent Analyses -->
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-clock"></i> Recent Chemical Analyses</h5>
                    </div>
                    <div class="card-body">
                        <?php if (!empty($recentAnalyses)): ?>
                        <div class="table-responsive">
                            <table class="table table-striped">
                                <thead>
                                    <tr>
                                        <th>SMILES</th>
                                        <th>Name</th>
                                        <th>Formula</th>
                                        <th>MW</th>
                                        <th>Analysis Date</th>
                                        <th>Actions</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <?php foreach ($recentAnalyses as $analysis): ?>
                                    <tr>
                                        <td>
                                            <code class="small"><?php echo htmlspecialchars(substr($analysis['smiles'], 0, 20)); ?><?php echo strlen($analysis['smiles']) > 20 ? '...' : ''; ?></code>
                                        </td>
                                        <td><?php echo htmlspecialchars($analysis['common_name'] ?: $analysis['iupac_name'] ?: 'Unknown'); ?></td>
                                        <td><?php echo htmlspecialchars($analysis['molecular_formula'] ?: '-'); ?></td>
                                        <td><?php echo $analysis['molecular_weight'] ? round($analysis['molecular_weight'], 1) . ' Da' : '-'; ?></td>
                                        <td><?php echo date('M j, Y H:i', strtotime($analysis['search_date'])); ?></td>
                                        <td>
                                            <a href="../tools/molecule_analyzer.php?smiles=<?php echo urlencode($analysis['smiles']); ?>" 
                                               class="btn btn-sm btn-outline-primary">
                                                <i class="fas fa-eye"></i> View
                                            </a>
                                        </td>
                                    </tr>
                                    <?php endforeach; ?>
                                </tbody>
                            </table>
                        </div>
                        <?php else: ?>
                        <p class="text-muted">No recent analyses found.</p>
                        <?php endif; ?>
                    </div>
                </div>
                
                <?php elseif ($action === 'popular'): ?>
                <!-- Popular Searches -->
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-fire"></i> Popular Searches (Last 30 Days)</h5>
                    </div>
                    <div class="card-body">
                        <?php if (!empty($popularSearches)): ?>
                        <div class="table-responsive">
                            <table class="table table-striped">
                                <thead>
                                    <tr>
                                        <th>Search Query</th>
                                        <th>Type</th>
                                        <th>Search Count</th>
                                        <th>Last Search</th>
                                        <th>Actions</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    <?php foreach ($popularSearches as $search): ?>
                                    <tr>
                                        <td>
                                            <code><?php echo htmlspecialchars(substr($search['search_query'], 0, 30)); ?><?php echo strlen($search['search_query']) > 30 ? '...' : ''; ?></code>
                                        </td>
                                        <td>
                                            <span class="badge bg-secondary"><?php echo htmlspecialchars($search['search_type']); ?></span>
                                        </td>
                                        <td>
                                            <span class="badge bg-primary"><?php echo $search['search_count']; ?></span>
                                        </td>
                                        <td><?php echo date('M j, Y H:i', strtotime($search['last_search'])); ?></td>
                                        <td>
                                            <?php if ($search['search_type'] === 'smiles'): ?>
                                            <a href="../tools/molecule_analyzer.php?smiles=<?php echo urlencode($search['search_query']); ?>" 
                                               class="btn btn-sm btn-outline-primary">
                                                <i class="fas fa-search"></i> Analyze
                                            </a>
                                            <?php endif; ?>
                                        </td>
                                    </tr>
                                    <?php endforeach; ?>
                                </tbody>
                            </table>
                        </div>
                        <?php else: ?>
                        <p class="text-muted">No popular searches found.</p>
                        <?php endif; ?>
                    </div>
                </div>
                <?php endif; ?>
            </div>
        </div>
    </main>
    
    <?php include '../includes/footer.php'; ?>
    
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
