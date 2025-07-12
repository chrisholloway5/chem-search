<?php
session_start();
require_once '../includes/config.php';

// Initialize variables
$query = '';
$searchResults = null;
$analysisResults = null;
$errorMessage = '';
$selectedSources = ['arxiv', 'pubmed', 'crossref'];
$maxResults = 30;
$searchPerformed = false;

// Handle form submission
if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    $query = trim($_POST['query'] ?? '');
    $selectedSources = $_POST['sources'] ?? ['arxiv'];
    $maxResults = min(max((int)($_POST['max_results'] ?? 30), 5), 100); // Limit between 5 and 100
    $includeAnalysis = isset($_POST['include_analysis']);
    
    if (!empty($query)) {
        try {
            $searchPerformed = true;
            
            // Prepare command for Python script
            $pythonPath = getPythonPath();
            $scriptPath = dirname(__FILE__) . '/../python/paper_search.py';
            
            // Escape arguments for command line
            $escapedQuery = escapeshellarg($query);
            $sourcesStr = implode(',', $selectedSources);
            $escapedSources = escapeshellarg($sourcesStr);
            $escapedMaxResults = escapeshellarg($maxResults);
            
            $command = "$pythonPath " . escapeshellarg($scriptPath) . " $escapedQuery $escapedSources $escapedMaxResults 2>&1";
            
            // Execute search
            $output = shell_exec($command);
            
            if ($output) {
                $searchResults = json_decode($output, true);
                
                if (json_last_error() === JSON_ERROR_NONE && $searchResults) {
                    // If analysis was requested and we have papers
                    if ($includeAnalysis && !empty($searchResults['papers'])) {
                        // Create temporary file with results for analysis
                        $tempFile = tempnam(sys_get_temp_dir(), 'papers_');
                        file_put_contents($tempFile, json_encode($searchResults));
                        
                        // Run analysis
                        $analyzerPath = dirname(__FILE__) . '/../python/paper_analyzer.py';
                        $analysisCommand = "$pythonPath " . escapeshellarg($analyzerPath) . " " . escapeshellarg($tempFile) . " $escapedQuery 2>&1";
                        
                        $analysisOutput = shell_exec($analysisCommand);
                        if ($analysisOutput) {
                            $analysisResults = json_decode($analysisOutput, true);
                        }
                        
                        // Clean up temp file
                        unlink($tempFile);
                    }
                } else {
                    $errorMessage = 'Invalid response from search service';
                }
            } else {
                $errorMessage = 'No response from search service';
            }
            
        } catch (Exception $e) {
            $errorMessage = 'Error performing search: ' . $e->getMessage();
        }
    } else {
        $errorMessage = 'Please enter a search query';
    }
}

// Helper function to get Python executable path
function getPythonPath() {
    // Try different common Python paths
    $paths = ['python', 'python3', 'py'];
    foreach ($paths as $path) {
        $test = shell_exec("$path --version 2>&1");
        if ($test && strpos($test, 'Python') !== false) {
            return $path;
        }
    }
    return 'python'; // fallback
}

function formatDate($dateString) {
    if (!$dateString) return 'Unknown';
    try {
        $date = new DateTime($dateString);
        return $date->format('M j, Y');
    } catch (Exception $e) {
        return $dateString;
    }
}

function truncateText($text, $maxLength = 300) {
    if (strlen($text) <= $maxLength) return $text;
    return substr($text, 0, $maxLength) . '...';
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Research Paper Search - Chem Search</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
    <link href="../assets/css/style.css" rel="stylesheet">
    <style>
        .search-form {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 15px;
            padding: 2rem;
            margin-bottom: 2rem;
        }
        
        .paper-card {
            transition: transform 0.2s ease-in-out;
            border-left: 4px solid #007bff;
        }
        
        .paper-card:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }
        
        .source-badge {
            font-size: 0.75rem;
            padding: 0.25rem 0.5rem;
        }
        
        .analysis-section {
            background: #f8f9fa;
            border-radius: 10px;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
        }
        
        .trend-item {
            background: white;
            border-radius: 8px;
            padding: 1rem;
            margin-bottom: 0.5rem;
            border-left: 3px solid #28a745;
        }
        
        .search-stats {
            background: linear-gradient(45deg, #667eea, #764ba2);
            color: white;
            border-radius: 10px;
            padding: 1rem;
        }
        
        .example-queries {
            background: #ffffff;
            border: 1px solid #e3e6f0;
            border-radius: 10px;
            padding: 1.5rem;
        }
        
        .query-tag {
            display: inline-block;
            background: #e3f2fd;
            color: #1976d2;
            padding: 0.25rem 0.75rem;
            border-radius: 15px;
            margin: 0.25rem;
            cursor: pointer;
            transition: all 0.2s;
        }
        
        .query-tag:hover {
            background: #1976d2;
            color: white;
        }
        
        .loading-overlay {
            position: fixed;
            top: 0;
            left: 0;
            width: 100%;
            height: 100%;
            background: rgba(0,0,0,0.7);
            display: flex;
            align-items: center;
            justify-content: center;
            z-index: 9999;
        }
        
        .citation-potential {
            border-radius: 50%;
            width: 40px;
            height: 40px;
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-weight: bold;
            font-size: 0.85rem;
        }
        
        .potential-high { background: #28a745; }
        .potential-medium { background: #ffc107; color: #000; }
        .potential-standard { background: #6c757d; }
    </style>
</head>
<body>
    <?php include '../includes/header.php'; ?>
    
    <main class="container my-5">
        <!-- Header -->
        <div class="row mb-4">
            <div class="col-lg-12">
                <h1 class="mb-3">
                    <i class="fas fa-file-alt text-primary"></i> Research Paper Search
                </h1>
                <p class="lead">AI-powered search across multiple academic databases for chemistry and materials science research.</p>
            </div>
        </div>

        <!-- Search Form -->
        <div class="row">
            <div class="col-lg-12">
                <div class="search-form">
                    <form method="POST" id="searchForm">
                        <div class="row align-items-end">
                            <div class="col-md-6">
                                <label for="query" class="form-label">
                                    <i class="fas fa-search"></i> Research Query
                                </label>
                                <input type="text" class="form-control form-control-lg" id="query" name="query" 
                                       value="<?php echo htmlspecialchars($query); ?>" 
                                       placeholder="Enter your research question or keywords..." required>
                            </div>
                            <div class="col-md-3">
                                <label for="max_results" class="form-label">
                                    <i class="fas fa-list-ol"></i> Max Results
                                </label>
                                <select class="form-select" id="max_results" name="max_results">
                                    <option value="20" <?php echo $maxResults == 20 ? 'selected' : ''; ?>>20 papers</option>
                                    <option value="30" <?php echo $maxResults == 30 ? 'selected' : ''; ?>>30 papers</option>
                                    <option value="50" <?php echo $maxResults == 50 ? 'selected' : ''; ?>>50 papers</option>
                                    <option value="100" <?php echo $maxResults == 100 ? 'selected' : ''; ?>>100 papers</option>
                                </select>
                            </div>
                            <div class="col-md-3">
                                <button type="submit" class="btn btn-light btn-lg w-100">
                                    <i class="fas fa-search"></i> Search Papers
                                </button>
                            </div>
                        </div>
                        
                        <div class="row mt-3">
                            <div class="col-md-6">
                                <label class="form-label">
                                    <i class="fas fa-database"></i> Data Sources
                                </label>
                                <div class="d-flex gap-3">
                                    <div class="form-check">
                                        <input class="form-check-input" type="checkbox" name="sources[]" value="arxiv" id="arxiv"
                                               <?php echo in_array('arxiv', $selectedSources) ? 'checked' : ''; ?>>
                                        <label class="form-check-label" for="arxiv">ArXiv</label>
                                    </div>
                                    <div class="form-check">
                                        <input class="form-check-input" type="checkbox" name="sources[]" value="pubmed" id="pubmed"
                                               <?php echo in_array('pubmed', $selectedSources) ? 'checked' : ''; ?>>
                                        <label class="form-check-label" for="pubmed">PubMed</label>
                                    </div>
                                    <div class="form-check">
                                        <input class="form-check-input" type="checkbox" name="sources[]" value="crossref" id="crossref"
                                               <?php echo in_array('crossref', $selectedSources) ? 'checked' : ''; ?>>
                                        <label class="form-check-label" for="crossref">CrossRef</label>
                                    </div>
                                </div>
                            </div>
                            <div class="col-md-6">
                                <label class="form-label">
                                    <i class="fas fa-brain"></i> AI Analysis
                                </label>
                                <div class="form-check">
                                    <input class="form-check-input" type="checkbox" name="include_analysis" id="include_analysis"
                                           <?php echo isset($_POST['include_analysis']) ? 'checked' : ''; ?>>
                                    <label class="form-check-label" for="include_analysis">
                                        Include comprehensive AI analysis (trends, gaps, insights)
                                    </label>
                                </div>
                            </div>
                        </div>
                    </form>
                </div>
            </div>
        </div>

        <!-- Example Queries -->
        <?php if (!$searchPerformed): ?>
        <div class="row">
            <div class="col-lg-12">
                <div class="example-queries">
                    <h5><i class="fas fa-lightbulb text-warning"></i> Example Research Queries</h5>
                    <p class="text-muted mb-3">Click on any example to search, or create your own query</p>
                    
                    <div class="row">
                        <div class="col-md-6">
                            <h6 class="text-primary">Synthesis & Materials</h6>
                            <div class="query-tag" onclick="setQuery('perovskite solar cells synthesis methods')">
                                Perovskite solar cell synthesis
                            </div>
                            <div class="query-tag" onclick="setQuery('metal-organic frameworks MOF synthesis')">
                                MOF synthesis methods
                            </div>
                            <div class="query-tag" onclick="setQuery('graphene oxide synthesis applications')">
                                Graphene oxide synthesis
                            </div>
                            <div class="query-tag" onclick="setQuery('nanoparticle synthesis green chemistry')">
                                Green nanoparticle synthesis
                            </div>
                        </div>
                        <div class="col-md-6">
                            <h6 class="text-success">Drug Discovery & Catalysis</h6>
                            <div class="query-tag" onclick="setQuery('COVID-19 drug targets therapeutic compounds')">
                                COVID-19 drug targets
                            </div>
                            <div class="query-tag" onclick="setQuery('heterogeneous catalysts CO2 reduction')">
                                CO2 reduction catalysts
                            </div>
                            <div class="query-tag" onclick="setQuery('machine learning drug discovery QSAR')">
                                ML in drug discovery
                            </div>
                            <div class="query-tag" onclick="setQuery('electrocatalysis water splitting hydrogen')">
                                Electrocatalytic water splitting
                            </div>
                        </div>
                    </div>
                    
                    <div class="row mt-3">
                        <div class="col-md-6">
                            <h6 class="text-info">Energy & Environment</h6>
                            <div class="query-tag" onclick="setQuery('lithium ion battery electrode materials')">
                                Battery electrode materials
                            </div>
                            <div class="query-tag" onclick="setQuery('photocatalysis environmental remediation')">
                                Photocatalytic remediation
                            </div>
                            <div class="query-tag" onclick="setQuery('fuel cells polymer electrolyte membrane')">
                                Fuel cell membranes
                            </div>
                        </div>
                        <div class="col-md-6">
                            <h6 class="text-danger">Advanced Materials</h6>
                            <div class="query-tag" onclick="setQuery('2D materials properties applications')">
                                2D materials applications
                            </div>
                            <div class="query-tag" onclick="setQuery('quantum dots synthesis optical properties')">
                                Quantum dots synthesis
                            </div>
                            <div class="query-tag" onclick="setQuery('biomaterials tissue engineering scaffold')">
                                Biomaterial scaffolds
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        <?php endif; ?>

        <!-- Error Message -->
        <?php if ($errorMessage): ?>
        <div class="row">
            <div class="col-lg-12">
                <div class="alert alert-danger">
                    <i class="fas fa-exclamation-triangle"></i>
                    <strong>Error:</strong> <?php echo htmlspecialchars($errorMessage); ?>
                </div>
            </div>
        </div>
        <?php endif; ?>

        <!-- Search Results -->
        <?php if ($searchResults && isset($searchResults['success']) && $searchResults['success']): ?>
        
        <!-- Search Statistics -->
        <div class="row">
            <div class="col-lg-12">
                <div class="search-stats">
                    <div class="row align-items-center">
                        <div class="col-md-8">
                            <h5 class="mb-1">
                                <i class="fas fa-chart-bar"></i> 
                                Found <?php echo $searchResults['total_found']; ?> papers for: 
                                "<?php echo htmlspecialchars($searchResults['query']); ?>"
                            </h5>
                            <p class="mb-0">
                                Sources: <?php echo implode(', ', array_map('ucfirst', $searchResults['sources_searched'])); ?>
                                <?php if (isset($searchResults['source_breakdown'])): ?>
                                    <br>
                                    <small>
                                        <?php foreach ($searchResults['source_breakdown'] as $source => $count): ?>
                                            <?php echo ucfirst($source); ?>: <?php echo $count; ?> papers<?php echo next($searchResults['source_breakdown']) !== false ? ', ' : ''; ?>
                                        <?php endforeach; ?>
                                    </small>
                                <?php endif; ?>
                            </p>
                        </div>
                        <div class="col-md-4 text-end">
                            <button class="btn btn-light" onclick="exportResults()">
                                <i class="fas fa-download"></i> Export Results
                            </button>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- AI Analysis Results -->
        <?php if ($analysisResults && !isset($analysisResults['error'])): ?>
        <div class="row mt-4">
            <div class="col-lg-12">
                <div class="analysis-section">
                    <h4><i class="fas fa-brain text-primary"></i> AI Analysis & Insights</h4>
                    
                    <!-- AI Summary -->
                    <?php if (isset($analysisResults['ai_summary']) && $analysisResults['ai_summary']): ?>
                    <div class="card mb-3">
                        <div class="card-header bg-primary text-white">
                            <i class="fas fa-robot"></i> AI-Generated Summary
                        </div>
                        <div class="card-body">
                            <p><?php echo nl2br(htmlspecialchars($analysisResults['ai_summary'])); ?></p>
                        </div>
                    </div>
                    <?php endif; ?>
                    
                    <div class="row">
                        <!-- Key Findings -->
                        <div class="col-lg-6">
                            <?php if (isset($analysisResults['key_findings']) && !empty($analysisResults['key_findings'])): ?>
                            <div class="card mb-3">
                                <div class="card-header">
                                    <i class="fas fa-lightbulb text-warning"></i> Key Findings
                                </div>
                                <div class="card-body">
                                    <?php foreach (array_slice($analysisResults['key_findings'], 0, 3) as $finding): ?>
                                    <div class="trend-item">
                                        <h6 class="mb-1"><?php echo htmlspecialchars($finding['paper_title']); ?></h6>
                                        <small class="text-muted">Score: <?php echo number_format($finding['score'], 1); ?></small>
                                        <ul class="mt-2 mb-0">
                                            <?php foreach ($finding['key_sentences'] as $sentence): ?>
                                            <li><small><?php echo htmlspecialchars($sentence); ?></small></li>
                                            <?php endforeach; ?>
                                        </ul>
                                    </div>
                                    <?php endforeach; ?>
                                </div>
                            </div>
                            <?php endif; ?>
                        </div>
                        
                        <!-- Research Areas -->
                        <div class="col-lg-6">
                            <?php if (isset($analysisResults['research_areas']) && !empty($analysisResults['research_areas'])): ?>
                            <div class="card mb-3">
                                <div class="card-header">
                                    <i class="fas fa-tags text-info"></i> Research Areas
                                </div>
                                <div class="card-body">
                                    <?php foreach ($analysisResults['research_areas'] as $area => $papers): ?>
                                    <div class="mb-2">
                                        <span class="badge bg-info"><?php echo ucfirst(str_replace('_', ' ', $area)); ?></span>
                                        <small class="text-muted"><?php echo count($papers); ?> papers</small>
                                    </div>
                                    <?php endforeach; ?>
                                </div>
                            </div>
                            <?php endif; ?>
                        </div>
                    </div>
                    
                    <!-- Research Gaps -->
                    <?php if (isset($analysisResults['research_gaps']) && !empty($analysisResults['research_gaps'])): ?>
                    <div class="card mb-3">
                        <div class="card-header bg-warning text-dark">
                            <i class="fas fa-exclamation-circle"></i> Identified Research Gaps
                        </div>
                        <div class="card-body">
                            <ul class="mb-0">
                                <?php foreach ($analysisResults['research_gaps'] as $gap): ?>
                                <li><?php echo htmlspecialchars($gap); ?></li>
                                <?php endforeach; ?>
                            </ul>
                        </div>
                    </div>
                    <?php endif; ?>
                    
                    <!-- Hot Topics -->
                    <?php if (isset($analysisResults['trends_analysis']['hot_topics']) && !empty($analysisResults['trends_analysis']['hot_topics'])): ?>
                    <div class="card">
                        <div class="card-header">
                            <i class="fas fa-fire text-danger"></i> Hot Topics
                        </div>
                        <div class="card-body">
                            <?php foreach (array_slice($analysisResults['trends_analysis']['hot_topics'], 0, 10) as $topic): ?>
                            <span class="badge bg-secondary me-1 mb-1"><?php echo htmlspecialchars($topic[0]); ?> (<?php echo $topic[1]; ?>)</span>
                            <?php endforeach; ?>
                        </div>
                    </div>
                    <?php endif; ?>
                </div>
            </div>
        </div>
        <?php endif; ?>

        <!-- Papers List -->
        <div class="row mt-4">
            <div class="col-lg-12">
                <h4><i class="fas fa-list"></i> Research Papers</h4>
                
                <?php foreach ($searchResults['papers'] as $index => $paper): ?>
                <div class="card paper-card mb-3">
                    <div class="card-header d-flex justify-content-between align-items-center">
                        <div class="d-flex align-items-center gap-2">
                            <span class="badge bg-primary">#<?php echo $index + 1; ?></span>
                            <span class="source-badge badge bg-<?php 
                                echo $paper['source'] === 'arxiv' ? 'info' : 
                                    ($paper['source'] === 'pubmed' ? 'success' : 'secondary'); 
                            ?>">
                                <?php echo strtoupper($paper['source']); ?>
                            </span>
                            
                            <!-- Citation Potential -->
                            <?php if ($analysisResults && isset($analysisResults['citation_analysis'])): ?>
                                <?php foreach ($analysisResults['citation_analysis'] as $citation): ?>
                                    <?php if ($citation['title'] === $paper['title']): ?>
                                        <div class="citation-potential <?php 
                                            echo $citation['recommendation'] === 'High potential' ? 'potential-high' : 
                                                ($citation['recommendation'] === 'Medium potential' ? 'potential-medium' : 'potential-standard'); 
                                        ?>" title="Citation Potential: <?php echo $citation['recommendation']; ?>">
                                            <?php echo $citation['citation_potential_score']; ?>
                                        </div>
                                        <?php break; ?>
                                    <?php endif; ?>
                                <?php endforeach; ?>
                            <?php endif; ?>
                        </div>
                        
                        <div class="text-end">
                            <?php if ($paper['score'] > 0): ?>
                            <small class="text-muted">Relevance: <?php echo number_format($paper['score'], 1); ?></small>
                            <?php endif; ?>
                        </div>
                    </div>
                    
                    <div class="card-body">
                        <h5 class="card-title">
                            <?php if ($paper['url']): ?>
                            <a href="<?php echo htmlspecialchars($paper['url']); ?>" target="_blank" class="text-decoration-none">
                                <?php echo htmlspecialchars($paper['title']); ?>
                                <i class="fas fa-external-link-alt fa-sm"></i>
                            </a>
                            <?php else: ?>
                            <?php echo htmlspecialchars($paper['title']); ?>
                            <?php endif; ?>
                        </h5>
                        
                        <!-- Authors -->
                        <p class="text-muted mb-2">
                            <i class="fas fa-users"></i>
                            <strong>Authors:</strong> 
                            <?php echo !empty($paper['authors']) ? htmlspecialchars(implode(', ', array_slice($paper['authors'], 0, 5))) : 'Unknown'; ?>
                            <?php if (count($paper['authors']) > 5): ?>
                                <em>et al.</em>
                            <?php endif; ?>
                        </p>
                        
                        <!-- Publication Info -->
                        <p class="text-muted mb-2">
                            <?php if ($paper['journal']): ?>
                            <i class="fas fa-journal-whills"></i>
                            <strong>Journal:</strong> <?php echo htmlspecialchars($paper['journal']); ?>
                            <?php endif; ?>
                            
                            <?php if ($paper['published_date']): ?>
                            <span class="ms-3">
                                <i class="fas fa-calendar"></i>
                                <strong>Published:</strong> <?php echo formatDate($paper['published_date']); ?>
                            </span>
                            <?php endif; ?>
                        </p>
                        
                        <!-- IDs -->
                        <p class="text-muted mb-2">
                            <?php if ($paper['doi']): ?>
                            <span class="badge bg-light text-dark">DOI: <?php echo htmlspecialchars($paper['doi']); ?></span>
                            <?php endif; ?>
                            <?php if ($paper['pmid']): ?>
                            <span class="badge bg-light text-dark">PMID: <?php echo htmlspecialchars($paper['pmid']); ?></span>
                            <?php endif; ?>
                            <?php if ($paper['arxiv_id']): ?>
                            <span class="badge bg-light text-dark">arXiv: <?php echo htmlspecialchars($paper['arxiv_id']); ?></span>
                            <?php endif; ?>
                        </p>
                        
                        <!-- Abstract -->
                        <p class="card-text">
                            <strong>Abstract:</strong> 
                            <?php echo nl2br(htmlspecialchars(truncateText($paper['abstract'], 400))); ?>
                        </p>
                        
                        <!-- Keywords -->
                        <?php if (!empty($paper['keywords'])): ?>
                        <div class="mt-2">
                            <small class="text-muted"><strong>Keywords:</strong></small><br>
                            <?php foreach (array_slice($paper['keywords'], 0, 8) as $keyword): ?>
                            <span class="badge bg-light text-dark me-1"><?php echo htmlspecialchars($keyword); ?></span>
                            <?php endforeach; ?>
                        </div>
                        <?php endif; ?>
                    </div>
                </div>
                <?php endforeach; ?>
            </div>
        </div>
        
        <?php elseif ($searchPerformed): ?>
        <!-- No Results -->
        <div class="row">
            <div class="col-lg-12">
                <div class="alert alert-info">
                    <i class="fas fa-info-circle"></i>
                    <strong>No papers found</strong> for your query. Try different keywords or broader search terms.
                </div>
            </div>
        </div>
        <?php endif; ?>

        <!-- ChemCrow Reference -->
        <div class="row mt-5">
            <div class="col-lg-12">
                <div class="card border-primary">
                    <div class="card-header bg-primary text-white">
                        <h5><i class="fas fa-graduation-cap"></i> ChemCrow Research Paper</h5>
                    </div>
                    <div class="card-body">
                        <div class="row align-items-center">
                            <div class="col-md-8">
                                <h6>ChemCrow: Augmenting large-language models with chemistry tools</h6>
                                <p class="text-muted mb-2">
                                    <strong>Authors:</strong> Andres M Bran, Sam Cox, Oliver Schilter, Carlo Baldassari, Andrew D White, Philippe Schwaller
                                </p>
                                <p class="text-muted mb-2">
                                    <strong>Abstract:</strong> This paper introduces ChemCrow, an LLM chemistry agent designed to accomplish tasks across organic synthesis, drug discovery, and materials design. ChemCrow integrates 18 expert-designed tools and uses GPT-4 as the LLM to augment it with chemical knowledge...
                                </p>
                            </div>
                            <div class="col-md-4 text-center">
                                <a href="https://arxiv.org/abs/2304.05376" target="_blank" class="btn btn-primary">
                                    <i class="fas fa-external-link-alt"></i> Read Paper
                                </a>
                                <br>
                                <small class="text-muted mt-2 d-block">arXiv:2304.05376</small>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </main>

    <?php include '../includes/footer.php'; ?>

    <!-- Loading Overlay -->
    <div id="loadingOverlay" class="loading-overlay" style="display: none;">
        <div class="text-center text-white">
            <div class="spinner-border spinner-border-lg mb-3" role="status">
                <span class="visually-hidden">Loading...</span>
            </div>
            <h5>Searching academic databases...</h5>
            <p>This may take a few moments</p>
        </div>
    </div>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script src="../assets/js/main.js"></script>
    
    <script>
        // Set query from example
        function setQuery(query) {
            document.getElementById('query').value = query;
            document.getElementById('query').focus();
        }
        
        // Show loading overlay on form submit
        document.getElementById('searchForm').addEventListener('submit', function() {
            document.getElementById('loadingOverlay').style.display = 'flex';
        });
        
        // Export results function
        function exportResults() {
            const results = <?php echo $searchResults ? json_encode($searchResults) : 'null'; ?>;
            if (results) {
                const dataStr = JSON.stringify(results, null, 2);
                const dataBlob = new Blob([dataStr], {type: 'application/json'});
                const url = URL.createObjectURL(dataBlob);
                const link = document.createElement('a');
                link.href = url;
                link.download = 'paper_search_results.json';
                link.click();
                URL.revokeObjectURL(url);
            }
        }
        
        // Auto-hide loading overlay if page loads with results
        window.addEventListener('load', function() {
            document.getElementById('loadingOverlay').style.display = 'none';
        });
    </script>
</body>
</html>
