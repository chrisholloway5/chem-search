<?php
session_start();
?>
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Chem Search - AI-Powered Chemistry Platform</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <link href="assets/css/style.css" rel="stylesheet">
</head>
<body>
    <?php include 'includes/header.php'; ?>
    
    <main>
        <!-- Hero Section -->
        <section class="hero-section py-5 bg-gradient">
            <div class="container">
                <div class="row align-items-center">
                    <div class="col-lg-6">
                        <h1 class="display-4 mb-4" style="color: black !important; text-shadow: 2px 2px 4px rgba(255,255,255,0.8);">
                            <i class="fas fa-flask" style="color: #ffc107 !important;"></i> Chem Search
                        </h1>
                        <p class="lead mb-4" style="color: #333 !important; text-shadow: 1px 1px 2px rgba(255,255,255,0.8);">
                            AI-powered chemistry platform for molecular analysis, property calculations, 
                            and chemical research. Built with Chem and RDKit.
                        </p>
                        <div class="d-grid gap-2 d-md-flex">
                            <a href="tools/molecule_analyzer.php" class="btn btn-warning btn-lg">
                                <i class="fas fa-atom"></i> Analyze Molecules
                            </a>
                            <a href="tools/ai_assistant.php" class="btn btn-outline-light btn-lg">
                                <i class="fas fa-robot"></i> AI Assistant
                            </a>
                        </div>
                    </div>
                    <div class="col-lg-6">
                        <div class="text-center">
                            <i class="fas fa-dna fa-10x" style="color: #ffc107 !important; opacity: 0.75;"></i>
                        </div>
                    </div>
                </div>
            </div>
        </section>

        <!-- Features Section -->
        <section class="py-5">
            <div class="container">
                <h2 class="text-center mb-5">Chemistry Tools & Features</h2>
                <div class="row g-4">
                    <div class="col-md-6 col-lg-4">
                        <div class="card h-100 shadow-sm">
                            <div class="card-body text-center">
                                <i class="fas fa-calculator fa-3x text-primary mb-3"></i>
                                <h5 class="card-title">Molecular Properties</h5>
                                <p class="card-text">Calculate molecular weight, LogP, polar surface area, and other important molecular descriptors.</p>
                                <a href="tools/molecule_analyzer.php" class="btn btn-primary">Try Now</a>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-4">
                        <div class="card h-100 shadow-sm">
                            <div class="card-body text-center">
                                <i class="fas fa-exchange-alt fa-3x text-success mb-3"></i>
                                <h5 class="card-title">Format Converter</h5>
                                <p class="card-text">Convert between SMILES, InChI, molecular formulas, and other chemical formats.</p>
                                <a href="tools/format_converter.php" class="btn btn-success">Convert</a>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-4">
                        <div class="card h-100 shadow-sm">
                            <div class="card-body text-center">
                                <i class="fas fa-robot fa-3x text-warning mb-3"></i>
                                <h5 class="card-title">AI Chemistry Assistant</h5>
                                <p class="card-text">Ask complex chemistry questions powered by Chem' AI agent system.</p>
                                <a href="tools/ai_assistant.php" class="btn btn-warning">Ask AI</a>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-4">
                        <div class="card h-100 shadow-sm">
                            <div class="card-body text-center">
                                <i class="fas fa-search fa-3x text-info mb-3"></i>
                                <h5 class="card-title">Database Search</h5>
                                <p class="card-text">Search chemical databases like PubChem for compound information.</p>
                                <a href="tools/database_search.php" class="btn btn-info">Search</a>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-4">
                        <div class="card h-100 shadow-sm">
                            <div class="card-body text-center">
                                <i class="fas fa-project-diagram fa-3x text-danger mb-3"></i>
                                <h5 class="card-title">Structure Viewer</h5>
                                <p class="card-text">Visualize molecular structures in 2D and 3D formats.</p>
                                <a href="tools/structure_viewer.php" class="btn btn-danger">Visualize</a>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-4">
                        <div class="card h-100 shadow-sm">
                            <div class="card-body text-center">
                                <i class="fas fa-file-alt fa-3x text-secondary mb-3"></i>
                                <h5 class="card-title">Research Papers</h5>
                                <p class="card-text">Find and analyze chemistry research papers using AI-powered search.</p>
                                <a href="tools/paper_search.php" class="btn btn-secondary">Research</a>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </section>

        <!-- Demo Section -->
        <section class="py-5 bg-light">
            <div class="container">
                <div class="row">
                    <div class="col-lg-8 mx-auto text-center">
                        <h2 class="mb-4">Quick Demo</h2>
                        <p class="lead mb-4">Try analyzing a common molecule like caffeine:</p>
                        <div class="card">
                            <div class="card-body">
                                <form action="tools/molecule_analyzer.php" method="GET" class="d-flex gap-2">
                                    <input type="text" name="smiles" class="form-control" 
                                           value="CN1C=NC2=C1C(=O)N(C(=O)N2C)C" 
                                           placeholder="Enter SMILES notation">
                                    <button type="submit" class="btn btn-primary">
                                        <i class="fas fa-search"></i> Analyze
                                    </button>
                                </form>
                                <small class="text-muted">Example: Caffeine SMILES notation</small>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </section>
    </main>

    <?php include 'includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script src="assets/js/main.js"></script>
</body>
</html>
