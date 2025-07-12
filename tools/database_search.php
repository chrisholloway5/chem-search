<?php
session_start();
require_once '../includes/config.php';
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
        <div class="row">
            <div class="col-lg-12">
                <h1 class="mb-4">
                    <i class="fas fa-search text-info"></i> Chemical Database Search
                </h1>
                <p class="lead">Search chemical databases like PubChem for compound information.</p>
            </div>
        </div>

        <!-- Coming Soon Notice -->
        <div class="row">
            <div class="col-lg-8 mx-auto">
                <div class="card border-info shadow">
                    <div class="card-header bg-info text-white text-center">
                        <h5><i class="fas fa-construction"></i> Feature Coming Soon</h5>
                    </div>
                    <div class="card-body text-center">
                        <i class="fas fa-database fa-5x text-info mb-4"></i>
                        <h4>Database Search Tool</h4>
                        <p class="lead">This feature will allow you to search chemical databases including:</p>
                        
                        <div class="row mt-4">
                            <div class="col-md-6">
                                <h6><i class="fas fa-check text-success"></i> Planned Features:</h6>
                                <ul class="list-unstyled text-start">
                                    <li>• PubChem compound search</li>
                                    <li>• ChEMBL database integration</li>
                                    <li>• Chemical property lookup</li>
                                    <li>• Biological activity data</li>
                                    <li>• Safety and toxicity information</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-search text-info"></i> Search Methods:</h6>
                                <ul class="list-unstyled text-start">
                                    <li>• By compound name</li>
                                    <li>• By SMILES notation</li>
                                    <li>• By InChI/InChI Key</li>
                                    <li>• By CAS number</li>
                                    <li>• By molecular formula</li>
                                </ul>
                            </div>
                        </div>
                        
                        <hr>
                        <p class="text-muted">In the meantime, you can use the other available tools:</p>
                        <div class="d-grid gap-2 d-md-flex justify-content-center">
                            <a href="molecule_analyzer.php" class="btn btn-primary">
                                <i class="fas fa-atom"></i> Molecule Analyzer
                            </a>
                            <a href="format_converter.php" class="btn btn-success">
                                <i class="fas fa-exchange-alt"></i> Format Converter
                            </a>
                            <a href="ai_assistant.php" class="btn btn-warning">
                                <i class="fas fa-robot"></i> AI Assistant
                            </a>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </main>

    <?php include '../includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script src="../assets/js/main.js"></script>
</body>
</html>
