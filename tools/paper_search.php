<?php
session_start();
require_once '../includes/config.php';
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Paper Search - Chem Search</title>
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
                    <i class="fas fa-file-alt text-secondary"></i> Research Paper Search
                </h1>
                <p class="lead">Find and analyze chemistry research papers using AI-powered search.</p>
            </div>
        </div>

        <!-- Coming Soon Notice -->
        <div class="row">
            <div class="col-lg-8 mx-auto">
                <div class="card border-secondary shadow">
                    <div class="card-header bg-secondary text-white text-center">
                        <h5><i class="fas fa-construction"></i> Feature Coming Soon</h5>
                    </div>
                    <div class="card-body text-center">
                        <i class="fas fa-file-alt fa-5x text-secondary mb-4"></i>
                        <h4>Research Paper Search Tool</h4>
                        <p class="lead">This feature will integrate with paper-qa and provide:</p>
                        
                        <div class="row mt-4">
                            <div class="col-md-6">
                                <h6><i class="fas fa-search text-primary"></i> Search Capabilities:</h6>
                                <ul class="list-unstyled text-start">
                                    <li>• AI-powered paper search</li>
                                    <li>• Question answering from papers</li>
                                    <li>• Chemistry-specific queries</li>
                                    <li>• Multi-paper analysis</li>
                                    <li>• Citation network exploration</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-database text-info"></i> Data Sources:</h6>
                                <ul class="list-unstyled text-start">
                                    <li>• ArXiv chemistry papers</li>
                                    <li>• PubMed/PMC articles</li>
                                    <li>• Chemical journal databases</li>
                                    <li>• Patent databases</li>
                                    <li>• Preprint servers</li>
                                </ul>
                            </div>
                        </div>
                        
                        <div class="row mt-4">
                            <div class="col-12">
                                <h6><i class="fas fa-brain text-warning"></i> AI Features:</h6>
                                <div class="row">
                                    <div class="col-md-6">
                                        <ul class="list-unstyled text-start">
                                            <li>• Semantic paper search</li>
                                            <li>• Key finding extraction</li>
                                            <li>• Method comparison</li>
                                        </ul>
                                    </div>
                                    <div class="col-md-6">
                                        <ul class="list-unstyled text-start">
                                            <li>• Literature summarization</li>
                                            <li>• Research gap identification</li>
                                            <li>• Trend analysis</li>
                                        </ul>
                                    </div>
                                </div>
                            </div>
                        </div>
                        
                        <hr>
                        
                        <!-- Example Queries -->
                        <div class="mb-4">
                            <h6><i class="fas fa-lightbulb text-warning"></i> Example Research Queries:</h6>
                            <div class="row">
                                <div class="col-md-6">
                                    <div class="card border-light">
                                        <div class="card-body">
                                            <h6 class="card-title">Synthesis Research</h6>
                                            <p class="card-text small">"What are the latest methods for synthesizing perovskite solar cells?"</p>
                                        </div>
                                    </div>
                                </div>
                                <div class="col-md-6">
                                    <div class="card border-light">
                                        <div class="card-body">
                                            <h6 class="card-title">Drug Discovery</h6>
                                            <p class="card-text small">"Recent advances in COVID-19 drug targets and therapeutic compounds"</p>
                                        </div>
                                    </div>
                                </div>
                            </div>
                            <div class="row mt-2">
                                <div class="col-md-6">
                                    <div class="card border-light">
                                        <div class="card-body">
                                            <h6 class="card-title">Catalysis</h6>
                                            <p class="card-text small">"Comparison of heterogeneous catalysts for CO2 reduction"</p>
                                        </div>
                                    </div>
                                </div>
                                <div class="col-md-6">
                                    <div class="card border-light">
                                        <div class="card-body">
                                            <h6 class="card-title">Materials Science</h6>
                                            <p class="card-text small">"Properties and applications of graphene-based materials"</p>
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                        
                        <p class="text-muted">In the meantime, you can use the other available tools:</p>
                        <div class="d-grid gap-2 d-md-flex justify-content-center">
                            <a href="ai_assistant.php" class="btn btn-warning">
                                <i class="fas fa-robot"></i> AI Assistant
                            </a>
                            <a href="molecule_analyzer.php" class="btn btn-primary">
                                <i class="fas fa-atom"></i> Molecule Analyzer
                            </a>
                            <a href="format_converter.php" class="btn btn-success">
                                <i class="fas fa-exchange-alt"></i> Format Converter
                            </a>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Chem Paper Reference -->
        <div class="row mt-5">
            <div class="col-lg-12">
                <div class="card border-primary">
                    <div class="card-header bg-primary text-white">
                        <h5><i class="fas fa-graduation-cap"></i> Chem Research Paper</h5>
                    </div>
                    <div class="card-body">
                        <div class="row align-items-center">
                            <div class="col-md-8">
                                <h6>Chem: Augmenting large-language models with chemistry tools</h6>
                                <p class="text-muted mb-2">
                                    <strong>Authors:</strong> Andres M Bran, Sam Cox, Oliver Schilter, Carlo Baldassari, Andrew D White, Philippe Schwaller
                                </p>
                                <p class="text-muted mb-2">
                                    <strong>Abstract:</strong> This paper introduces Chem, an LLM chemistry agent designed to accomplish tasks across organic synthesis, drug discovery, and materials design...
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

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script src="../assets/js/main.js"></script>
</body>
</html>
