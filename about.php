<?php
session_start();
require_once 'includes/config.php';
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>About - Chem Search</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
    <link href="assets/css/style.css" rel="stylesheet">
</head>
<body>
    <?php include 'includes/header.php'; ?>
    
    <main class="container my-5">
        <!-- Header -->
        <div class="row mb-5">
            <div class="col-lg-12 text-center">
                <h1 class="mb-4">
                    <i class="fas fa-flask text-primary"></i> About Chem Search
                </h1>
                <p class="lead">AI-powered chemistry platform built with modern web technologies</p>
            </div>
        </div>

        <!-- About Chem Search -->
        <div class="row mb-5">
            <div class="col-lg-8 mx-auto">
                <div class="card shadow">
                    <div class="card-body">
                        <h3 class="text-center mb-4">What is Chem Search?</h3>
                        <p>
                            Chem is an open-source chemistry-focused agent designed to accomplish tasks across 
                            organic synthesis, drug discovery, and materials design. It integrates large language models (LLMs) 
                            with chemistry tools to provide expert-level chemical analysis and reasoning.
                        </p>
                        <p>
                            Built with LangChain, Chem Search uses a collection of chemical tools including RDKit, paper-qa, 
                            as well as relevant databases in chemistry like PubChem and chem-space.
                        </p>
                        <div class="text-center mt-4">
                            <a href="https://github.com/chrisholloway5/chem-search" target="_blank" class="btn btn-primary me-2">
                                <i class="fab fa-github"></i> GitHub Repository
                            </a>
                            <a href="https://arxiv.org/abs/2304.05376" target="_blank" class="btn btn-outline-primary">
                                <i class="fas fa-file-pdf"></i> Research Paper
                            </a>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Features -->
        <div class="row mb-5">
            <div class="col-lg-12">
                <h3 class="text-center mb-4">Key Features</h3>
                <div class="row g-4">
                    <div class="col-md-6 col-lg-3">
                        <div class="card h-100 text-center">
                            <div class="card-body">
                                <i class="fas fa-calculator fa-3x text-primary mb-3"></i>
                                <h6>Molecular Analysis</h6>
                                <p class="small">Calculate molecular properties, descriptors, and chemical characteristics.</p>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-3">
                        <div class="card h-100 text-center">
                            <div class="card-body">
                                <i class="fas fa-robot fa-3x text-warning mb-3"></i>
                                <h6>AI Assistant</h6>
                                <p class="small">Ask complex chemistry questions and get expert-level answers.</p>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-3">
                        <div class="card h-100 text-center">
                            <div class="card-body">
                                <i class="fas fa-exchange-alt fa-3x text-success mb-3"></i>
                                <h6>Format Conversion</h6>
                                <p class="small">Convert between SMILES, InChI, molecular formulas, and more.</p>
                            </div>
                        </div>
                    </div>
                    <div class="col-md-6 col-lg-3">
                        <div class="card h-100 text-center">
                            <div class="card-body">
                                <i class="fas fa-search fa-3x text-info mb-3"></i>
                                <h6>Database Integration</h6>
                                <p class="small">Access chemical databases and research literature.</p>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Technology Stack -->
        <div class="row mb-5">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header bg-primary text-white">
                        <h5><i class="fas fa-code"></i> Technology Stack</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h6><i class="fas fa-server text-primary"></i> Backend</h6>
                                <ul class="list-unstyled">
                                    <li><i class="fab fa-python text-info"></i> Python (ChemCrow, RDKit)</li>
                                    <li><i class="fab fa-php text-purple"></i> PHP (Web Interface)</li>
                                    <li><i class="fas fa-brain text-warning"></i> OpenAI GPT Models</li>
                                    <li><i class="fas fa-link text-success"></i> LangChain Framework</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-desktop text-success"></i> Frontend</h6>
                                <ul class="list-unstyled">
                                    <li><i class="fab fa-bootstrap text-purple"></i> Bootstrap 5</li>
                                    <li><i class="fab fa-js text-warning"></i> JavaScript (ES6+)</li>
                                    <li><i class="fab fa-css3 text-info"></i> CSS3 with Custom Styling</li>
                                    <li><i class="fas fa-mobile-alt text-danger"></i> Responsive Design</li>
                                </ul>
                            </div>
                        </div>
                        <div class="row mt-3">
                            <div class="col-md-6">
                                <h6><i class="fas fa-flask text-danger"></i> Chemistry Libraries</h6>
                                <ul class="list-unstyled">
                                    <li><i class="fas fa-atom text-primary"></i> RDKit (Cheminformatics)</li>
                                    <li><i class="fas fa-database text-info"></i> PubChem Integration</li>
                                    <li><i class="fas fa-file-alt text-secondary"></i> Paper-QA (Literature)</li>
                                    <li><i class="fas fa-search text-warning"></i> Chemical Databases</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-cloud text-info"></i> APIs & Services</h6>
                                <ul class="list-unstyled">
                                    <li><i class="fas fa-robot text-success"></i> OpenAI API</li>
                                    <li><i class="fas fa-search text-warning"></i> SerpAPI (Optional)</li>
                                    <li><i class="fas fa-globe text-primary"></i> Web-based Interface</li>
                                    <li><i class="fas fa-mobile text-danger"></i> Mobile-friendly</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Research Paper -->
        <div class="row mb-5">
            <div class="col-lg-12">
                <div class="card border-warning">
                    <div class="card-header bg-warning text-dark">
                        <h5><i class="fas fa-graduation-cap"></i> Research & Citation</h5>
                    </div>
                    <div class="card-body">
                        <h6>ChemCrow: Augmenting large-language models with chemistry tools</h6>
                        <p class="text-muted">
                            <strong>Authors:</strong> Andres M Bran, Sam Cox, Oliver Schilter, Carlo Baldassari, Andrew D White, Philippe Schwaller<br>
                            <strong>Published:</strong> arXiv preprint arXiv:2304.05376 (2023)
                        </p>
                        
                        <div class="bg-light p-3 rounded">
                            <h6>BibTeX Citation:</h6>
                            <code class="small">
@article{bran2023chemcrow,<br>
&nbsp;&nbsp;title={ChemCrow: Augmenting large-language models with chemistry tools},<br>
&nbsp;&nbsp;author={Andres M Bran and Sam Cox and Oliver Schilter and Carlo Baldassari and Andrew D White and Philippe Schwaller},<br>
&nbsp;&nbsp;year={2023},<br>
&nbsp;&nbsp;eprint={2304.05376},<br>
&nbsp;&nbsp;archivePrefix={arXiv},<br>
&nbsp;&nbsp;primaryClass={physics.chem-ph},<br>
&nbsp;&nbsp;publisher={arXiv}<br>
}
                            </code>
                        </div>
                        
                        <div class="mt-3">
                            <a href="https://arxiv.org/abs/2304.05376" target="_blank" class="btn btn-outline-warning">
                                <i class="fas fa-external-link-alt"></i> Read Full Paper
                            </a>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Developer -->
        <div class="row mb-5">
            <div class="col-lg-12">
                <div class="card border-primary">
                    <div class="card-header bg-primary text-white">
                        <h5><i class="fas fa-user-cog"></i> Developer</h5>
                    </div>
                    <div class="card-body">
                        <div class="row align-items-center">
                            <div class="col-md-8">
                                <h6><i class="fas fa-code text-primary"></i> Platform Developer</h6>
                                <p class="mb-2">
                                    <strong><i class="fas fa-user text-warning"></i> Christopher Holloway</strong><br>
                                    <span class="text-muted"><i class="fas fa-building text-info"></i> Progressive Robot Ltd</span>
                                </p>
                                <p class="small">
                                    Full-stack developer specializing in AI-powered web applications and chemistry platforms. 
                                    Passionate about bridging the gap between complex scientific tools and user-friendly interfaces.
                                </p>
                                <div class="small text-muted">
                                    <i class="fas fa-envelope"></i> Available for custom chemistry platform development<br>
                                    <i class="fas fa-globe"></i> Expertise in Python, PHP, AI integration, and cheminformatics
                                </div>
                            </div>
                            <div class="col-md-4 text-center">
                                <div class="bg-light p-3 rounded">
                                    <i class="fab fa-github fa-3x text-dark mb-2"></i>
                                    <h6>Open Source</h6>
                                    <a href="https://github.com/chrisholloway5/chem-search" target="_blank" class="btn btn-outline-primary">
                                        <i class="fab fa-github"></i> GitHub Profile
                                    </a>
                                </div>
                            </div>
                        </div>
                        
                        <hr class="my-3">
                        
                        <div class="row">
                            <div class="col-md-6">
                                <h6><i class="fas fa-tools text-success"></i> Development Stack</h6>
                                <ul class="list-unstyled small">
                                    <li><i class="fab fa-python text-info"></i> Python & Scientific Computing</li>
                                    <li><i class="fab fa-php text-purple"></i> PHP Web Development</li>
                                    <li><i class="fas fa-brain text-warning"></i> AI/ML Integration</li>
                                    <li><i class="fas fa-flask text-danger"></i> Chemistry Platform Design</li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-handshake text-info"></i> Services</h6>
                                <ul class="list-unstyled small">
                                    <li><i class="fas fa-laptop-code text-primary"></i> Custom Chemistry Platforms</li>
                                    <li><i class="fas fa-robot text-warning"></i> AI Tool Integration</li>
                                    <li><i class="fas fa-mobile-alt text-success"></i> Responsive Web Design</li>
                                    <li><i class="fas fa-database text-info"></i> Scientific Data Management</li>
                                </ul>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- License & Disclaimer -->
        <div class="row">
            <div class="col-lg-12">
                <div class="card border-secondary">
                    <div class="card-header bg-secondary text-white">
                        <h5><i class="fas fa-balance-scale"></i> License & Disclaimer</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h6><i class="fas fa-certificate text-success"></i> License</h6>
                                <p class="small">
                                    Chem Search is released under the MIT License. This web interface is built just because I could.
                                </p>
                            </div>
                            <div class="col-md-6">
                                <h6><i class="fas fa-exclamation-triangle text-warning"></i> Disclaimer</h6>
                                <p class="small">
                                    This tool is for educational and research purposes only. Results should be 
                                    verified independently before use in any critical applications.
                                </p>
                            </div>
                        </div>
                        
                        <hr>
                        
                        <div class="text-center">
                            <p class="mb-0">
                                <strong>Chem Search Web</strong> - Version <?php echo APP_VERSION; ?> | 
                                Built with ❤️ for the chemistry community
                            </p>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    </main>

    <?php include 'includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script src="assets/js/main.js"></script>
</body>
</html>
