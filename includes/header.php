<?php
// Determine the base path relative to current file
$current_dir = dirname($_SERVER['SCRIPT_NAME']);
$base_path = '';

// If we're in a subdirectory (like tools/), set base path to go up one level
if (basename($current_dir) === 'tools') {
    $base_path = '../';
}
?>
<nav class="navbar navbar-expand-lg navbar-dark bg-dark">
    <div class="container">
        <a class="navbar-brand" href="<?php echo $base_path; ?>index.php">
            <i class="fas fa-flask text-warning"></i> Chem Search
        </a>
        <button class="navbar-toggler" type="button" data-bs-toggle="collapse" data-bs-target="#navbarNav">
            <span class="navbar-toggler-icon"></span>
        </button>
        <div class="collapse navbar-collapse" id="navbarNav">
            <ul class="navbar-nav me-auto">
                <li class="nav-item">
                    <a class="nav-link" href="<?php echo $base_path; ?>index.php">
                        <i class="fas fa-home"></i> Home
                    </a>
                </li>
                <li class="nav-item dropdown">
                    <a class="nav-link dropdown-toggle" href="#" role="button" data-bs-toggle="dropdown">
                        <i class="fas fa-tools"></i> Tools
                    </a>
                    <ul class="dropdown-menu">
                        <li><a class="dropdown-item" href="<?php echo $base_path; ?>tools/molecule_analyzer.php">
                            <i class="fas fa-atom"></i> Molecule Analyzer
                        </a></li>
                        <li><a class="dropdown-item" href="<?php echo $base_path; ?>tools/format_converter.php">
                            <i class="fas fa-exchange-alt"></i> Format Converter
                        </a></li>
                        <li><a class="dropdown-item" href="<?php echo $base_path; ?>tools/ai_assistant.php">
                            <i class="fas fa-robot"></i> AI Assistant
                        </a></li>
                        <li><a class="dropdown-item" href="<?php echo $base_path; ?>tools/database_search.php">
                            <i class="fas fa-search"></i> Database Search
                        </a></li>
                        <li><a class="dropdown-item" href="<?php echo $base_path; ?>tools/structure_viewer.php">
                            <i class="fas fa-project-diagram"></i> Structure Viewer
                        </a></li>
                        <li><a class="dropdown-item" href="<?php echo $base_path; ?>tools/paper_search.php">
                            <i class="fas fa-file-alt"></i> Paper Search
                        </a></li>
                    </ul>
                </li>
                <li class="nav-item">
                    <a class="nav-link" href="<?php echo $base_path; ?>about.php">
                        <i class="fas fa-info-circle"></i> About
                    </a>
                </li>
            </ul>
            <ul class="navbar-nav">
                <li class="nav-item">
                    <a class="nav-link" href="<?php echo $base_path; ?>api_status.php">
                        <i class="fas fa-server"></i> API Status
                    </a>
                </li>
            </ul>
        </div>
    </div>
</nav>
