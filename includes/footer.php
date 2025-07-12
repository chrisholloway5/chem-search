<?php
// Determine the base path relative to current file
$current_dir = dirname($_SERVER['SCRIPT_NAME']);
$base_path = '';

// If we're in a subdirectory (like tools/), set base path to go up one level
if (basename($current_dir) === 'tools') {
    $base_path = '../';
}
?>
<footer class="bg-dark text-light py-4 mt-5">
    <div class="container">
        <div class="row">
            <div class="col-md-6">
                <h5 class="text-white"><i class="fas fa-flask text-warning"></i> Chem Search</h5>
                <p class="text-white">AI-powered chemistry platform built with Chem, RDKit, and modern web technologies.</p>
                <p class="text-white">&copy; <?php echo date('Y'); ?> Chem Search. Built for chemistry research and education.</p>
            </div>
            <div class="col-md-3">
                <h6>Tools</h6>
                <ul class="list-unstyled">
                    <li><a href="<?php echo $base_path; ?>tools/molecule_analyzer.php" class="text-light">Molecule Analyzer</a></li>
                    <li><a href="<?php echo $base_path; ?>tools/format_converter.php" class="text-light">Format Converter</a></li>
                    <li><a href="<?php echo $base_path; ?>tools/ai_assistant.php" class="text-light">AI Assistant</a></li>
                    <li><a href="<?php echo $base_path; ?>tools/database_search.php" class="text-light">Database Search</a></li>
                </ul>
            </div>
            <div class="col-md-3">
                <h6>Resources</h6>
                <ul class="list-unstyled">
                    <li><a href="<?php echo $base_path; ?>about.php" class="text-light">About Chem</a></li>
                    <li><a href="<?php echo $base_path; ?>api_status.php" class="text-light">API Status</a></li>
                    <li><a href="https://github.com/ur-whitelab/chemcrow-public" target="_blank" class="text-light">
                        <i class="fab fa-github"></i> GitHub
                    </a></li>
                    <li><a href="https://arxiv.org/abs/2304.05376" target="_blank" class="text-light">
                        <i class="fas fa-file-pdf"></i> Research Paper
                    </a></li>
                </ul>
            </div>
        </div>
    </div>
</footer>
