<?php
session_start();
require_once '../includes/config.php';

// Handle form submission
$result = null;
$error = null;

if ($_SERVER['REQUEST_METHOD'] === 'POST') {
    $question = $_POST['question'] ?? '';
    
    if (!empty($question)) {
        try {
            $result = askChem($question);
        } catch (Exception $e) {
            $error = $e->getMessage();
        }
    } else {
        $error = "Please enter a chemistry question.";
    }
}

function askChem($question) {
    $pythonScript = realpath('../python/chem_assistant.py');
    $pythonPath = 'C:/xampp/htdocs/chemcrow/project/.venv/Scripts/python.exe';
    
    $escapedQuestion = escapeshellarg($question);
    
    $command = "$pythonPath \"$pythonScript\" $escapedQuestion";
    $output = shell_exec($command);
    
    if ($output === null) {
        throw new Exception("Failed to execute Chem assistant script");
    }
    
    $data = json_decode($output, true);
    
    if (json_last_error() !== JSON_ERROR_NONE) {
        throw new Exception("Invalid JSON response from Chem assistant");
    }
    
    if (isset($data['error'])) {
        throw new Exception($data['error']);
    }
    
    return $data;
}
?>

<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AI Chemistry Assistant - Chem Search</title>
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
                    <i class="fas fa-robot text-warning"></i> AI Chemistry Assistant
                </h1>
                <p class="lead">Ask complex chemistry questions powered by Chem' AI agent system.</p>
            </div>
        </div>

        <!-- Question Form -->
        <div class="row mb-4">
            <div class="col-lg-10 mx-auto">
                <div class="card shadow">
                    <div class="card-body">
                        <form method="POST" class="needs-validation" novalidate>
                            <div class="mb-3">
                                <label for="question" class="form-label">Your Chemistry Question</label>
                                <textarea class="form-control" id="question" name="question" rows="4" 
                                          placeholder="Ask me anything about chemistry..."
                                          required><?php echo htmlspecialchars($_POST['question'] ?? ''); ?></textarea>
                                <div class="invalid-feedback">
                                    Please enter a chemistry question.
                                </div>
                            </div>
                            <div class="d-grid">
                                <button type="submit" class="btn btn-warning btn-lg">
                                    <i class="fas fa-brain"></i> Ask Chem AI
                                </button>
                            </div>
                        </form>
                    </div>
                </div>
            </div>
        </div>

        <!-- Example Questions -->
        <div class="row mb-4">
            <div class="col-lg-12">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-lightbulb text-warning"></i> Example Questions</h5>
                    </div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-6">
                                <h6>Molecular Properties:</h6>
                                <ul class="list-unstyled">
                                    <li><a href="#" class="text-decoration-none" onclick="fillQuestion('What is the molecular weight of caffeine?')">• What is the molecular weight of caffeine?</a></li>
                                    <li><a href="#" class="text-decoration-none" onclick="fillQuestion('Calculate the LogP of aspirin')">• Calculate the LogP of aspirin</a></li>
                                    <li><a href="#" class="text-decoration-none" onclick="fillQuestion('How many hydrogen bond donors does ibuprofen have?')">• How many hydrogen bond donors does ibuprofen have?</a></li>
                                </ul>
                            </div>
                            <div class="col-md-6">
                                <h6>Chemical Analysis:</h6>
                                <ul class="list-unstyled">
                                    <li><a href="#" class="text-decoration-none" onclick="fillQuestion('Compare the toxicity of methanol vs ethanol')">• Compare the toxicity of methanol vs ethanol</a></li>
                                    <li><a href="#" class="text-decoration-none" onclick="fillQuestion('What are the side effects of acetaminophen?')">• What are the side effects of acetaminophen?</a></li>
                                    <li><a href="#" class="text-decoration-none" onclick="fillQuestion('Explain the mechanism of aspirin')">• Explain the mechanism of aspirin</a></li>
                                </ul>
                            </div>
                        </div>
                        
                        <!-- API Status Warning -->
                        <div class="alert alert-warning mt-3" role="alert">
                            <i class="fas fa-exclamation-triangle"></i> 
                            <strong>Note:</strong> The AI assistant requires a valid OpenAI API key with available credits. 
                            <a href="../api_status.php" class="alert-link">Check API status</a>.
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <?php if ($error): ?>
        <div class="alert alert-danger" role="alert">
            <i class="fas fa-exclamation-triangle"></i> <strong>Error:</strong> <?php echo htmlspecialchars($error); ?>
            
            <?php if (strpos($error, 'quota') !== false || strpos($error, 'billing') !== false): ?>
            <hr>
            <h6>API Quota Issue:</h6>
            <p class="mb-0">Your OpenAI API key has exceeded its quota or has billing issues. Please:</p>
            <ul class="mb-0">
                <li>Check your OpenAI account at <a href="https://platform.openai.com/" target="_blank">platform.openai.com</a></li>
                <li>Add billing information or purchase credits</li>
                <li>Wait for quota reset if using free tier</li>
            </ul>
            <?php elseif (strpos($error, 'access') !== false): ?>
            <hr>
            <h6>Model Access Issue:</h6>
            <p class="mb-0">Your API key doesn't have access to the requested model. Consider upgrading your OpenAI plan.</p>
            <?php endif; ?>
        </div>
        <?php endif; ?>

        <?php if ($result): ?>
        <!-- AI Response -->
        <div class="row">
            <div class="col-lg-12">
                <div class="card border-warning shadow fade-in">
                    <div class="card-header bg-warning text-dark">
                        <h5><i class="fas fa-robot"></i> Chem AI Response</h5>
                    </div>
                    <div class="card-body">
                        <div class="mb-3">
                            <h6><i class="fas fa-question-circle text-primary"></i> Your Question:</h6>
                            <div class="bg-light p-3 rounded">
                                <?php echo nl2br(htmlspecialchars($_POST['question'])); ?>
                            </div>
                        </div>
                        
                        <div class="mb-3">
                            <h6><i class="fas fa-brain text-warning"></i> AI Answer:</h6>
                            <div class="bg-warning bg-opacity-10 p-3 rounded">
                                <?php echo nl2br(htmlspecialchars($result['answer'])); ?>
                            </div>
                        </div>
                        
                        <?php if (isset($result['tools_used']) && !empty($result['tools_used'])): ?>
                        <div class="mb-3">
                            <h6><i class="fas fa-tools text-info"></i> Tools Used:</h6>
                            <div class="d-flex flex-wrap gap-2">
                                <?php foreach ($result['tools_used'] as $tool): ?>
                                <span class="badge bg-info"><?php echo htmlspecialchars($tool); ?></span>
                                <?php endforeach; ?>
                            </div>
                        </div>
                        <?php endif; ?>
                        
                        <?php if (isset($result['processing_time'])): ?>
                        <div class="text-muted small">
                            <i class="fas fa-clock"></i> Processing time: <?php echo number_format($result['processing_time'], 2); ?> seconds
                        </div>
                        <?php endif; ?>
                    </div>
                </div>
            </div>
        </div>
        <?php endif; ?>
    </main>

    <?php include '../includes/footer.php'; ?>

    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>
    <script>
        // Form validation
        (function() {
            'use strict';
            window.addEventListener('load', function() {
                var forms = document.getElementsByClassName('needs-validation');
                var validation = Array.prototype.filter.call(forms, function(form) {
                    form.addEventListener('submit', function(event) {
                        if (form.checkValidity() === false) {
                            event.preventDefault();
                            event.stopPropagation();
                        }
                        form.classList.add('was-validated');
                    }, false);
                });
            }, false);
        })();

        // Fill question function
        function fillQuestion(question) {
            document.getElementById('question').value = question;
            document.getElementById('question').focus();
        }
    </script>
</body>
</html>
