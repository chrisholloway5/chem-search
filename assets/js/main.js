// Subchems Search - Main JavaScript Functions

document.addEventListener('DOMContentLoaded', function() {
    // Initialize tooltips
    var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="tooltip"]'));
    var tooltipList = tooltipTriggerList.map(function (tooltipTriggerEl) {
        return new bootstrap.Tooltip(tooltipTriggerEl);
    });

    // Initialize popovers
    var popoverTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle="popover"]'));
    var popoverList = popoverTriggerList.map(function (popoverTriggerEl) {
        return new bootstrap.Popover(popoverTriggerEl);
    });

    // Smooth scrolling for anchor links
    document.querySelectorAll('a[href^="#"]').forEach(anchor => {
        anchor.addEventListener('click', function (e) {
            e.preventDefault();
            const target = document.querySelector(this.getAttribute('href'));
            if (target) {
                target.scrollIntoView({
                    behavior: 'smooth',
                    block: 'start'
                });
            }
        });
    });

    // Auto-hide alerts after 5 seconds
    const alerts = document.querySelectorAll('.alert');
    alerts.forEach(alert => {
        if (alert.classList.contains('alert-success') || alert.classList.contains('alert-info')) {
            setTimeout(() => {
                const alertInstance = new bootstrap.Alert(alert);
                alertInstance.close();
            }, 5000);
        }
    });
});

// Utility Functions

// Copy text to clipboard
function copyToClipboard(text, button = null) {
    if (navigator.clipboard) {
        navigator.clipboard.writeText(text).then(function() {
            showCopySuccess(button);
        }, function(err) {
            console.error('Could not copy text: ', err);
            fallbackCopyTextToClipboard(text, button);
        });
    } else {
        fallbackCopyTextToClipboard(text, button);
    }
}

// Fallback copy function for older browsers
function fallbackCopyTextToClipboard(text, button = null) {
    const textArea = document.createElement("textarea");
    textArea.value = text;
    textArea.style.top = "0";
    textArea.style.left = "0";
    textArea.style.position = "fixed";
    
    document.body.appendChild(textArea);
    textArea.focus();
    textArea.select();
    
    try {
        const successful = document.execCommand('copy');
        if (successful) {
            showCopySuccess(button);
        }
    } catch (err) {
        console.error('Fallback: Oops, unable to copy', err);
    }
    
    document.body.removeChild(textArea);
}

// Show copy success feedback
function showCopySuccess(button) {
    if (button) {
        const originalText = button.innerHTML;
        const originalClass = button.className;
        
        button.innerHTML = '<i class="fas fa-check"></i> Copied!';
        button.className = button.className.replace(/btn-outline-\w+/, 'btn-success');
        
        setTimeout(() => {
            button.innerHTML = originalText;
            button.className = originalClass;
        }, 2000);
    } else {
        showToast('Copied to clipboard!', 'success');
    }
}

// Show toast notification
function showToast(message, type = 'info') {
    const toastContainer = getOrCreateToastContainer();
    const toastId = 'toast-' + Date.now();
    
    const toastHtml = `
        <div id="${toastId}" class="toast align-items-center text-white bg-${type} border-0" role="alert" aria-live="assertive" aria-atomic="true">
            <div class="d-flex">
                <div class="toast-body">
                    ${message}
                </div>
                <button type="button" class="btn-close btn-close-white me-2 m-auto" data-bs-dismiss="toast" aria-label="Close"></button>
            </div>
        </div>
    `;
    
    toastContainer.insertAdjacentHTML('beforeend', toastHtml);
    
    const toastElement = document.getElementById(toastId);
    const toast = new bootstrap.Toast(toastElement);
    toast.show();
    
    // Remove toast element after it's hidden
    toastElement.addEventListener('hidden.bs.toast', function() {
        toastElement.remove();
    });
}

// Get or create toast container
function getOrCreateToastContainer() {
    let container = document.getElementById('toast-container');
    if (!container) {
        container = document.createElement('div');
        container.id = 'toast-container';
        container.className = 'toast-container position-fixed bottom-0 end-0 p-3';
        container.style.zIndex = '1055';
        document.body.appendChild(container);
    }
    return container;
}

// Validate SMILES notation (basic validation)
function validateSMILES(smiles) {
    if (!smiles || smiles.trim() === '') {
        return false;
    }
    
    // Basic SMILES validation - check for valid characters
    const validChars = /^[A-Za-z0-9\[\]()@+\-=#$:\/\\\.%]+$/;
    return validChars.test(smiles.trim());
}

// Format molecular weight
function formatMolecularWeight(weight) {
    return parseFloat(weight).toFixed(2) + ' g/mol';
}

// Format LogP value
function formatLogP(logp) {
    return parseFloat(logp).toFixed(2);
}

// Auto-resize textarea
function autoResizeTextarea(textarea) {
    textarea.style.height = 'auto';
    textarea.style.height = textarea.scrollHeight + 'px';
}

// Initialize auto-resize for all textareas with class 'auto-resize'
document.addEventListener('DOMContentLoaded', function() {
    const textareas = document.querySelectorAll('textarea.auto-resize');
    textareas.forEach(textarea => {
        textarea.addEventListener('input', function() {
            autoResizeTextarea(this);
        });
        // Initial resize
        autoResizeTextarea(textarea);
    });
});

// Loading animation helpers
function showLoading(element) {
    if (element) {
        element.disabled = true;
        const originalText = element.innerHTML;
        element.setAttribute('data-original-text', originalText);
        element.innerHTML = '<span class="loading"></span> Loading...';
    }
}

function hideLoading(element) {
    if (element) {
        element.disabled = false;
        const originalText = element.getAttribute('data-original-text');
        if (originalText) {
            element.innerHTML = originalText;
            element.removeAttribute('data-original-text');
        }
    }
}

// Form validation helpers
function markFieldAsValid(field) {
    field.classList.remove('is-invalid');
    field.classList.add('is-valid');
}

function markFieldAsInvalid(field, message = '') {
    field.classList.remove('is-valid');
    field.classList.add('is-invalid');
    
    const feedback = field.parentNode.querySelector('.invalid-feedback');
    if (feedback && message) {
        feedback.textContent = message;
    }
}

// Example molecules data
const exampleMolecules = {
    'water': {
        name: 'Water',
        smiles: 'O',
        formula: 'H2O'
    },
    'ethanol': {
        name: 'Ethanol',
        smiles: 'CCO',
        formula: 'C2H6O'
    },
    'caffeine': {
        name: 'Caffeine',
        smiles: 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
        formula: 'C8H10N4O2'
    },
    'aspirin': {
        name: 'Aspirin',
        smiles: 'CC(=O)OC1=CC=CC=C1C(=O)O',
        formula: 'C9H8O4'
    },
    'tylenol': {
        name: 'Acetaminophen (Tylenol)',
        smiles: 'CC(=O)NC1=CC=C(C=C1)O',
        formula: 'C8H9NO2'
    },
    'ibuprofen': {
        name: 'Ibuprofen',
        smiles: 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O',
        formula: 'C13H18O2'
    }
};

// Get example molecule by key
function getExampleMolecule(key) {
    return exampleMolecules[key] || null;
}

// Fill form with example molecule
function fillExampleMolecule(key, smilesFieldId = 'smiles') {
    const molecule = getExampleMolecule(key);
    if (molecule) {
        const field = document.getElementById(smilesFieldId);
        if (field) {
            field.value = molecule.smiles;
            field.focus();
            showToast(`Filled with ${molecule.name}`, 'info');
        }
    }
}

// Export functions for global use
window.SubchemsWeb = {
    copyToClipboard,
    showToast,
    validateSMILES,
    formatMolecularWeight,
    formatLogP,
    autoResizeTextarea,
    showLoading,
    hideLoading,
    markFieldAsValid,
    markFieldAsInvalid,
    getExampleMolecule,
    fillExampleMolecule,
    exampleMolecules
};
