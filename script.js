// Preview data loading removed - section no longer exists

// Load and display scores from task_scores.csv or calculate from submission.csv
async function loadScores() {
    try {
        // Try to load from task_scores.csv first
        const response = await fetch('task_scores.csv');
        const text = await response.text();
        const lines = text.split('\n').filter(line => line.trim());
        
        // Skip header
        const dataRows = lines.slice(1);
        const scores = [];
        let totalScore = 0;
        let totalItems = 0;
        let taskCount = 0;
        
        dataRows.forEach(line => {
            const [taskId, itemsCount, score] = line.split(',');
            if (taskId && taskId !== 'TOTAL') {
                scores.push({
                    taskId: taskId,
                    itemsCount: parseInt(itemsCount),
                    score: parseFloat(score)
                });
                totalScore += parseFloat(score);
                totalItems += parseInt(itemsCount);
                taskCount++;
            } else if (taskId === 'TOTAL') {
                // Use total from CSV if available
                if (!totalScore && score) totalScore = parseFloat(score);
                if (!totalItems && itemsCount) totalItems = parseInt(itemsCount);
            }
        });
        
        // Display scores
        displayScores(scores, totalScore, totalItems, taskCount);
        
    } catch (error) {
        console.log('task_scores.csv not found, calculating from submission.csv...');
        // If task_scores.csv doesn't exist, calculate from submission.csv
        calculateScoresFromSubmission();
    }
}

// Calculate scores from submission.csv
async function calculateScoresFromSubmission() {
    try {
        const response = await fetch('submission.csv');
        const text = await response.text();
        const lines = text.split('\n').filter(line => line.trim());
        
        // Skip header
        const dataRows = lines.slice(1);
        const taskScores = {};
        let totalItems = 0;
        
        // Group by task and count items
        dataRows.forEach(line => {
            const [id, x, y, deg] = line.split(',');
            if (id && id.includes('_')) {
                const taskId = id.split('_')[0];
                if (!taskScores[taskId]) {
                    taskScores[taskId] = 0;
                }
                taskScores[taskId]++;
                totalItems++;
            }
        });
        
        // Convert to array format
        const scores = Object.keys(taskScores)
            .sort((a, b) => parseInt(a) - parseInt(b))
            .map(taskId => ({
                taskId: taskId,
                itemsCount: taskScores[taskId],
                score: taskScores[taskId] // Simple metric: count = score
            }));
        
        const totalScore = scores.reduce((sum, s) => sum + s.score, 0);
        const taskCount = scores.length;
        
        displayScores(scores, totalScore, totalItems, taskCount);
        
    } catch (error) {
        console.error('Error calculating scores:', error);
        displayScoreError();
    }
}

// Display scores in the UI
function displayScores(scores, totalScore, totalItems, taskCount) {
    // Update summary stats
    const totalScoreDisplay = document.getElementById('total-score-display');
    
    // Format score with 12 decimal places
    const formatScore12 = (score) => {
        return score.toFixed(12);
    };
    
    if (totalScoreDisplay) {
        totalScoreDisplay.textContent = formatScore12(totalScore);
    }
    
    // Display scores table (show all tasks)
    const scoresBody = document.getElementById('scores-body');
    if (scoresBody) {
        scoresBody.innerHTML = '';
        
        // Store all scores for filtering
        scoresBody.dataset.allScores = JSON.stringify(scores);
        
        // Display all tasks
        renderScores(scores, totalScore, scoresBody);
        
        // Set up filter
        const filterInput = document.getElementById('task-filter');
        if (filterInput) {
            filterInput.addEventListener('input', (e) => {
                const filterValue = e.target.value.toLowerCase().trim();
                const allScores = JSON.parse(scoresBody.dataset.allScores || '[]');
                
                if (filterValue === '') {
                    renderScores(allScores, totalScore, scoresBody, true);
                } else {
                    const filtered = filterScores(allScores, filterValue);
                    // Calculate filtered total
                    const filteredTotal = filtered.reduce((sum, s) => sum + s.score, 0);
                    renderScores(filtered, filteredTotal, scoresBody, true);
                }
            });
        }
    }
}

// Render scores in the table
function renderScores(scores, totalScore, scoresBody, showTotal = true) {
    scoresBody.innerHTML = '';
    
    // Display all filtered tasks
    for (let i = 0; i < scores.length; i++) {
        const score = scores[i];
        const row = document.createElement('tr');
        row.className = 'task-row';
        row.setAttribute('data-task-id', score.taskId);
        row.innerHTML = `
            <td><strong>${score.taskId}</strong></td>
            <td>${score.score.toFixed(12)}</td>
        `;
        scoresBody.appendChild(row);
    }
    
    // Add total row if showing all tasks (not filtered)
    if (showTotal) {
        const totalRow = document.createElement('tr');
        totalRow.className = 'total-row';
        totalRow.style.backgroundColor = '#f0f0f0';
        totalRow.style.fontWeight = 'bold';
        totalRow.innerHTML = `
            <td><strong>TOTAL</strong></td>
            <td><strong>${totalScore.toFixed(12)}</strong></td>
        `;
        scoresBody.appendChild(totalRow);
    }
}

// Filter scores based on filter input
function filterScores(scores, filterValue) {
    if (!filterValue) return scores;
    
    // Check for range (e.g., "100-200")
    const rangeMatch = filterValue.match(/(\d+)-(\d+)/);
    if (rangeMatch) {
        const start = parseInt(rangeMatch[1]);
        const end = parseInt(rangeMatch[2]);
        return scores.filter(score => {
            const taskId = parseInt(score.taskId);
            return taskId >= start && taskId <= end;
        });
    }
    
    // Check for multiple IDs separated by comma or space
    const ids = filterValue.split(/[,\s]+/).filter(id => id);
    if (ids.length > 1) {
        return scores.filter(score => {
            return ids.some(id => score.taskId.toLowerCase().includes(id));
        });
    }
    
    // Single filter - match task ID
    return scores.filter(score => {
        return score.taskId.toLowerCase().includes(filterValue);
    });
}

// Display error message
function displayScoreError() {
    const scoresBody = document.getElementById('scores-body');
    if (scoresBody) {
        scoresBody.innerHTML = `
            <tr>
                <td colspan="2" style="text-align: center; padding: 20px; color: #666;">
                    Unable to load scores. Please check the files.
                </td>
            </tr>
        `;
    }
}

// Initialize when page loads
document.addEventListener('DOMContentLoaded', () => {
    loadScores();
});
