// Basic demo interactions: not connected to backend.
// Put this file in frontend/script.js

document.addEventListener("DOMContentLoaded", () => {
  const runBtn = document.getElementById("runBtn");
  const explainBtn = document.getElementById("explainBtn");
  const explainText = document.getElementById("explainText");
  const mechanism = document.getElementById("mechanism");
  const targetField = document.getElementById("target");
  const startingField = document.getElementById("starting");
  const reagentsField = document.getElementById("reagents");

  runBtn.addEventListener("click", () => {
    const target = targetField.value || "C1=CC=CC=C1";
    const starting = startingField.value || "CCO";
    const reagents = reagentsField.value || "H+, heat";

    // Simple demo output - replace with backend fetch when ready
    mechanism.innerHTML = `
      <h3>Generated (demo) Pathways</h3>
      <ol>
        <li><strong>Step 1:</strong> ${starting} → intermediate A (via ${reagents})</li>
        <li><strong>Step 2:</strong> intermediate A → intermediate B</li>
        <li><strong>Step 3:</strong> intermediate B → <em>${target}</em></li>
      </ol>
      <p style="color:#666">This is a static demo. Connect the frontend to your backend to generate real pathways.</p>
    `;

    explainText.textContent = "Click 'Explain Step' to get a short explanation of Step 1.";
  });

  explainBtn.addEventListener("click", () => {
    explainText.textContent = "Step 1: Protonation of the carbonyl followed by nucleophilic addition. (Demo text — replace with backend explanation.)";
  });

  // Better UX: run demo if fields are blank on first load
  // runBtn.click();
});
