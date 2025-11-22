
// Global instances
window.engine = null;
window.renderer = null;

document.addEventListener('DOMContentLoaded', () => {
    console.log("Initializing Spacecraft Editor...");

    // 1. Initialize Engine
    window.engine = new SpaceCraftEngine();

    // 2. Initialize Renderer
    window.renderer = new SpaceCraftRenderer(window.engine);
    window.renderer.init(document.getElementById('canvas-container'));

    // 3. UI Bindings
    const runBtn = document.getElementById('run-script-btn');
    const scriptInput = document.getElementById('script-input');
    const verbositySelect = document.getElementById('verbosity-select');
    const loadScriptBtn = document.getElementById('load-script-btn');
    const fileInput = document.getElementById('script-file-input');

    runBtn.addEventListener('click', () => {
        const code = scriptInput.value;
        console.log("Running script...");
        window.engine.runScript(code);
    });

    verbositySelect.addEventListener('change', (e) => {
        const v = parseInt(e.target.value);
        window.engine.verbosity = v;
        console.log(`Verbosity set to ${v}`);
    });

    loadScriptBtn.addEventListener('click', () => {
        fileInput.click();
    });

    fileInput.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (!file) return;
        const reader = new FileReader();
        reader.onload = (e) => {
            scriptInput.value = e.target.result;
        };
        reader.readAsText(file);
    });

    // 4. Start Loop
    animate();
});

function animate() {
    requestAnimationFrame(animate);

    // Update Engine (Physics/Logic)
    // window.engine.update(); // TODO: Physics

    // Update Renderer
    window.renderer.update();
    window.renderer.render();
}

// Helper for logging to UI
window.logToUI = function (msg) {
    const logContent = document.getElementById('log-content');
    const line = document.createElement('div');
    line.textContent = `> ${msg}`;
    logContent.appendChild(line);
    logContent.scrollTop = logContent.scrollHeight;
    console.log(msg);
};
