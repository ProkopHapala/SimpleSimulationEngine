class GUI {
    constructor(engine, renderer) {
        this.engine = engine;
        this.renderer = renderer;

        this.sidebar = document.getElementById('side-panel');
        this.resizer = document.getElementById('gui-resizer');
        this.canvasContainer = document.getElementById('canvas-container');

        this.initResize();
        this.initControls();
    }

    initResize() {
        if (!this.resizer || !this.sidebar) return;

        let isResizing = false;

        this.resizer.addEventListener('mousedown', (e) => {
            isResizing = true;
            document.body.style.cursor = 'col-resize';
            e.preventDefault();
        });

        window.addEventListener('mousemove', (e) => {
            if (!isResizing) return;

            // Right-side panel logic
            const newWidth = window.innerWidth - e.clientX;

            if (newWidth < 200 || newWidth > 800) return; // Limits

            this.sidebar.style.width = `${newWidth}px`;
            this.resizer.style.right = `${newWidth}px`; // Move resizer with panel

            // Adjust canvas if needed (optional, but good for layout)
            // this.canvasContainer.style.width = `calc(100vw - ${newWidth}px)`;

            if (this.renderer) {
                this.renderer.resize();
            }
        });

        window.addEventListener('mouseup', () => {
            if (isResizing) {
                isResizing = false;
                document.body.style.cursor = 'default';
                if (this.renderer) this.renderer.resize();
            }
        });
    }

    initControls() {
        // --- Scripting ---
        const runBtn = document.getElementById('run-script-btn');
        const scriptInput = document.getElementById('script-input');
        const loadScriptBtn = document.getElementById('load-script-btn');
        const fileInput = document.getElementById('script-file-input');

        if (runBtn) {
            runBtn.addEventListener('click', () => {
                const code = scriptInput.value;
                window.logger.info("Running script...");
                this.engine.runScript(code);
            });
        }

        if (loadScriptBtn && fileInput) {
            loadScriptBtn.addEventListener('click', () => fileInput.click());
            fileInput.addEventListener('change', (e) => {
                const file = e.target.files[0];
                if (!file) return;
                const reader = new FileReader();
                reader.onload = (ev) => {
                    scriptInput.value = ev.target.result;
                };
                reader.readAsText(file);
            });
        }

        // --- View Controls ---
        const zoomSlider = document.getElementById('zoomSlider');
        if (zoomSlider) {
            zoomSlider.addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                const zoom = Math.pow(10, val);
                if (this.renderer && this.renderer.camera) {
                    this.renderer.camera.zoom = zoom;
                    this.renderer.camera.updateProjectionMatrix();
                }
            });
        }

        const setView = (pos, up) => {
            if (this.renderer && this.renderer.camera && this.renderer.controls) {
                const cam = this.renderer.camera;
                cam.position.set(pos[0], pos[1], pos[2]);
                cam.up.set(up[0], up[1], up[2]);
                cam.lookAt(0, 0, 0);
                cam.updateProjectionMatrix();
                this.renderer.controls.update();
            }
        };

        document.getElementById('btnViewXY')?.addEventListener('click', () => setView([0, 0, 20], [0, 1, 0]));
        document.getElementById('btnViewXZ')?.addEventListener('click', () => setView([0, 20, 0], [0, 0, -1]));
        document.getElementById('btnViewYZ')?.addEventListener('click', () => setView([20, 0, 0], [0, 1, 0]));

        // --- Tests ---
        const selTest = document.getElementById('selTest');
        const btnRunTest = document.getElementById('btnRunTest');

        if (selTest && window.ConstructionBlockTests) {
            // Clear existing options first if any
            selTest.innerHTML = '';
            for (const testName in ConstructionBlockTests.tests) {
                const opt = document.createElement('option');
                opt.value = testName;
                opt.textContent = testName;
                selTest.appendChild(opt);
            }
        }

        if (btnRunTest) {
            btnRunTest.addEventListener('click', () => {
                const testName = selTest.value;
                const testFunc = ConstructionBlockTests.tests[testName];
                if (testFunc) {
                    window.logger.info(`Running test: ${testName}`);
                    testFunc(this.engine);
                }
            });
        }

        // --- Verbosity ---
        const numVerbosity = document.getElementById('numVerbosity');
        if (numVerbosity) {
            numVerbosity.value = window.VERBOSITY_LEVEL;
            numVerbosity.addEventListener('change', () => {
                const val = parseInt(numVerbosity.value);
                if (!isNaN(val)) {
                    window.VERBOSITY_LEVEL = val;
                    window.logger.info(`Verbosity set to ${val}`);
                }
            });
        }

        // --- Log ---
        const btnClearLog = document.getElementById('btnClearLog');
        if (btnClearLog) {
            btnClearLog.addEventListener('click', () => {
                if (window.logger) window.logger.clear();
            });
        }
    }
}
