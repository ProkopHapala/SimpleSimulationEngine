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
                logger.info("Running script...");

                // Heuristic: Check if script uses 'mesh' or 'ConstructionBlockTests'
                if (code.includes('mesh.') || code.includes('ConstructionBlockTests')) {
                    try {
                        // Check if 'mesh' is declared in the code to avoid "Identifier 'mesh' has already been declared"
                        const declaresMesh = /^\s*(const|let|var)\s+mesh\b/m.test(code);
                        const meshArgName = declaresMesh ? '_mesh_ignored' : 'mesh';

                        // Run on Main Thread
                        const func = new Function('engine', meshArgName, 'Vec3', 'ConstructionBlockTests', 'window',
                            `"use strict";\n${code}`
                        );

                        // Clear mesh before running if it's a mesh test (usually desired)
                        // But maybe the script wants to add to it? 
                        // The tests usually call setup() which clears it.

                        func(this.engine, this.engine.mesh, Vec3, ConstructionBlockTests, window);

                        if (this.renderer) {
                            this.renderer.updateGeometry(this.engine.mesh);
                        }
                        logger.info("Script executed on Main Thread.");
                    } catch (e) {
                        logger.error(`Script Error: ${e.message}`);
                        logger.error(e);
                    }
                } else {
                    // Run in Worker (Simulation Script)
                    this.engine.runScript(code);
                }
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
        const zoomInput = document.getElementById('zoomInput');
        if (zoomInput) {
            zoomInput.addEventListener('input', (e) => {
                const val = parseFloat(e.target.value);
                if (!isNaN(val) && this.renderer && this.renderer.camera) {
                    this.renderer.camera.zoom = val;
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

        // --- Label Controls ---
        const selMode = document.getElementById('selLabelMode');
        const colLabel = document.getElementById('colLabel');
        const numSize = document.getElementById('numLabelSize');
        const chkFixed = document.getElementById('chkLabelFixed');

        if (selMode && colLabel && numSize) {
            selMode.addEventListener('change', () => {
                if (this.renderer) this.renderer.setLabelMode(selMode.value);
            });

            const updateStyle = () => {
                if (this.renderer) {
                    this.renderer.setLabelStyle(colLabel.value, parseFloat(numSize.value));
                }
            };

            colLabel.addEventListener('input', updateStyle);
            numSize.addEventListener('input', updateStyle);

            if (chkFixed) {
                chkFixed.addEventListener('change', () => {
                    if (this.renderer) this.renderer.setLabelScreenSpace(chkFixed.checked);
                });
            }

            // Set initial state
            if (this.renderer) {
                this.renderer.setLabelMode(selMode.value);
                updateStyle();
                if (chkFixed) this.renderer.setLabelScreenSpace(chkFixed.checked);
            }
        }

        // --- Camera Controls ---
        const selCam = document.getElementById('selCameraMode');
        if (selCam) {
            selCam.addEventListener('change', () => {
                if (this.renderer) this.renderer.setCameraMode(selCam.value);
            });
        }

        // --- Tests ---
        const selTest = document.getElementById('selTest');
        const btnRunTest = document.getElementById('btnRunTest');

        if (selTest && window.ConstructionBlockTests) {
            // Clear existing options first if any
            selTest.innerHTML = '';
            const defaultOpt = document.createElement('option');
            defaultOpt.value = "";
            defaultOpt.textContent = "-- Select Test --";
            selTest.appendChild(defaultOpt);

            logger.info(`Loaded ${Object.keys(ConstructionBlockTests.tests).length} tests: ${Object.keys(ConstructionBlockTests.tests).join(', ')}`);

            for (const testName in ConstructionBlockTests.tests) {
                const opt = document.createElement('option');
                opt.value = testName;
                opt.textContent = testName;
                selTest.appendChild(opt);
            }

            selTest.addEventListener('change', () => {
                const testName = selTest.value;
                logger.info(`Selected test: ${testName}`);
                if (!testName) return;

                const testFunc = ConstructionBlockTests.tests[testName];
                if (testFunc) {
                    try {
                        // Extract function body
                        let code = testFunc.toString();
                        // Remove function signature "function (engine) {" or "(engine) => {"
                        const bodyStart = code.indexOf('{');
                        const bodyEnd = code.lastIndexOf('}');

                        if (bodyStart !== -1 && bodyEnd !== -1) {
                            code = code.substring(bodyStart + 1, bodyEnd);

                            // Normalize indentation
                            const lines = code.split('\n');
                            // Find minimum indentation (ignoring empty lines)
                            let minIndent = Infinity;
                            for (const line of lines) {
                                if (line.trim().length > 0) {
                                    const indent = line.match(/^\s*/)[0].length;
                                    if (indent < minIndent) minIndent = indent;
                                }
                            }

                            if (minIndent !== Infinity && minIndent > 0) {
                                code = lines.map(line => {
                                    if (line.trim().length > 0) {
                                        return line.substring(minIndent);
                                    }
                                    return line;
                                }).join('\n');
                            }

                            // Trim and clean up
                            code = code.trim();
                            // Add standard header
                            code = `// Test: ${testName}\n// Context: engine, mesh, Vec3, ConstructionBlockTests\n\n${code}`;

                            scriptInput.value = code;
                            logger.info(`Populated script area for ${testName}`);
                        } else {
                            logger.error(`Could not parse function body for ${testName}`);
                        }
                    } catch (e) {
                        logger.error(`Error populating script: ${e.message}`);
                    }
                } else {
                    logger.error(`Test function not found for ${testName}`);
                }
            });
        }

        if (btnRunTest) {
            // btnRunTest is now redundant if we run from script, but let's keep it as a shortcut
            // Or better, make it run the CURRENT script in the box? 
            // The user request implies "put it into the SCRIPTING test area, and from there we can run it".
            // So the primary run button is 'run-script-btn'. 
            // We can hide btnRunTest or make it just populate and run.
            // Let's keep it simple: selecting populates. 'Run Script' runs.
            // We can remove btnRunTest listener or make it trigger the run-script logic.
            btnRunTest.style.display = 'none'; // Hide it as requested workflow is via Scripting area
        }



        // --- Verbosity ---
        const numVerbosityUI = document.getElementById('numVerbosityUI');
        const numVerbosityCon = document.getElementById('numVerbosityCon');

        if (numVerbosityUI && logger) {
            numVerbosityUI.value = logger.uiVerbosity;
            numVerbosityUI.addEventListener('change', () => {
                const val = parseInt(numVerbosityUI.value);
                if (!isNaN(val)) {
                    logger.setUIVerbosity(val);
                    logger.info(`UI Verbosity set to ${val}`);
                }
            });
        }

        if (numVerbosityCon && logger) {
            numVerbosityCon.value = logger.consoleVerbosity;
            numVerbosityCon.addEventListener('change', () => {
                const val = parseInt(numVerbosityCon.value);
                if (!isNaN(val)) {
                    logger.setConsoleVerbosity(val);
                    logger.info(`Console Verbosity set to ${val}`);
                }
            });
        }

        // --- Log ---
        const btnClearLog = document.getElementById('btnClearLog');
        if (btnClearLog) {
            btnClearLog.addEventListener('click', () => {
                if (logger) logger.clear();
            });
        }
    }
}
