import { logger } from '../../common_js/Logger.js';
import { Vec3 } from '../../common_js/Vec3.js';
import { SDfuncs } from '../../common_js/SDfuncs.js';
import { GUIutils } from '../../common_js/GUIutils.js';
import { ConstructionBlockTests } from './constructionBlockTests.js';
import { MeshGenTestGUI } from './MeshGenTestGUI.js';

export class GUI {
    constructor(engine, renderer) {
        this.engine = engine;
        this.renderer = renderer;

        this.sidebar = document.getElementById('side-panel');
        this.resizer = document.getElementById('gui-resizer');
        this.canvasContainer = document.getElementById('canvas-container');

        // Create screen-space selection rectangle (yellow box) for box selection.
        // IMPORTANT: do this BEFORE wiring event handlers so initControls /
        // initWidgetCallbacks capture a valid element in their closures.
        this.selectionBox = document.createElement('div');
        this.selectionBox.className = 'selection-box';
        document.body.appendChild(this.selectionBox);

        this.initResize();
        this.initControls();

        // DEBUG: draw a fixed test rectangle so we can verify rendering at all.
        // This should show a 100x100 yellow box near the top-left on load.
        // this.selectionBox.style.display = 'block';
        // this.selectionBox.style.width = '100px';
        // this.selectionBox.style.height = '100px';
        // this.selectionBox.style.transform = 'translate(100px, 100px)';

        // this.initMeshGenControls(); // Refactored to MeshGenTestGUI
        if (window.MeshGenTestGUI) {
            console.log("GUI: Instantiating MeshGenTestGUI...");
            this.meshGenGUI = new MeshGenTestGUI(this.engine, this.renderer);
        } else {
            console.error("GUI: MeshGenTestGUI not found in window!");
        }

        // --- Render Layers (vertices, edges, labels, selection) + sizes ---
        const chkRenderVerts     = document.getElementById('chkRenderVerts');
        const chkRenderEdges     = document.getElementById('chkRenderEdges');
        const chkRenderLabelsVis = document.getElementById('chkRenderLabels');
        const chkRenderSelection = document.getElementById('chkRenderSelection');
        const numVertSize        = document.getElementById('numVertSize');
        const numSelSize         = document.getElementById('numSelSize');

        if (this.renderer) {
            if (chkRenderVerts && typeof this.renderer.setAtomsVisible === 'function') {
                this.renderer.setAtomsVisible(chkRenderVerts.checked);
                chkRenderVerts.addEventListener('change', () => {
                    this.renderer.setAtomsVisible(chkRenderVerts.checked);
                });
            }

            if (chkRenderEdges && typeof this.renderer.setBondsVisible === 'function') {
                this.renderer.setBondsVisible(chkRenderEdges.checked);
                chkRenderEdges.addEventListener('change', () => {
                    this.renderer.setBondsVisible(chkRenderEdges.checked);
                });
            }

            if (chkRenderLabelsVis && typeof this.renderer.setLabelsVisible === 'function') {
                this.renderer.setLabelsVisible(chkRenderLabelsVis.checked);
                chkRenderLabelsVis.addEventListener('change', () => {
                    this.renderer.setLabelsVisible(chkRenderLabelsVis.checked);
                });
            }

            if (chkRenderSelection && typeof this.renderer.setSelectionVisible === 'function') {
                this.renderer.setSelectionVisible(chkRenderSelection.checked);
                chkRenderSelection.addEventListener('change', () => {
                    this.renderer.setSelectionVisible(chkRenderSelection.checked);
                });
            }

            if (numVertSize && typeof this.renderer.setNodeScale === 'function') {
                const applyNodeScale = () => {
                    const v = parseFloat(numVertSize.value);
                    if (!Number.isNaN(v) && v > 0) {
                        this.renderer.setNodeScale(v);
                    }
                };
                applyNodeScale();
                numVertSize.addEventListener('change', applyNodeScale);
            }

            if (numSelSize && typeof this.renderer.setSelectionScale === 'function') {
                const applySelScale = () => {
                    const v = parseFloat(numSelSize.value);
                    if (!Number.isNaN(v) && v > 0) {
                        this.renderer.setSelectionScale(v);
                    }
                };
                applySelScale();
                numSelSize.addEventListener('change', applySelScale);
            }
        }
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
                if (this.renderer && typeof this.renderer.resize === 'function') {
                    this.renderer.resize();
                }
            }
        });
    }

    runScript(code) {
        if (!code) return;
        logger.info("Running script...");

        // Simple rule: if script touches mesh or ConstructionBlockTests, run on main thread
        // with explicit access to engine + mesh. We do NOT try to hide errors if the user
        // redeclares 'mesh' locally; that should fail loudly.
        if (code.includes('mesh.') || code.includes('ConstructionBlockTests')) {
            const func = new Function('engine', 'Vec3', 'ConstructionBlockTests', 'window',
                `"use strict";\n${code}`
            );

            func(this.engine, Vec3, ConstructionBlockTests, window);

            if (this.renderer) {
                this.renderer.updateGeometry(this.engine.mesh);
            }
            logger.info("Script executed on Main Thread.");
        } else {
            // Run in Worker (Simulation Script)
            this.engine.runScript(code);
        }
    }

    /**
     * Configure mouse-related camera controls (OrbitControls buttons).
     * LMB is reserved for selection; MMB = dolly, RMB = rotate.
     */
    initMouseControls() {
        if (this.renderer && this.renderer.controls) {
            const controls = this.renderer.controls;
            controls.enableDamping = false;
            controls.dampingFactor = 0.05;

            // Match MolGUI-style controls:
            // LMB: reserved for selection (handled elsewhere), so disable for camera.
            // MMB: zoom/dolly, RMB: rotate.
            controls.mouseButtons = {
                LEFT: null,
                MIDDLE: THREE.MOUSE.DOLLY,
                RIGHT: THREE.MOUSE.ROTATE
            };
        }
    }

    /**
     * Configure global key handlers (e.g. Shift+RMB -> pan for OrbitControls).
     */
    initKeyControls() {
        if (this.renderer && this.renderer.controls) {
            const controls = this.renderer.controls;

            window.addEventListener('keydown', (e) => {
                if (e.key === 'Shift') {
                    controls.mouseButtons.RIGHT = THREE.MOUSE.PAN;
                }
            });

            window.addEventListener('keyup', (e) => {
                if (e.key === 'Shift') {
                    controls.mouseButtons.RIGHT = THREE.MOUSE.ROTATE;
                }
            });
        }
    }

    /**
     * Wire up all DOM widgets (buttons, inputs, dropdowns) to the engine/renderer.
     */
    initWidgetCallbacks() {
        // --- Scripting ---
        const runBtn = document.getElementById('run-script-btn');
        const scriptInput = document.getElementById('script-input');
        const loadScriptBtn = document.getElementById('load-script-btn');
        const fileInput = document.getElementById('script-file-input');

        if (runBtn) { runBtn.addEventListener('click', () => { this.runScript(scriptInput.value); }); }

        if (loadScriptBtn && fileInput) {
            loadScriptBtn.addEventListener('click', () => fileInput.click());
            fileInput.addEventListener('change', (e) => {
                const file = e.target.files[0];
                if (!file) return;
                const reader = new FileReader();
                reader.onload = (ev) => { scriptInput.value = ev.target.result; };
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

        // View Controls
        const btnCamOrtho = document.getElementById('btnCamOrtho');
        const btnCamPersp = document.getElementById('btnCamPersp');
        const btnViewXY = document.getElementById('btnViewXY');
        const btnViewXZ = document.getElementById('btnViewXZ');
        const btnViewYZ = document.getElementById('btnViewYZ');
        
        const chkShowGrid = document.getElementById('chkShowGrid');
        const chkShowAxis = document.getElementById('chkShowAxis');

        if (chkShowGrid) {
            // Apply initial state
            if (this.renderer) this.renderer.setGridVisible(chkShowGrid.checked);
            chkShowGrid.addEventListener('change', () => this.renderer?.setGridVisible(chkShowGrid.checked));
        }
        if (chkShowAxis) {
            // Apply initial state
            if (this.renderer) this.renderer.setAxisVisible(chkShowAxis.checked);
            chkShowAxis.addEventListener('change', () => this.renderer?.setAxisVisible(chkShowAxis.checked));
        }

        if (btnCamOrtho) btnCamOrtho.addEventListener('click', () => {
            if (this.renderer) this.renderer.setCameraMode('ortho');
        });
        if (btnCamPersp) btnCamPersp.addEventListener('click', () => {
            if (this.renderer) this.renderer.setCameraMode('perspective');
        });

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
        if (selTest) {
            this.initTestDropdown(selTest, btnRunTest);
            selTest.addEventListener('change', () => {
                this.onTestSelected(selTest, scriptInput);
            });
        }


        // --- Verbosity ---
        const numVerbosityUI = document.getElementById('numVerbosityUI');
        const numVerbosityCon = document.getElementById('numVerbosityCon');

        if (numVerbosityUI && logger) {
            // Initialize Logger from HTML value
            const initialVal = parseInt(numVerbosityUI.value);
            if (!isNaN(initialVal)) {
                logger.setUIVerbosity(initialVal);
            } else {
                numVerbosityUI.value = logger.uiVerbosity;
            }

            numVerbosityUI.addEventListener('change', () => {
                const val = parseInt(numVerbosityUI.value);
                if (!isNaN(val)) {
                    logger.setUIVerbosity(val);
                    logger.info(`UI Verbosity set to ${val}`);
                }
            });
        }

        if (numVerbosityCon && logger) {
            // Initialize Logger from HTML value
            const initialVal = parseInt(numVerbosityCon.value);
            if (!isNaN(initialVal)) {
                logger.setConsoleVerbosity(initialVal);
            } else {
                numVerbosityCon.value = logger.consoleVerbosity;
            }

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
        // --- Selection Panel ---
        const selBank      = document.getElementById('selSelectionBank');
        const btnSelSave   = document.getElementById('btnSelSave');
        const btnSelDelete = document.getElementById('btnSelDelete');
        const selKind      = document.getElementById('selSelectionKind');
        const chkAdditive  = document.getElementById('chkSelectionAdditive');
        const lblCount     = document.getElementById('lblSelectionCount');
        const txtIndices   = document.getElementById('txtSelectionIndices');
        const btnScanDupVerts   = document.getElementById('btnScanDupVerts');
        const btnScanDupEdges   = document.getElementById('btnScanDupEdges');
        const chkTopoCheckVerts = document.getElementById('chkTopoCheckVerts');
        const chkTopoVertError  = document.getElementById('chkTopoVertError');
        const chkTopoVertSkip   = document.getElementById('chkTopoVertSkip');
        const numRvertCollapse  = document.getElementById('numRvertCollapse');
        const chkTopoCheckEdges = document.getElementById('chkTopoCheckEdges');
        const chkTopoEdgeError  = document.getElementById('chkTopoEdgeError');
        const chkTopoEdgeSkip   = document.getElementById('chkTopoEdgeSkip');

        const mesh = this.engine ? this.engine.mesh : null;

        const updateSelectionUI = () => {
            if (!mesh || !lblCount || !txtIndices) return;
            const arr = mesh.getSelectionArray();
            lblCount.textContent = arr.length.toString();
            txtIndices.value = arr.join(', ');
            if (selKind) {
                selKind.value = mesh.selection.kind;
            }
            // Update selection rendering in SpaceCraftRenderer
            if (this.renderer && typeof this.renderer.updateSelection === 'function') {
                this.renderer.updateSelection(arr);
            }
        };

        if (mesh) {
            updateSelectionUI();

            if (chkTopoCheckVerts)  chkTopoCheckVerts.checked  = !!mesh.bCheckVertExist;
            if (chkTopoVertError)   chkTopoVertError.checked   = !!mesh.bVertExistError;
            if (chkTopoVertSkip)    chkTopoVertSkip.checked    = !!mesh.bVertExistSkip;
            if (numRvertCollapse && typeof mesh.RvertCollapse === 'number') {
                numRvertCollapse.value = mesh.RvertCollapse;
            }
            if (chkTopoCheckEdges)  chkTopoCheckEdges.checked  = !!mesh.bCheckEdgeExist;
            if (chkTopoEdgeError)   chkTopoEdgeError.checked   = !!mesh.bEdgeExistError;
            if (chkTopoEdgeSkip)    chkTopoEdgeSkip.checked    = !!mesh.bEdgeExistSkip;
        }

        if (selKind && mesh) {
            selKind.addEventListener('change', () => {
                mesh.setSelectionKind(selKind.value);
                updateSelectionUI();
            });
        }

        if (selBank && mesh) {
            selBank.addEventListener('change', () => {
                const i = parseInt(selBank.value) || 0;
                mesh.loadSelectionFromBank(i);
                updateSelectionUI();
            });
        }

        if (btnSelSave && mesh && selBank) {
            btnSelSave.addEventListener('click', () => {
                const i = parseInt(selBank.value) || 0;
                mesh.saveSelectionToBank(i);
                logger.info(`Selection saved to bank ${i}.`);
            });
        }

        if (btnSelDelete && mesh && selBank) {
            btnSelDelete.addEventListener('click', () => {
                const i = parseInt(selBank.value) || 0;
                mesh.clearSelectionBank(i);
                if (mesh.currentSelectionBank === i) {
                    mesh.clearSelection();
                }
                updateSelectionUI();
                logger.info(`Selection bank ${i} cleared.`);
            });
        }

        if (txtIndices && mesh) {
            txtIndices.addEventListener('change', () => {
                const text = txtIndices.value;
                const parts = text.split(',');
                const indices = [];
                for (const p of parts) {
                    const v = parseInt(p.trim());
                    if (!Number.isNaN(v) && v >= 0) {
                        indices.push(v);
                    }
                }
                const additive = chkAdditive ? chkAdditive.checked : false;
                const kind = selKind ? selKind.value : 'vert';
                mesh.applySelection(indices, kind, additive);
                updateSelectionUI();
            });
        }

        if (btnScanDupVerts && mesh) {
            btnScanDupVerts.addEventListener('click', () => {
                const count = mesh.scanDuplicateVerts();
                logger.info(`GUI: scanDuplicateVerts finished, count=${count}`);
            });
        }

        if (btnScanDupEdges && mesh) {
            btnScanDupEdges.addEventListener('click', () => {
                const count = mesh.scanDuplicateEdges();
                logger.info(`GUI: scanDuplicateEdges finished, count=${count}`);
            });
        }

        if (mesh && chkTopoCheckVerts) {
            chkTopoCheckVerts.addEventListener('change', () => {
                mesh.bCheckVertExist = chkTopoCheckVerts.checked;
            });
        }
        if (mesh && chkTopoVertError) {
            chkTopoVertError.addEventListener('change', () => {
                mesh.bVertExistError = chkTopoVertError.checked;
            });
        }
        if (mesh && chkTopoVertSkip) {
            chkTopoVertSkip.addEventListener('change', () => {
                mesh.bVertExistSkip = chkTopoVertSkip.checked;
            });
        }
        if (mesh && numRvertCollapse) {
            numRvertCollapse.addEventListener('change', () => {
                const v = parseFloat(numRvertCollapse.value);
                if (!Number.isNaN(v) && v > 0) {
                    mesh.RvertCollapse = v;
                }
            });
        }
        if (mesh && chkTopoCheckEdges) {
            chkTopoCheckEdges.addEventListener('change', () => {
                mesh.bCheckEdgeExist = chkTopoCheckEdges.checked;
            });
        }
        if (mesh && chkTopoEdgeError) {
            chkTopoEdgeError.addEventListener('change', () => {
                mesh.bEdgeExistError = chkTopoEdgeError.checked;
            });
        }
        if (mesh && chkTopoEdgeSkip) {
            chkTopoEdgeSkip.addEventListener('change', () => {
                mesh.bEdgeExistSkip = chkTopoEdgeSkip.checked;
            });
        }

        // --- Mouse box selection (verts/edges), modeled after molgui_web Editor ---
        if (this.renderer && this.renderer.renderer && this.renderer.renderer.domElement) {
            const canvas = this.renderer.renderer.domElement;
            const selectionBox = this.selectionBox;
            let isDragging = false;
            let dragStart = null;

            canvas.addEventListener('pointerdown', (e) => {
                console.log('GUI box-select canvas pointerdown RAW', { button: e.button, x: e.clientX, y: e.clientY, target: e.target && e.target.id });
                if (logger && logger.info) {
                    logger.info(`GUI box-select canvas pointerdown RAW button=${e.button} at (${e.clientX}, ${e.clientY})`);
                }
                // Only left button for selection
                if (e.button !== 0) return;
                isDragging = true;
                dragStart = { x: e.clientX, y: e.clientY };
                if (selectionBox) {
                    selectionBox.style.display = 'none';
                    selectionBox.style.width = '0px';
                    selectionBox.style.height = '0px';
                }
                console.log('GUI box-select pointerdown', dragStart);
                if (logger && logger.info) {
                    logger.info(`Box selection started at (${dragStart.x}, ${dragStart.y})`);
                }
            });

            window.addEventListener('pointermove', (e) => {
                if (!isDragging || !dragStart || !selectionBox) return;

                const dx = e.clientX - dragStart.x;
                const dy = e.clientY - dragStart.y;

                // Only show box when drag distance is significant
                if (Math.abs(dx) > 5 || Math.abs(dy) > 5) {
                    if (selectionBox.style.display !== 'block') {
                        selectionBox.style.display = 'block';
                    }

                    const x = dx < 0 ? e.clientX : dragStart.x;
                    const y = dy < 0 ? e.clientY : dragStart.y;
                    const w = Math.abs(dx);
                    const h = Math.abs(dy);

                    selectionBox.style.transform = `translate(${x}px, ${y}px)`;
                    selectionBox.style.width = `${w}px`;
                    selectionBox.style.height = `${h}px`;

                    // DEBUG: log current rectangle coordinates/size to verify wiring
                    //console.log('GUI box-select move rect', { x, y, w, h });
                    //if (logger && logger.debug) {logger.debug(`Box-select rect x=${x}, y=${y}, w=${w}, h=${h}`);}
                }
            });

            window.addEventListener('pointerup', (e) => {
                console.log('GUI box-select window pointerup RAW', { isDragging, hasDragStart: !!dragStart, button: e.button, x: e.clientX, y: e.clientY });
                if (!isDragging || !dragStart) return;
                isDragging = false;
                if (selectionBox) {
                    selectionBox.style.display = 'none';
                }
                const dragEnd = { x: e.clientX, y: e.clientY };

                const minX = Math.min(dragStart.x, dragEnd.x);
                const maxX = Math.max(dragStart.x, dragEnd.x);
                const minY = Math.min(dragStart.y, dragEnd.y);
                const maxY = Math.max(dragStart.y, dragEnd.y);

                const rect = canvas.getBoundingClientRect();
                const cam = this.renderer.camera;
                const width = rect.width;
                const height = rect.height;

                const meshNow = this.engine ? this.engine.mesh : null;
                if (!meshNow) {
                    console.warn('GUI box-select: engine.mesh is not available at pointerup; skipping selection.');
                    if (logger && logger.warn) {
                        logger.warn('GUI box-select: engine.mesh is not available at pointerup; skipping selection.');
                    }
                    return;
                }

                const kind = selKind ? selKind.value : 'vert';

                // Determine mode: replace / add / subtract (like MolGUI Editor)
                let mode = 'replace';
                if (e.shiftKey || (chkAdditive && chkAdditive.checked)) {
                    mode = 'add';
                } else if (e.ctrlKey || e.altKey) {
                    mode = 'subtract';
                }

                console.log('GUI box-select pointerup', { start: dragStart, end: dragEnd, kind, mode });
                if (logger && logger.info) {
                    logger.info(`Box selection drag from (${dragStart.x}, ${dragStart.y}) to (${dragEnd.x}, ${dragEnd.y}), kind=${kind}, mode=${mode}`);
                }

                const v = new THREE.Vector3();
                const indices = [];

                if (kind === 'vert') {
                    // SDF-based vertex selection using SDfuncs.SDF_ScreenRectWorld and MeshBuilder.selectVertsBySDF.
                    // If SDfuncs or selectVertsBySDF are not wired correctly, this should throw loudly.
                    const sdf = SDfuncs.SDF_ScreenRectWorld(cam, rect, minX, maxX, minY, maxY);

                    if (mode === 'replace') {
                        meshNow.selectVertsBySDF(sdf, 0.0, true);
                    } else if (mode === 'add') {
                        meshNow.selectVertsBySDF(sdf, 0.0, false);
                    } else if (mode === 'subtract') {
                        // Subtract mode: delegate to MeshBuilder helper that uses Selection.subtract.
                        const nRemoved = meshNow.subtractVertsBySDF(sdf, 0.0);
                        // For logging we only care about count, so reflect that in indices.length.
                        indices.length = nRemoved;
                    }
                } 
                /*
                else if (kind === 'vert') {
                    // Fallback: local SDF-based implementation (kept for backup).
                    const width = rect.width;
                    const height = rect.height;
                    const toNDC = (sx, sy) => {
                        const xndc = ((sx - rect.left) / width) * 2.0 - 1.0;
                        const yndc = 1.0 - 2.0 * ((sy - rect.top) / height);
                        return { x: xndc, y: yndc };
                    };

                    const p00 = toNDC(minX, minY);
                    const p11 = toNDC(maxX, maxY);
                    const minXN = Math.min(p00.x, p11.x);
                    const maxXN = Math.max(p00.x, p11.x);
                    const minYN = Math.min(p00.y, p11.y);
                    const maxYN = Math.max(p00.y, p11.y);

                    const tmp = new THREE.Vector3();
                    const sdfScreenRect = (pos) => {
                        tmp.set(pos.x, pos.y, pos.z).project(cam);
                        const x = tmp.x;
                        const y = tmp.y;
                        const dx = (x < minXN) ? (minXN - x) : (x > maxXN ? x - maxXN : 0.0);
                        const dy = (y < minYN) ? (minYN - y) : (y > maxYN ? y - maxYN : 0.0);
                        if (dx === 0.0 && dy === 0.0) return -1.0;
                        return Math.hypot(dx, dy);
                    };

                    for (let i = 0; i < meshNow.verts.length; i++) {
                        const p = meshNow.verts[i].pos;
                        if (sdfScreenRect(p) < 0.0) {
                            indices.push(i);
                        }
                    }
                
                }
                    */ 
                else if (kind === 'edge') {
                    for (let ie = 0; ie < meshNow.edges.length; ie++) {
                        const eEdge = meshNow.edges[ie];
                        const pa = meshNow.verts[eEdge.x].pos;
                        const pb = meshNow.verts[eEdge.y].pos;

                        v.set(pa.x, pa.y, pa.z);
                        v.project(cam);
                        const sax = rect.left + (v.x + 1) * 0.5 * width;
                        const say = rect.top + (1 - (v.y + 1) * 0.5) * height;

                        const vb = new THREE.Vector3(pb.x, pb.y, pb.z);
                        vb.project(cam);
                        const sbx = rect.left + (vb.x + 1) * 0.5 * width;
                        const sby = rect.top + (1 - (vb.y + 1) * 0.5) * height;

                        if (sax >= minX && sax <= maxX && say >= minY && say <= maxY &&
                            sbx >= minX && sbx <= maxX && sby >= minY && sby <= maxY) {
                            indices.push(ie);
                        }
                    }

                    // Apply edge selection based on mode. Vertices are handled directly via SDF above.
                    if (mode === 'replace') {
                        meshNow.applySelection(indices, 'edge', false);
                    } else if (mode === 'add') {
                        meshNow.applySelection(indices, 'edge', true);
                    } else if (mode === 'subtract') {
                        meshNow.subtractSelection(indices);
                    }
                }
                updateSelectionUI();
                logger.info(`Box-selected ${indices.length} ${kind === 'edge' ? 'edges' : 'vertices'} with mode=${mode}.`);
            });
        }
    }

    /**
     * Populate the Tests dropdown from ConstructionBlockTests.tests.
     * This only fills the <select> and hides btnRunTest; it does not attach
     * the on-change handler (see onTestSelected for that logic).
     */
    initTestDropdown(selTest, btnRunTest) {
        // Clear existing options first if any
        selTest.innerHTML = '';
        const defaultOpt = document.createElement('option');
        defaultOpt.value = "";
        defaultOpt.textContent = "-- Select Test --";
        selTest.appendChild(defaultOpt);

        const tests = ConstructionBlockTests.tests;
        logger.info(`Loaded ${Object.keys(tests).length} tests: ${Object.keys(tests).join(', ')}`);

        for (const testName in tests) {
            const opt = document.createElement('option');
            opt.value = testName;
            opt.textContent = testName;
            selTest.appendChild(opt);
        }

        if (btnRunTest) {
            // btnRunTest is redundant; main flow is via script textarea + Run Script.
            btnRunTest.style.display = 'none';
        }
    }

    /**
     * Handle selection of a test in the Tests dropdown.
     * Extracts the body of the chosen ConstructionBlockTests.tests[name] function,
     * normalizes indentation, and writes it into the scripting textarea with a
     * standard header documenting the context.
     */
    onTestSelected(selTest, scriptInput) {
        const testName = selTest.value;
        logger.info(`Selected test: ${testName}`);
        if (!testName) return;

        const tests = ConstructionBlockTests.tests;
        const testFunc = tests[testName];
        if (!testFunc) {
            logger.error(`Test function not found for ${testName}`);
            return;
        }

        try {
            // Extract function body
            let code = testFunc.toString();
            // Remove function signature "function (engine) {" or "(engine) => {"
            const bodyStart = code.indexOf('{');
            const bodyEnd = code.lastIndexOf('}');

            if (bodyStart === -1 || bodyEnd === -1) {
                logger.error(`Could not parse function body for ${testName}`);
                return;
            }

            code = code.substring(bodyStart + 1, bodyEnd);

            // Normalize indentation
            const lines = code.split('\n');
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

            // Trim and clean up, then add standard header
            code = code.trim();
            code = `// Test: ${testName}\n// Context: engine, mesh, Vec3, ConstructionBlockTests\n\n${code}`;

            scriptInput.value = code;
            logger.info(`Populated script area for ${testName}`);
        } catch (e) {
            logger.error(`Error populating script: ${e.message}`);
        }
    }

    initControls() {
        this.initMouseControls();
        this.initKeyControls();
        this.initWidgetCallbacks();
    }
}

