class GUI {
    constructor(io) {
        this.io = io;
        this.init();
    }

    init() {
        // Create Sidebar
        const sidebar = document.createElement('div');
        sidebar.id = 'gui-sidebar';
        document.body.appendChild(sidebar);

        // Create Resizer
        const resizer = document.createElement('div');
        resizer.id = 'gui-resizer';
        document.body.appendChild(resizer);

        // Resizer Logic
        let isResizing = false;
        resizer.addEventListener('mousedown', (e) => {
            isResizing = true;
            document.body.style.cursor = 'col-resize';
            e.preventDefault(); // Prevent selection
        });

        window.addEventListener('mousemove', (e) => {
            if (!isResizing) return;
            const newWidth = e.clientX;
            if (newWidth < 150 || newWidth > 600) return; // Limits

            sidebar.style.width = `${newWidth}px`;
            resizer.style.left = `${newWidth}px`;

            const canvasContainer = document.getElementById('canvas-container');
            if (canvasContainer) {
                const resizerWidth = 5;
                canvasContainer.style.left = `${newWidth + resizerWidth}px`;
                canvasContainer.style.width = `calc(100vw - ${newWidth + resizerWidth}px)`;
            }

            // Trigger Resize for Renderer
            if (window.app) {
                window.app.onWindowResize();
            }
        });

        window.addEventListener('mouseup', () => {
            if (isResizing) {
                isResizing = false;
                document.body.style.cursor = 'default';
            }
        });

        // --- Section: Selection ---
        this.createSection(sidebar, 'Selection', (container) => {
            const row = document.createElement('div');
            row.className = 'gui-row';

            const lbl = document.createElement('span');
            lbl.textContent = 'Count: ';
            lbl.style.fontSize = '0.9em';
            lbl.style.marginRight = '5px';
            row.appendChild(lbl);

            this.lblCount = document.createElement('span');
            this.lblCount.textContent = '0';
            this.lblCount.style.fontWeight = 'bold';
            row.appendChild(this.lblCount);

            container.appendChild(row);

            this.inpSelection = document.createElement('input');
            this.inpSelection.type = 'text';
            this.inpSelection.className = 'gui-input';
            this.inpSelection.placeholder = 'IDs (e.g. 1,5)';
            this.inpSelection.onchange = (e) => this.onSelectionInputChange(e.target.value);
            container.appendChild(this.inpSelection);
        });

        // --- Section: View ---
        this.createSection(sidebar, 'View', (container) => {
            // Zoom Slider
            const zoomRow = document.createElement('div');
            zoomRow.className = 'gui-row';
            zoomRow.style.flexDirection = 'column';
            zoomRow.style.alignItems = 'flex-start';

            const zoomLabel = document.createElement('label');
            zoomLabel.textContent = 'Zoom Level (Log10)';
            zoomLabel.style.fontSize = '0.9em';
            zoomRow.appendChild(zoomLabel);

            const zoomSlider = document.createElement('input');
            zoomSlider.type = 'range';
            zoomSlider.min = '-2.0';
            zoomSlider.max = '3.0';
            zoomSlider.step = '0.1';
            zoomSlider.value = '1.0'; // Default
            zoomSlider.style.width = '100%';
            zoomSlider.oninput = (e) => {
                const val = parseFloat(e.target.value);
                this.setZoom(Math.pow(10, val));
            };
            zoomRow.appendChild(zoomSlider);
            container.appendChild(zoomRow);

            // View Buttons
            const viewRow = document.createElement('div');
            viewRow.className = 'gui-row';
            viewRow.style.marginTop = '10px';

            const views = [
                { name: 'XY', pos: [0, 0, 20], up: [0, 1, 0] },
                { name: 'XZ', pos: [0, 20, 0], up: [0, 0, -1] },
                { name: 'YZ', pos: [20, 0, 0], up: [0, 1, 0] }
            ];

            views.forEach(v => {
                const btn = document.createElement('button');
                btn.textContent = v.name;
                btn.className = 'gui-btn';
                btn.style.marginRight = '2px';
                btn.onclick = () => this.setView(v.pos, v.up);
                viewRow.appendChild(btn);
            });
            container.appendChild(viewRow);

            // Axis Toggle
            const axisRow = document.createElement('div');
            axisRow.className = 'gui-row';
            axisRow.style.marginTop = '10px';

            const axisLabel = document.createElement('label');
            axisLabel.className = 'gui-checkbox-label';
            const axisChk = document.createElement('input');
            axisChk.type = 'checkbox';
            axisChk.checked = false;
            axisChk.onchange = (e) => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.toggleAxes(e.target.checked);
                }
            };
            axisLabel.appendChild(axisChk);
            axisLabel.appendChild(document.createTextNode('Show Axes'));
            axisRow.appendChild(axisLabel);
            container.appendChild(axisRow);

            // Label Mode Dropdown
            const labelRow = document.createElement('div');
            labelRow.className = 'gui-row';
            labelRow.style.marginTop = '10px';

            const labelLbl = document.createElement('span');
            labelLbl.textContent = 'Labels: ';
            labelRow.appendChild(labelLbl);

            const labelSel = document.createElement('select');
            labelSel.className = 'gui-select';
            labelSel.style.flexGrow = '1';

            const modes = [
                { value: 'none', text: 'None' },
                { value: 'id', text: 'Atom ID' },
                { value: 'element', text: 'Element' },
                { value: 'type', text: 'Atom Type' }
            ];

            modes.forEach(m => {
                const opt = document.createElement('option');
                opt.value = m.value;
                opt.textContent = m.text;
                labelSel.appendChild(opt);
            });

            labelSel.onchange = (e) => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.setLabelMode(e.target.value);
                }
            };
            labelRow.appendChild(labelSel);
            container.appendChild(labelRow);

            // Label Color & Size
            const styleRow = document.createElement('div');
            styleRow.className = 'gui-row';
            styleRow.style.marginTop = '5px';

            const colorLbl = document.createElement('span');
            colorLbl.textContent = 'Color: ';
            styleRow.appendChild(colorLbl);

            const colorInput = document.createElement('input');
            colorInput.type = 'color';
            colorInput.value = '#ffffff';
            colorInput.style.flexGrow = '0';
            colorInput.style.width = '40px';
            styleRow.appendChild(colorInput);

            const sizeLbl = document.createElement('span');
            sizeLbl.textContent = ' Size: ';
            sizeLbl.style.marginLeft = '10px';
            styleRow.appendChild(sizeLbl);

            const sizeInput = document.createElement('input');
            sizeInput.type = 'number';
            sizeInput.value = '0.5';
            sizeInput.step = '0.1';
            sizeInput.min = '0.1';
            sizeInput.style.width = '50px';
            styleRow.appendChild(sizeInput);

            const updateLabelStyle = () => {
                if (window.app && window.app.molRenderer) {
                    window.app.molRenderer.setLabelStyle(colorInput.value, sizeInput.value);
                }
            };

            colorInput.oninput = updateLabelStyle;
            sizeInput.oninput = updateLabelStyle;

            container.appendChild(styleRow);
        });

        // --- Section: Gizmo ---
        this.createSection(sidebar, 'Gizmo', (container) => {
            // Enable Checkbox
            const label = document.createElement('label');
            label.className = 'gui-checkbox-label';
            const chk = document.createElement('input');
            chk.type = 'checkbox';
            chk.checked = true; // Default on
            chk.onchange = (e) => {
                if (window.app && window.app.editor) {
                    window.app.editor.toggleGizmo(e.target.checked);
                }
            };
            label.appendChild(chk);
            label.appendChild(document.createTextNode('Enable Gizmo'));
            container.appendChild(label);

            // Lock Selection Checkbox
            const labelLock = document.createElement('label');
            labelLock.className = 'gui-checkbox-label';
            labelLock.style.marginTop = '5px';
            const chkLock = document.createElement('input');
            chkLock.type = 'checkbox';
            chkLock.checked = false;
            chkLock.onchange = (e) => {
                if (window.app && window.app.editor) {
                    window.app.editor.selectionLocked = e.target.checked;
                    window.logger.info(`Selection ${e.target.checked ? 'Locked' : 'Unlocked'}`);
                }
            };
            labelLock.appendChild(chkLock);
            labelLock.appendChild(document.createTextNode('Lock Selection'));
            container.appendChild(labelLock);

            // Modes
            const modes = ['translate', 'rotate', 'scale'];
            const modeContainer = document.createElement('div');
            modeContainer.style.marginTop = '10px';

            modes.forEach(mode => {
                const mLabel = document.createElement('label');
                mLabel.className = 'gui-checkbox-label';
                mLabel.style.marginBottom = '5px';

                const radio = document.createElement('input');
                radio.type = 'radio';
                radio.name = 'gizmo-mode';
                radio.value = mode;
                radio.checked = (mode === 'translate');
                radio.onchange = () => {
                    if (window.app && window.app.editor) {
                        window.app.editor.setGizmoMode(mode);
                    }
                };

                mLabel.appendChild(radio);
                mLabel.appendChild(document.createTextNode(mode.charAt(0).toUpperCase() + mode.slice(1)));
                modeContainer.appendChild(mLabel);
            });
            container.appendChild(modeContainer);
        });

        // --- Section: Structure ---
        this.createSection(sidebar, 'Structure', (container) => {
            // Recalculate Bonds
            const btnRecalc = document.createElement('button');
            btnRecalc.textContent = 'Recalculate Bonds';
            btnRecalc.className = 'gui-btn';
            btnRecalc.onclick = () => {
                if (window.app && window.app.editor) {
                    window.app.editor.recalculateBonds();
                }
            };
            container.appendChild(btnRecalc);

            const hr = document.createElement('hr');
            hr.style.borderColor = '#444';
            hr.style.margin = '10px 0';
            container.appendChild(hr);

            // Add Atom Controls
            const lblAdd = document.createElement('div');
            lblAdd.textContent = 'Add Atom Settings:';
            lblAdd.style.fontSize = '0.9em';
            lblAdd.style.marginBottom = '5px';
            container.appendChild(lblAdd);

            // Element Dropdown
            const rowEl = document.createElement('div');
            rowEl.className = 'gui-row';
            const lblEl = document.createElement('span');
            lblEl.textContent = 'Element: ';
            rowEl.appendChild(lblEl);

            this.selElement = document.createElement('select');
            this.selElement.className = 'gui-select';
            this.selElement.style.flexGrow = '1';
            this.selElement.onchange = (e) => this.onElementChange(e.target.value);
            rowEl.appendChild(this.selElement);
            container.appendChild(rowEl);

            // Atom Type Dropdown
            const rowType = document.createElement('div');
            rowType.className = 'gui-row';
            rowType.style.marginTop = '5px';
            const lblType = document.createElement('span');
            lblType.textContent = 'Type: ';
            rowType.appendChild(lblType);

            this.selAtomType = document.createElement('select');
            this.selAtomType.className = 'gui-select';
            this.selAtomType.style.flexGrow = '1';
            this.selAtomType.onchange = (e) => this.onAtomTypeChange(e.target.value);
            rowType.appendChild(this.selAtomType);
            container.appendChild(rowType);

            // Populate initially (delayed to ensure MMParams loaded)
            setTimeout(() => this.populateStructureControls(), 1000);
        });

        // --- Section: Geometry ---
        this.createSection(sidebar, 'Geometry', (container) => {
            // Load
            const btnLoad = document.createElement('button');
            btnLoad.textContent = 'Load XYZ File...';
            btnLoad.className = 'gui-btn';
            container.appendChild(btnLoad);

            const fileInput = document.createElement('input');
            fileInput.type = 'file';
            fileInput.accept = '.xyz';
            fileInput.style.display = 'none';
            document.body.appendChild(fileInput);

            btnLoad.onclick = () => fileInput.click();
            fileInput.onchange = (e) => {
                if (e.target.files.length > 0) {
                    this.io.loadXYZ(e.target.files[0]);
                    fileInput.value = '';
                }
            };

            // Save
            const btnSave = document.createElement('button');
            btnSave.textContent = 'Save XYZ';
            btnSave.className = 'gui-btn';
            btnSave.style.marginTop = '5px';
            btnSave.onclick = () => this.io.saveFile();
            container.appendChild(btnSave);

            // Clear
            const btnClear = document.createElement('button');
            btnClear.textContent = 'Clear Scene';
            btnClear.className = 'gui-btn';
            btnClear.style.marginTop = '5px';
            btnClear.onclick = () => {
                this.io.system.clear();
                this.io.renderer.update();
                window.logger.info("Scene cleared.");
            };
            container.appendChild(btnClear);

            // Separator
            const hr = document.createElement('hr');
            hr.style.borderColor = '#444';
            hr.style.margin = '10px 0';
            container.appendChild(hr);

            // Manual Edit
            const btnToggle = document.createElement('button');
            btnToggle.textContent = 'Edit XYZ Manually';
            btnToggle.className = 'gui-btn';
            container.appendChild(btnToggle);

            const textArea = document.createElement('textarea');
            textArea.className = 'gui-textarea';
            textArea.style.display = 'none';
            textArea.style.height = '150px';
            textArea.style.width = '100%';
            textArea.style.marginTop = '5px';
            textArea.style.fontSize = '0.8em';
            textArea.style.fontFamily = 'monospace';
            textArea.placeholder = 'Paste XYZ content here...';
            container.appendChild(textArea);

            btnToggle.onclick = () => {
                const isHidden = textArea.style.display === 'none';
                textArea.style.display = isHidden ? 'block' : 'none';
                if (isHidden) {
                    // Populate with current system state if showing
                    textArea.value = this.io.generateXYZ();
                }
            };

            const btnApply = document.createElement('button');
            btnApply.textContent = 'Apply XYZ';
            btnApply.className = 'gui-btn';
            btnApply.style.marginTop = '5px';
            btnApply.onclick = () => {
                if (textArea.value.trim()) {
                    this.io.loadXYZString(textArea.value);
                }
            };
            container.appendChild(btnApply);
        });

        // --- Section: Parameters ---
        this.createSection(sidebar, 'Parameters', (container) => {
            // Element Types
            this.createParamControl(container, 'Element Types', 'common_resources/ElementTypes.dat',
                (content) => {
                    if (window.app && window.app.mmParams) {
                        window.app.mmParams.parseElementTypes(content);
                        window.app.molRenderer.updateStructure();
                    }
                }
            );

            // Atom Types
            this.createParamControl(container, 'Atom Types', 'common_resources/AtomTypes.dat',
                (content) => {
                    if (window.app && window.app.mmParams) {
                        window.app.mmParams.parseAtomTypes(content);
                        // Atom types might not affect rendering directly yet (unless we use them for something)
                        // But good to update.
                    }
                }
            );
        });

        // --- Section: Help ---
        this.createSection(sidebar, 'Help', (container) => {
            const btnHelp = document.createElement('button');
            btnHelp.textContent = 'Show Controls (?)';
            btnHelp.className = 'gui-btn';
            btnHelp.onclick = () => {
                const help = document.getElementById('help-overlay');
                help.style.display = help.style.display === 'none' ? 'block' : 'none';
            };
            container.appendChild(btnHelp);
        });

        // --- Section: System Log ---
        this.createSection(sidebar, 'System Log', (container) => {
            // Verbosity
            const row = document.createElement('div');
            row.className = 'gui-row';

            const lbl = document.createElement('span');
            lbl.textContent = 'Verbosity: ';
            lbl.style.fontSize = '0.9em';
            row.appendChild(lbl);

            const input = document.createElement('input');
            input.type = 'number';
            input.className = 'gui-input';
            input.style.width = '50px';
            input.style.flexGrow = '0';
            input.value = window.VERBOSITY_LEVEL;
            input.min = 0;
            input.max = 4;
            input.onchange = (e) => {
                window.VERBOSITY_LEVEL = parseInt(e.target.value);
                window.logger.info(`Verbosity set to ${window.VERBOSITY_LEVEL}`);
            };
            row.appendChild(input);

            const btnClear = document.createElement('button');
            btnClear.textContent = 'Clear';
            btnClear.className = 'gui-btn';
            btnClear.style.marginLeft = '5px';
            btnClear.style.flexGrow = '0';
            btnClear.onclick = () => {
                if (window.logger) window.logger.clear();
            };
            row.appendChild(btnClear);
            container.appendChild(row);

            // Log Output
            const logOut = document.createElement('div');
            logOut.className = 'gui-log-output';
            container.appendChild(logOut);

            if (window.logger) window.logger.setContainer(logOut);
        });
    }

    createSection(parent, title, contentFn) {
        const section = document.createElement('div');
        section.className = 'gui-section';

        const titleEl = document.createElement('div');
        titleEl.className = 'gui-section-title';
        titleEl.textContent = title;
        section.appendChild(titleEl);

        contentFn(section);
        parent.appendChild(section);
    }

    updateSelectionCount() {
        const count = this.io.system.selection.size;
        const ids = Array.from(this.io.system.selection).sort((a, b) => a - b).join(', ');
        if (this.lblCount) this.lblCount.textContent = count;
        if (this.inpSelection) this.inpSelection.value = ids;
        window.logger.info(`Selection Updated: ${count} atoms`);
    }

    onSelectionInputChange(value) {
        const parts = value.split(',');
        this.io.system.clearSelection();
        let count = 0;
        for (const part of parts) {
            const id = parseInt(part.trim());
            if (!isNaN(id) && id >= 0 && id < this.io.system.nAtoms) {
                this.io.system.select(id, 'add');
                count++;
            }
        }
        this.io.renderer.update();
        if (this.onSelectionChanged) this.onSelectionChanged();
        // Also update editor gizmo if possible
        if (window.app && window.app.editor) {
            window.app.editor.updateGizmo();
        }
        window.logger.info(`Selection set from input: ${count} atoms.`);
    }

    setZoom(zoomValue) {
        if (window.app && window.app.camera) {
            const cam = window.app.camera;
            // For OrthographicCamera, zoom is controlled by 'zoom' property
            cam.zoom = zoomValue;
            cam.updateProjectionMatrix();
        }
    }

    createParamControl(container, label, defaultPath, onApply) {
        const wrapper = document.createElement('div');
        wrapper.style.marginBottom = '10px';
        wrapper.style.borderBottom = '1px solid #444';
        wrapper.style.paddingBottom = '5px';

        const btnToggle = document.createElement('button');
        btnToggle.textContent = `Show/Hide ${label}`;
        btnToggle.className = 'gui-btn';
        wrapper.appendChild(btnToggle);

        const textArea = document.createElement('textarea');
        textArea.className = 'gui-textarea';
        textArea.style.display = 'none';
        textArea.style.height = '100px';
        textArea.style.width = '100%';
        textArea.style.marginTop = '5px';
        textArea.style.fontSize = '0.8em';
        textArea.style.fontFamily = 'monospace';
        textArea.placeholder = `Paste ${label} content here...`;
        wrapper.appendChild(textArea);

        // Load Default (Fetch)
        btnToggle.onclick = async () => {
            const isHidden = textArea.style.display === 'none';
            textArea.style.display = isHidden ? 'block' : 'none';
            if (isHidden && !textArea.value) {
                try {
                    const res = await fetch(defaultPath);
                    if (res.ok) {
                        textArea.value = await res.text();
                    }
                } catch (e) {
                    console.warn("Failed to fetch default params:", e);
                }
            }
        };

        // Load File Button
        const btnLoad = document.createElement('button');
        btnLoad.textContent = 'Load File...';
        btnLoad.className = 'gui-btn';
        btnLoad.style.marginTop = '5px';
        wrapper.appendChild(btnLoad);

        const fileInput = document.createElement('input');
        fileInput.type = 'file';
        fileInput.style.display = 'none';
        wrapper.appendChild(fileInput);

        btnLoad.onclick = () => fileInput.click();
        fileInput.onchange = (e) => {
            if (e.target.files.length > 0) {
                const reader = new FileReader();
                reader.onload = (ev) => {
                    textArea.value = ev.target.result;
                    textArea.style.display = 'block'; // Show it
                };
                reader.readAsText(e.target.files[0]);
                fileInput.value = '';
            }
        };

        // Apply Button
        const btnApply = document.createElement('button');
        btnApply.textContent = 'Apply';
        btnApply.className = 'gui-btn';
        btnApply.style.marginTop = '5px';
        btnApply.onclick = () => {
            if (textArea.value.trim()) {
                onApply(textArea.value);
                window.logger.info(`${label} updated.`);
            }
        };
        wrapper.appendChild(btnApply);

        container.appendChild(wrapper);
    }

    setView(pos, up) {
        if (window.app && window.app.camera && window.app.controls) {
            const cam = window.app.camera;
            cam.position.set(pos[0], pos[1], pos[2]);
            cam.up.set(up[0], up[1], up[2]);
            cam.lookAt(0, 0, 0);
            cam.updateProjectionMatrix();
            window.app.controls.update();
        }
    }

    populateStructureControls() {
        if (!window.app || !window.app.mmParams) return;
        const params = window.app.mmParams;

        // Populate Elements
        this.selElement.innerHTML = '';
        const elements = Object.values(params.elementTypes).sort((a, b) => a.iZ - b.iZ);

        // Add common ones first if needed, or just all
        for (const el of elements) {
            const opt = document.createElement('option');
            opt.value = el.iZ;
            opt.textContent = `${el.name} (${el.iZ})`;
            if (el.name === 'C') opt.selected = true;
            this.selElement.appendChild(opt);
        }

        // Trigger update for Atom Types
        this.onElementChange(this.selElement.value);
    }

    onElementChange(iZ) {
        iZ = parseInt(iZ);
        if (!window.app || !window.app.mmParams) return;
        const params = window.app.mmParams;

        // Update Editor
        if (window.app.editor) {
            window.app.editor.selectedElement = iZ;
        }

        // Filter Atom Types by Element Name
        const el = params.byAtomicNumber[iZ];
        if (!el) return;
        const elName = el.name;

        this.selAtomType.innerHTML = '';
        const types = Object.values(params.atomTypes).filter(t => t.element_name === elName);

        if (types.length === 0) {
            // Fallback if no specific types
            const opt = document.createElement('option');
            opt.value = elName;
            opt.textContent = elName;
            this.selAtomType.appendChild(opt);
        } else {
            for (const t of types) {
                const opt = document.createElement('option');
                opt.value = t.name;
                opt.textContent = t.name;
                this.selAtomType.appendChild(opt);
            }
        }

        // Trigger Atom Type Update
        if (this.selAtomType.options.length > 0) {
            this.onAtomTypeChange(this.selAtomType.options[0].value);
        }
    }

    onAtomTypeChange(typeName) {
        if (window.app && window.app.editor) {
            window.app.editor.selectedAtomType = typeName;
            window.logger.info(`Selected Atom Type: ${typeName}`);
        }
    }
}
