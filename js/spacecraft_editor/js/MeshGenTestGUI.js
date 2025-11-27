"use strict";

class MeshGenTestGUI {
    constructor(engine, renderer, containerId = 'mesh-gen-panel') {
        this.engine = engine;
        this.renderer = renderer;
        this.containerId = containerId;
        console.log("MeshGenTestGUI: Constructor called");
        this.init();
    }

    init() {
        console.log("MeshGenTestGUI: init() called");
        const panel = document.getElementById(this.containerId);
        if (!panel) {
            console.error(`MeshGenTestGUI: Container '${this.containerId}' not found.`);
            return;
        }
        console.log("MeshGenTestGUI: Panel found", panel);

        // Shared argument builder to reduce boilerplate in run() functions
        const buildCommonArgs = (vals) => {
            const n     = { x: vals.nx, y: vals.ny };
            const UVmin = { x: 0, y: 0 };
            const UVmax = { x: 1, y: 1 };
            const mask  = parseInt(vals.dirMask, 2);
            const args = { n, UVmin, UVmax, mask };
            if ('R' in vals) { args.Rs = { x: vals.R, y: vals.R }; }
            if ('thick' in vals || 'offsetX' in vals || 'offsetY' in vals) {
                args.up = {
                    x: vals.offsetX || 0,
                    y: vals.offsetY || 0,
                    z: vals.thick   || 0
                };
            }
            if ('scale' in vals) { args.scale = vals.scale; }
            if ('R_major' in vals && 'R_minor' in vals) { args.RsTorus = { x: vals.R_minor, y: vals.R_major }; }
            return args;
        };

        const square = (s) => {
            const p00 = { x: -s, y: -s, z: 0 };
            const p01 = { x: -s, y: s, z: 0 };
            const p10 = { x: s, y: -s, z: 0 };
            const p11 = { x: s, y: s, z: 1 }; // twisted
            return [p00, p01, p10, p11];
        };


        // Shared default params (can be extended inline per shape using spread)
        const defaultParams = {
            'nx':     { group: 'Grid',     widget: 'int',    value: 3,    range: [2, 100],  step: 1   },
            'ny':     { group: 'Grid',     widget: 'int',    value: 4,    range: [2, 100],  step: 1   },
            'L':      { group: 'Geometry', widget: 'double', value: 10.0,range: [0.1, 50.0], step: 0.1 },
            'dirMask': { group: 'Topology', widget: 'text', value: "1011" }
            // for slab dirMask should be : 10101000111   with offset (0.5,0.5)
        };

        // Tube-specific params (shared between TubeSheet and SlabTube)
        const params_tube = {
            'R':     { group: 'Geometry', widget: 'double', value: 5.0, range: [0.1, 50.0], step: 0.1 },
            'twist': { group: 'Geometry', widget: 'double', value: 0.5, range: [-5.0, 5.0], step: 0.1 },
        };

        // Slab-specific params (shared between SlabTube and QuadSlab)
        const params_slab = {
            'thick':   { group: 'Geometry', widget: 'double', value: 1.0, range: [0.1, 10.0], step: 0.1 },
            'offsetX': { group: 'Geometry', widget: 'double', value: 0.0, range: [-1.0, 1.0], step: 0.1 },
            'offsetY': { group: 'Geometry', widget: 'double', value: 0.0, range: [-1.0, 1.0], step: 0.1 },
        };

        // --- Mesh Function Definitions ---
        const meshFuncs = {
            'TubeSheet': {
                params: { ...defaultParams, ...params_tube, 'offset': { group: 'Geometry', widget: 'double', value: 0.0, range: [-1.0, 1.0], step: 0.1 }, },
                run: (vals) => {
                    const { n, UVmin, UVmax, Rs, mask } = buildCommonArgs(vals);
                    this.engine.mesh.TubeSheet(n, UVmin, UVmax, Rs, vals.L, mask, vals.twist);
                }
            },
            'TubeSheetPost': {
                params: {  ...defaultParams, ...params_tube,
                    'xskip': { group: 'Post', widget: 'int', value: -1, range: [-10, 10], step: 1 },
                    'yskip': { group: 'Post', widget: 'int', value:  2, range: [-10, 10], step: 1 },
                },
                run: (vals) => {
                    const { n, UVmin, UVmax, Rs, mask } = buildCommonArgs(vals);
                    const mesh = this.engine.mesh;
                    const iv0 = mesh.verts.length;
                    //mesh.TubeSheet(n, UVmin, UVmax, Rs, vals.L, mask, vals.twist);
                    mesh.TubeSheet_swapped(n, UVmin, UVmax, Rs, vals.L, mask, vals.twist);
                    const xskip = vals.xskip | 0;
                    const yskip = vals.yskip | 0;
                    // Use matID=1 for post-processed rails so they can be colored differently.
                    mesh.addSkipEdges(n, iv0, xskip, yskip, 1);
                }
            },
            'SlabTube': {
                params: { ...defaultParams, ...params_tube, ...params_slab },
                run: (vals) => {
                    const { n, UVmin, UVmax, Rs, up, mask } = buildCommonArgs(vals);
                    this.engine.mesh.SlabTube(n, UVmin, UVmax, Rs, vals.L, up, mask, vals.twist);
                }
            },
            'QuadSheet': {
                params: { ...defaultParams,  'offset': { group: 'Geometry', widget: 'double', value: 0.0,  range: [-1.0, 1.0], step: 0.1 }, },
                run: (vals) => {
                    const { n, UVmin, UVmax, mask } = buildCommonArgs(vals);
                    const s = vals.L;
                    const [p00, p01, p10, p11] = square(s);
                    this.engine.mesh.QuadSheet(n, UVmin, UVmax, p00, p01, p10, p11, mask);
                }
            },
            'QuadSlab': {
                params: { ...defaultParams, ...params_slab },
                run: (vals) => {
                    const { n, UVmin, UVmax, up, mask } = buildCommonArgs(vals);
                    const [p00, p01, p10, p11] = square(vals.L);
                    this.engine.mesh.QuadSlab(n, UVmin, UVmax, p00, p01, p10, p11, up, mask);
                }
            },
            'TorusSheet': {
                params: { ...defaultParams, ...params_tube,  'thick':{ group: 'Geometry', widget: 'double', value: 2.0, range: [0.0, 20.0], step: 0.1 },},
                run: (vals) => {
                    const { n, UVmin, UVmax, Rs, mask } = buildCommonArgs(vals);
                    const RsTorus = { x: Rs.x, y: Rs.x + vals.thick };
                    this.engine.mesh.TorusSheet(n, UVmin, UVmax, RsTorus, mask, vals.twist);
                }
            }
        };

        // --- GUI Construction ---
        
        // 1. Selector
        const rowSel = GUIutils.box(panel, false);
        rowSel.style.marginBottom = '5px';
        GUIutils.label(rowSel, "Func:");
        
        const funcNames = {};
        for(let k in meshFuncs) funcNames[k] = k;
        
        const selFunc = GUIutils.select(rowSel, funcNames, 'TubeSheet', (val) => {
            rebuildParams(val);
            if (chkAuto.checked) runMeshGen();
        });

        // 2. Controls (Auto, Update)
        const rowCtrl = GUIutils.box(panel, false);
        rowCtrl.style.marginBottom = '5px';
        const { input: chkAuto } = GUIutils.checkBox(rowCtrl, "Auto", true);
        GUIutils.button(rowCtrl, "Update", () => runMeshGen());

        // 3. Parameter Container
        const paramContainer = document.createElement('div');
        panel.appendChild(paramContainer);

        let currentGui = null;
        let currentFunc = null;

        const runMeshGen = () => {
            if (currentGui && currentFunc) {
                const vals = currentGui.getValues();
                this.engine.mesh.clear();
                currentFunc.run(vals);
                if (this.renderer) this.renderer.updateGeometry(this.engine.mesh);
                if (typeof logger !== 'undefined') logger.info(`Generated ${selFunc.value}`);
            }
        };

        const rebuildParams = (funcName) => {
            paramContainer.innerHTML = ''; // Clear
            const funcData = meshFuncs[funcName];
            if (!funcData) return;
            currentFunc = funcData;
            currentGui = GUIutils.create_gui(paramContainer, funcData.params, (name, val) => { if (chkAuto.checked) runMeshGen(); });
        };

        // Initialize
        rebuildParams(selFunc.value);
    }
}

if (typeof window !== 'undefined') {
    window.MeshGenTestGUI = MeshGenTestGUI;
}
