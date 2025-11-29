"use strict";

import { GUIutils } from '../../common_js/GUIutils.js';
import { logger } from '../../common_js/Logger.js';

export class MeshGenTestGUI {
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
            'nx':     { group: 'Grid',     widget: 'int',    value: 8,    range: [2, 100],  step: 1   },
            'ny':     { group: 'Grid',     widget: 'int',    value: 8,    range: [2, 100],  step: 1   },
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

        const params_angled_ropes = {
            'angle':   { group: 'Geometry', widget: 'double', value: 20.0, range: [0.0, 89.0],  step: 1.0 },
            'nx2':    { group: 'Grid',      widget: 'int',    value: 15,   range: [2, 100],    step: 1   },
            'ropeN1':  { group: 'Grid',     widget: 'int',    value: 10,   range: [2, 100],    step: 1   },
            'ropeN2':  { group: 'Grid',     widget: 'int',    value: 10,   range: [2, 100],    step: 1   },
            'ropeL1':  { group: 'Geometry', widget: 'double', value: 10.0, range: [0.1, 100.0],step: 0.1 },
            'ropeL2':  { group: 'Geometry', widget: 'double', value: 10.0, range: [0.1, 100.0],step: 0.1 },
            'offsetX': { group: 'Geometry', widget: 'double', value: 5.0,  range: [0.0, 100.0],step: 0.1 },
            'inset':   { group: 'Geometry', widget: 'double', value: 0.5,  range: [-10.0,10.0], step: 0.05 },
            'attachR': { group: 'Topology', widget: 'double', value: 1.0,  range: [0.01, 10.0], step: 0.05 },
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
            'QuadParametric': {
                params: { ...defaultParams },
                run: (vals) => {
                    const mesh = this.engine.mesh;
                    const s = vals.L;
                    const [p00, p01, p10, p11] = square(s);
                    const nTop = vals.nx | 0;
                    const nBottom = vals.nx2 | 0;
                    const nRows = vals.ny | 0;
                    mesh.ParametricQuadPatch(nTop, nBottom, nRows, p00, p01, p10, p11);
                }
            },
            'QuadParametricRopes': {
                params: { ...defaultParams, ...params_angled_ropes },
                run: (vals) => {
                    const mesh = this.engine.mesh;
                    const angleRad = vals.angle * Math.PI / 180.0;
                    // Symmetric V around Y axis: left rope rotates -angle, right rope +angle from +Y
                    const dirLeft  = { x: -Math.sin(angleRad), y: Math.cos(angleRad), z: 0 };
                    const dirRight = { x:  Math.sin(angleRad), y: Math.cos(angleRad), z: 0 };

                    const ropeN1 = vals.ropeN1 | 0;
                    const ropeN2 = vals.ropeN2 | 0;
                    const ropeL1 = vals.ropeL1;
                    const ropeL2 = vals.ropeL2;
                    const offsetX = vals.offsetX;
                    const inset = vals.inset;   // can be negative
                    const attachR = vals.attachR;
                    const attachR2 = attachR * attachR;

                    // Rope roots: symmetric around X=0
                    const A0 = { x: -offsetX, y: 0, z: 0 };
                    const B0 = { x:  offsetX, y: 0, z: 0 };
                    const A1 = { x: A0.x + dirLeft.x  * ropeL1, y: A0.y + dirLeft.y  * ropeL1, z: 0 };
                    const B1 = { x: B0.x + dirRight.x * ropeL2, y: B0.y + dirRight.y * ropeL2, z: 0 };

                    // Quad corners inside between ropes, inset towards X=0
                    const p00 = { x: A0.x + inset, y: A0.y, z: 0 };
                    const p01 = { x: A1.x + inset, y: A1.y, z: 0 };
                    const p10 = { x: B0.x - inset, y: B0.y, z: 0 };
                    const p11 = { x: B1.x - inset, y: B1.y, z: 0 };

                    const nTop = vals.nx | 0;
                    const nBottom = vals.nx2 | 0;
                    const nRows = vals.ny | 0;

                    const ivPatch0 = mesh.verts.length;
                    mesh.ParametricQuadPatch(nTop, nBottom, nRows, p00, p01, p10, p11);

                    const rope1 = mesh.rope(A0, dirLeft,  ropeL1, ropeN1);
                    const rope2 = mesh.rope(B0, dirRight, ropeL2, ropeN2);

                    const ivPatch1 = mesh.verts.length;
                    for (const rv of [...rope1, ...rope2]) {
                        const rp = mesh.verts[rv];
                        for (let iv = ivPatch0; iv < ivPatch1; iv++) {
                            const pv = mesh.verts[iv];
                            const dx = rp.x - pv.x;
                            const dy = rp.y - pv.y;
                            const dz = rp.z - pv.z;
                            const d2 = dx * dx + dy * dy + dz * dz;
                            if (d2 <= attachR2) {
                                mesh.edge(rv, iv);
                            }
                        }
                    }
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
