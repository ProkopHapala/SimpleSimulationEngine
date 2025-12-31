import { SpaceCraft } from './SpaceCraft.js';
import { MeshBuilder } from '../../common_js/MeshBuilder.js';
import { logger } from '../../common_js/Logger.js';
import { BuildCraft_blocks_js, BuildCraft_aux_js } from './SpaceCraft2Mesh.js';
import { extendMeshBuilder } from '../../common_js/MeshesUV.js';
import { extendMeshBuilderWithGenerators } from './MeshGenerators.js';

// Ensure MeshBuilder prototype is extended with UV and generator helpers
// when used via ES modules (the old global auto-extend path relied on
// window.MeshBuilder, which we no longer depend on here).
if (!MeshBuilder.__extendedForSpacecraft) {
    extendMeshBuilder(MeshBuilder);
    extendMeshBuilderWithGenerators(MeshBuilder);
    MeshBuilder.__extendedForSpacecraft = true;
}

export class SpaceCraftEngine {
    constructor() {
        // Configuration
        this.maxVerts = 10000;
        // this.verbosity = 2; // Use global window.VERBOSITY_LEVEL

        // Pipeline Objects
        this.craft = new SpaceCraft();
        this.mesh = new MeshBuilder();
        this.auxMesh = new MeshBuilder();
        this.renderer = null; // Linked in main.js
        // Global registry for unambiguous references
        this.uidCounter = 0;
        this.globalRegistry = {};

        // GPU Data Texture (still used for positions)
        this.dataTexture = null;

        // Worker
        this.initWorker();
    }

    initWorker() {
        this.worker = new Worker('js/SpaceCraftWorker.js');

        this.worker.onmessage = (e) => {
            const msg = e.data;
            switch (msg.type) {
                case 'CMD_BATCH':
                    this.processCommands(msg.cmds);
                    break;
                case 'LOG':
                    logger.info(msg.payload);
                    break;
                case 'ERROR':
                    logger.error(msg.payload);
                    break;
            }
        };
    }

    runScript(code) {
        // Reset state before running
        this.reset();
        this.worker.postMessage({ type: 'RUN', code: code });
    }

    reset() {
        this.craft.clear();
        this.mesh.clear();
        this.auxMesh.clear();
        this.uidCounter = 0;
        this.globalRegistry = {};
        logger.info("Engine reset.");
    }

    rebuildMesh() {
        if (!this.mesh || !this.craft) return;
        BuildCraft_blocks_js(this.mesh, this.craft);
        this.updateAux();
        if (this.renderer) {
            this.renderer.updateGeometry(this.mesh);
        }
    }

    updateAux() {
        if (!this.auxMesh || !this.craft || !this.mesh) return;
        BuildCraft_aux_js(this.auxMesh, this.craft, this.mesh);
        if (this.renderer) {
            this.renderer.updateAuxGeometry(this.auxMesh);
        }
    }

    processCommands(cmds) {
        logger.info(`Processing ${cmds.length} commands...`);

        // 1. Populate Abstract SpaceCraft
        // We need a temporary map to resolve Shadow IDs from worker to real objects
        const idMap = { Node: [], Girder: [], Rope: [], Plate: [], Slider: [], Ring: [] };
        const seq = { Node: 0, Girder: 0, Rope: 0, Plate: 0, Slider: 0, Ring: 0 };
        const logV = (lvl, msg) => { try { if (typeof logger !== 'undefined' && logger.uiVerbosity >= lvl) logger.info(msg); } catch (e) {} };
        const addToRegistry = (type, obj, uid) => {
            obj.uid = uid;
            this.globalRegistry[uid] = { type, obj };
        };
        const resolveRef = (uid) => {
            if (uid === null || uid === undefined) return null;
            return this.globalRegistry[uid]?.obj || idMap.Node[uid] || idMap.Girder[uid] || idMap.Ring[uid] || idMap.Rope[uid] || idMap.Plate[uid] || idMap.Slider[uid] || null;
        };

        for (const cmd of cmds) {
            switch (cmd.method) {
                case 'Node': {
                    // Args: [pos, size]
                    const idx = (cmd.id !== undefined) ? cmd.id : seq.Node;
                    const n = this.craft.addNode(cmd.args[0], cmd.args[1]);
                    idMap.Node[idx] = n;
                    addToRegistry('Node', n, idx);
                    seq.Node++;
                    break;
                }
                case 'Girder': {
                    // Args: [id1, id2, nseg, matName]
                    const n1 = resolveRef(cmd.args[0]);
                    const n2 = resolveRef(cmd.args[1]);
                    if (n1 && n2) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Girder;
                        const g = this.craft.addGirder(n1, n2, cmd.args[2], cmd.args[3]);
                        idMap.Girder[idx] = g;
                        addToRegistry('Girder', g, idx);
                        seq.Girder++;
                    } else {
                        logV(1, `[Engine] Girder skipped: missing nodes ${cmd.args[0]}, ${cmd.args[1]}`);
                    }
                    break;
                }
                case 'Rope': {
                    const n1 = resolveRef(cmd.args[0]);
                    const n2 = resolveRef(cmd.args[1]);
                    if (n1 && n2) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Rope;
                        const r = this.craft.addRope(n1, n2, cmd.args[2], cmd.args[3]);
                        idMap.Rope[idx] = r;
                        addToRegistry('Rope', r, idx);
                        seq.Rope++;
                    } else {
                        logV(1, `[Engine] Rope skipped: missing nodes ${cmd.args[0]}, ${cmd.args[1]}`);
                    }
                    break;
                }
                case 'Plate': {
                    const b1 = resolveRef(cmd.args[0]);
                    const b2 = resolveRef(cmd.args[1]);
                    if (b1 && b2) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Plate;
                        const plate = this.craft.addPlate(
                            b1, b2,
                            cmd.args[2], cmd.args[3],             // spans
                            cmd.args[4], cmd.args[5],             // type/kind
                            cmd.args[6], cmd.args[7], cmd.args[8], // nx,ny,nz
                            cmd.args[9], cmd.args[10],            // upA, upB
                            cmd.args[11], cmd.args[12]            // sideOffset, weldDist
                        );
                        idMap.Plate[idx] = plate;
                        addToRegistry('Plate', plate, idx);
                        seq.Plate++;
                    } else {
                        logV(1, `[Engine] Plate skipped: missing bounds ${cmd.args[0]}, ${cmd.args[1]}`);
                    }
                    break;
                }
                case 'Slider': {
                    const railUid = cmd.args[0];
                    const slidingUid = cmd.args[1]; // sliding component comes second now
                    const rail = resolveRef(railUid);
                    const sliding = resolveRef(slidingUid);
                    
                    logger.info(`[Engine] Slider processing: railUid=${railUid} railFound=${!!rail} (${rail?.constructor?.name}) slidingUid=${slidingUid} slidingFound=${!!sliding} (${sliding?.constructor?.name})`);
                    
                    if (rail) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Slider;
                        const s = this.craft.addSlider(rail, sliding, cmd.args[2], cmd.args[3], cmd.args[4], cmd.args[5], cmd.args[6]);
                        idMap.Slider[idx] = s;
                        addToRegistry('Slider', s, idx);
                        seq.Slider++;
                    } else {
                        logger.warn(`[Engine] Slider skipped: missing rail uid=${railUid}. Available Rings: ${Object.keys(idMap.Ring)}, Girders: ${Object.keys(idMap.Girder)}`);
                    }
                    break;
                }
                case 'Ring': {
                    const pos = cmd.args[0];
                    const dir = cmd.args[1];
                    const up = cmd.args[2];
                    const idx = (cmd.id !== undefined) ? cmd.id : seq.Ring;
                    const R = cmd.args[3];
                    const nseg = cmd.args[4] || 8;
                    const wh = cmd.args[5] || 0.1;
                    const matName = cmd.args[6];
                    const st = cmd.args[7];
                    const ring = this.craft.addRing(pos, dir, up, R, nseg, wh, matName, st);
                    idMap.Ring[idx] = ring;
                    addToRegistry('Ring', ring, idx);
                    seq.Ring++;
                    break;
                }
                case 'Ring3P': {
                    const p1 = cmd.args[0], p2 = cmd.args[1], p3 = cmd.args[2];
                    const nseg = cmd.args[3] || 16, wh = cmd.args[4], matName = cmd.args[5], st = cmd.args[6];
                    const ring = this.craft.addRing3P(p1, p2, p3, nseg, wh, matName, st);
                    const idx = (cmd.id !== undefined) ? cmd.id : seq.Ring;
                    idMap.Ring[idx] = ring;
                    addToRegistry('Ring', ring, idx);
                    seq.Ring++;
                    break;
                }
                case 'Material':
                    // TODO: Store materials
                    break;
            }
        }

        logV(1, `[Engine] Craft counts: nodes=${this.craft.nodes.length} girders=${this.craft.girders.length} ropes=${this.craft.ropes.length} plates=${this.craft.plates.length} rings=${this.craft.rings.length} sliders=${this.craft.sliders.length}`);

        // 2. Generate Concrete Mesh (BuildCraft blocks already computes slider paths)
        BuildCraft_blocks_js(this.mesh, this.craft);
        this.updateAux();

        logger.info(`Generated Mesh: ${this.mesh.verts.length} verts, ${this.mesh.edges.length} edges.`);

        // 3. Notify Renderer to update
        if (this.renderer) {
            this.renderer.updateGeometry(this.mesh);
        }

        // 4. Update GUI components list
        if (window.gui) {
            window.gui.updateSliderList();
        }
    }
}
