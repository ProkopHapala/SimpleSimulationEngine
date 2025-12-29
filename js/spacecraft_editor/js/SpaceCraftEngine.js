import { SpaceCraft } from './SpaceCraft.js';
import { MeshBuilder } from '../../common_js/MeshBuilder.js';
import { logger } from '../../common_js/Logger.js';
import { BuildCraft_blocks_js } from './SpaceCraft2Mesh.js';
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
        logger.info("Engine reset.");
    }

    processCommands(cmds) {
        logger.info(`Processing ${cmds.length} commands...`);

        // 1. Populate Abstract SpaceCraft
        // We need a temporary map to resolve Shadow IDs from worker to real objects
        const idMap = { Node: [], Girder: [], Rope: [], Plate: [], Slider: [], Ring: [] };
        const seq = { Node: 0, Girder: 0, Rope: 0, Plate: 0, Slider: 0, Ring: 0 };
        const logV = (lvl, msg) => { try { if (typeof logger !== 'undefined' && logger.uiVerbosity >= lvl) logger.info(msg); } catch (e) {} };

        for (const cmd of cmds) {
            switch (cmd.method) {
                case 'Node': {
                    // Args: [pos, size]
                    const idx = (cmd.id !== undefined) ? cmd.id : seq.Node;
                    const n = this.craft.addNode(cmd.args[0], cmd.args[1]);
                    idMap.Node[idx] = n;
                    seq.Node++;
                    break;
                }
                case 'Girder': {
                    // Args: [id1, id2, nseg, matName]
                    const n1 = idMap.Node[cmd.args[0]];
                    const n2 = idMap.Node[cmd.args[1]];
                    if (n1 && n2) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Girder;
                        const g = this.craft.addGirder(n1, n2, cmd.args[2], cmd.args[3]);
                        idMap.Girder[idx] = g;
                        seq.Girder++;
                    } else {
                        logV(1, `[Engine] Girder skipped: missing nodes ${cmd.args[0]}, ${cmd.args[1]}`);
                    }
                    break;
                }
                case 'Rope': {
                    const n1 = idMap.Node[cmd.args[0]];
                    const n2 = idMap.Node[cmd.args[1]];
                    if (n1 && n2) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Rope;
                        const r = this.craft.addRope(n1, n2, cmd.args[2], cmd.args[3]);
                        idMap.Rope[idx] = r;
                        seq.Rope++;
                    } else {
                        logV(1, `[Engine] Rope skipped: missing nodes ${cmd.args[0]}, ${cmd.args[1]}`);
                    }
                    break;
                }
                case 'Plate': {
                    const b1 = idMap.Girder[cmd.args[0]] || idMap.Rope[cmd.args[0]];
                    const b2 = idMap.Girder[cmd.args[1]] || idMap.Rope[cmd.args[1]];
                    if (b1 && b2) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Plate;
                        const plate = this.craft.addPlate(b1, b2, cmd.args[2], cmd.args[3], cmd.args[4], cmd.args[5]);
                        idMap.Plate[idx] = plate;
                        seq.Plate++;
                    } else {
                        logV(1, `[Engine] Plate skipped: missing bounds ${cmd.args[0]}, ${cmd.args[1]}`);
                    }
                    break;
                }
                case 'Slider': {
                    const boundTo = idMap.Girder[cmd.args[0]] || idMap.Ring[cmd.args[0]] || idMap.Rope[cmd.args[0]];
                    if (boundTo) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Slider;
                        this.craft.addSlider(boundTo, cmd.args[1], cmd.args[2]);
                        seq.Slider++;
                    } else {
                        logV(1, `[Engine] Slider skipped: missing bound ${cmd.args[0]}`);
                    }
                    break;
                }
                case 'Ring': {
                    const idx = (cmd.id !== undefined) ? cmd.id : seq.Ring;
                    const r = this.craft.addRing(cmd.args[0], cmd.args[1], cmd.args[2], cmd.args[3]);
                    idMap.Ring[idx] = r;
                    seq.Ring++;
                    break;
                }
                case 'Material':
                    // TODO: Store materials
                    break;
            }
        }

        logV(1, `[Engine] Craft counts: nodes=${this.craft.nodes.length} girders=${this.craft.girders.length} ropes=${this.craft.ropes.length} plates=${this.craft.plates.length} rings=${this.craft.rings.length} sliders=${this.craft.sliders.length}`);

        // 2. Generate Concrete Mesh
        BuildCraft_blocks_js(this.mesh, this.craft);

        logger.info(`Generated Mesh: ${this.mesh.verts.length / 3} verts, ${this.mesh.edges.length / 4} edges.`);

        // 3. Notify Renderer to update
        if (window.renderer) {
            window.renderer.updateGeometry(this.mesh);
        }
    }
}
