import { SpaceCraft } from './SpaceCraft.js';
import { MeshBuilder } from '../../common_js/MeshBuilder.js';
import { logger } from '../../common_js/Logger.js';
import { BuildCraft_blocks_js, BuildCraft_aux_js } from './SpaceCraft2Mesh.js';
import { extendMeshBuilder } from '../../common_js/MeshesUV.js';
//import { extendMeshBuilderWithGenerators } from './MeshGenerators.js';

// Ensure MeshBuilder prototype is extended with UV and generator helpers
// when used via ES modules (the old global auto-extend path relied on
// window.MeshBuilder, which we no longer depend on here).
// if (!MeshBuilder.__extendedForSpacecraft) {
//     extendMeshBuilder(MeshBuilder);
//     extendMeshBuilderWithGenerators(MeshBuilder);
//     MeshBuilder.__extendedForSpacecraft = true;
// }

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
        const idMap = { Node: [], Girder: [], Rope: [], Plate: [], Slider: [], Ring: [], Path: [] };
        const seq = { Node: 0, Girder: 0, Rope: 0, Plate: 0, Slider: 0, Ring: 0, Path: 0 };
        const logV = (lvl, msg) => { try { if (typeof logger !== 'undefined' && logger.uiVerbosity >= lvl) logger.info(msg); } catch (e) {} };
        const addToRegistry = (type, obj, uid) => {
            obj.uid = uid;
            this.globalRegistry[uid] = { type, obj };
        };
        const resolveRef = (uid) => {
            if (uid === null || uid === undefined) return null;
            return this.globalRegistry[uid]?.obj || idMap.Node[uid] || idMap.Girder[uid] || idMap.Ring[uid] || idMap.Rope[uid] || idMap.Plate[uid] || idMap.Slider[uid] || idMap.Path[uid] || null;
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
                    const pathUid = cmd.args[0];
                    const slidingUid = cmd.args[1];
                    const path = resolveRef(pathUid);
                    const sliding = resolveRef(slidingUid);
                    
                    logger.info(`[Engine] Slider processing: pathUid=${pathUid} pathFound=${!!path} slidingUid=${slidingUid} slidingFound=${!!sliding}`);
                    
                    if (path) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Slider;
                        const s = this.craft.addSlider(path, sliding, cmd.args[2], cmd.args[3], cmd.args[4]);
                        idMap.Slider[idx] = s;
                        addToRegistry('Slider', s, idx);
                        seq.Slider++;
                    } else {
                        logger.warn(`[Engine] Slider skipped: missing path uid=${pathUid}`);
                    }
                    break;
                }
                case 'Path': {
                    const railUid = cmd.args[0];
                    const rail = resolveRef(railUid);
                    if (rail) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Path;
                        const p = this.craft.addPath(rail, cmd.args[1], cmd.args[2], cmd.args[3]);
                        idMap.Path[idx] = p;
                        addToRegistry('Path', p, idx);
                        seq.Path++;
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
                    const phase = cmd.args[8] || 0.0;
                    const ring = this.craft.addRing(pos, dir, up, R, nseg, wh, matName, st, phase);
                    idMap.Ring[idx] = ring;
                    addToRegistry('Ring', ring, idx);
                    seq.Ring++;
                    break;
                }
                case 'Ring3P': {
                    const p1 = cmd.args[0], p2 = cmd.args[1], p3 = cmd.args[2];
                    const nseg = cmd.args[3] || 16, wh = cmd.args[4], matName = cmd.args[5], st = cmd.args[6], phase = cmd.args[7] || 0.0;
                    const ring = this.craft.addRing3P(p1, p2, p3, nseg, wh, matName, st, phase);
                    const idx = (cmd.id !== undefined) ? cmd.id : seq.Ring;
                    idMap.Ring[idx] = ring;
                    addToRegistry('Ring', ring, idx);
                    seq.Ring++;
                    break;
                }
                case 'GetGirderCircleIntersection': {
                    const girderUid = cmd.args[0];
                    const center = new Vec3().setV(cmd.args[1]);
                    const axis = new Vec3().setV(cmd.args[2]);
                    const radius = cmd.args[3];
                    const girder = resolveRef(girderUid);
                    if (girder) {
                        const p1 = new Vec3().setV(girder.nodeA.pos);
                        const p2 = new Vec3().setV(girder.nodeB.pos);
                        const hits = Vec3.lineCylinderIntersect(p1, p2, center, axis, radius);
                        logger.info(`[Engine] GetGirderCircleIntersection girder=${girderUid} hits=${hits.length}`);
                        hits.forEach((h, i) => {
                            logger.info(`  Hit ${i}: ${h.toString()}`);
                            this.mesh.pointCross(h);
                        });
                    }
                    break;
                }
                case 'RingAttached': {
                    const girderConfigs = cmd.args[0]; // Array of {uid, t, side, sideDir}
                    const params = cmd.args[1];       // {nseg, wh, matName, st, phase, axis}
                    
                    const points = [];
                    const sliders = [];

                    for (const cfg of girderConfigs) {
                        const girder = resolveRef(cfg.uid);
                        if (!girder) continue;

                        const pA = new Vec3().setV(girder.nodeA.pos);
                        const pB = new Vec3().setV(girder.nodeB.pos);
                        const sA = girder.nodeA.size * 0.5;
                        const sB = girder.nodeB.size * 0.5;
                        const t = cfg.t ?? 0.5;
                        let side = cfg.side ?? 0;

                        // Calculate side offset logic similar to Radiator
                        const dir = new Vec3().setSub(pB, pA);
                        const len = dir.normalize();
                        const upHint = params.axis ? new Vec3().setV(params.axis) : new Vec3(0, 1, 0);
                        let sideVec = new Vec3().setCross(dir, upHint);
                        if (sideVec.norm() < 1e-6) sideVec = dir.getSomePerp();
                        else sideVec.normalize();
                        
                        let upVec = new Vec3().setCross(sideVec, dir);
                        upVec.normalize();

                        const cornerDirs = [ {a:1, b:1}, {a:1, b:-1}, {a:-1, b:1}, {a:-1, b:-1} ];

                        // Auto-pick side using either explicit sideDir or fallback to params.axis
                        const sideDirInput = cfg.sideDir || params.axis;
                        if (sideDirInput) {
                            const targetDir = new Vec3().setV(sideDirInput);
                            const lt = targetDir.normalize();
                            if (lt > 1e-9) {
                                let maxDot = -Infinity;
                                for (let i = 0; i < 4; i++) {
                                    const cd = cornerDirs[i];
                                    const off = new Vec3().setV(upVec).mulScalar(cd.b).addMul(sideVec, cd.a);
                                    const lo = off.normalize();
                                    if (lo < 1e-9) continue;
                                    const d = off.dot(targetDir);
                                    if (d > maxDot) {
                                        maxDot = d;
                                        side = i;
                                    }
                                }
                            }
                        }

                        const cd = cornerDirs[side % 4];
                        const offMag = (1-t)*sA + t*sB;
                        const offset = new Vec3().setV(upVec).mulScalar(cd.b * offMag).addMul(sideVec, cd.a * offMag);
                        
                        const pAttach = new Vec3().setAddMul(pA, new Vec3().setSub(pB, pA), t).add(offset);
                        points.push(pAttach);
                        sliders.push({ girder, t, side });
                    }

                    if (points.length >= 3) {
                        const ring = this.craft.addRing3P(points[0], points[1], points[2], params.nseg, params.wh, params.matName, params.st, params.phase);
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Ring;
                        idMap.Ring[idx] = ring;
                        addToRegistry('Ring', ring, idx);
                        seq.Ring++;

                        // 3. SHARED PATH: Create ONE independent Path object for this ring.
                        // If user passed params.ringSide (>=0) we honor it; otherwise auto-pick using inverted axis sorting:
                        // if avg attachment is "above" ring center along ring axis, choose the OPPOSITE edge (bottom) to face girders.
                        const center = new Vec3().setV(ring.pos);
                        let ringAxis = new Vec3().setV(ring.dir);
                        if (ringAxis.normalize() < 1e-9) {
                            ringAxis = params.axis ? new Vec3().setV(params.axis).normalize() : new Vec3(0, 0, 1);
                        }
                        let avgAttach = new Vec3();
                        for (const p of points) avgAttach.add(p);
                        avgAttach.mulScalar(1.0 / points.length);
                        const distAlongAxis = new Vec3().setSub(avgAttach, center).dot(ringAxis);
                        let ringPathSide = (distAlongAxis > 0) ? 3 : 2; // inverted selection relative to attachment half-space
                        if (typeof params.ringSide === 'number' && params.ringSide >= 0) {
                            ringPathSide = params.ringSide;
                        }
                        const ringPath = this.craft.addPath(ring, ringPathSide, true, 0.35);

                        // 4. Attach Sliders.
                        // IMPORTANT: Each slider's vertex is created as a NEW vertex in the mesh
                        // located exactly at the interpolation point (points[i]) used for ring creation.
                        // This ensures the wheel fits perfectly without strain.
                        for (let i = 0; i < points.length; i++) {
                            const sData = sliders[i];
                            // Create slider attaching the shared 'ringPath' to the respective 'girder'
                            const s = this.craft.addSlider(ringPath, sData.girder, 0, params.matName, -1);
                            
                            // This pAttach property triggers NEW vertex creation at exact coordinates in BuildCraft_blocks_js
                            s.pAttach = points[i]; 
                            // The weldDist controls snapping this new vertex to the girder's mesh vertices
                            s.weldDist = params.weldDist || 0.5;
                        }
                    }
                    break;
                }
                case 'Material':
                    // TODO: Store materials
                    break;
                case 'Slider': {
                    const pathUid = cmd.args[0];
                    const slidingUid = cmd.args[1];
                    const path = resolveRef(pathUid);
                    const sliding = resolveRef(slidingUid);
                    
                    if (path) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Slider;
                        const s = this.craft.addSlider(path, sliding, cmd.args[2], cmd.args[3], cmd.args[4]);
                        idMap.Slider[idx] = s;
                        addToRegistry('Slider', s, idx);
                        seq.Slider++;
                    } else {
                        logger.warn(`[Engine] Slider skipped: missing path uid=${pathUid}`);
                    }
                    break;
                }
                case 'Path': {
                    const railUid = cmd.args[0];
                    const rail = resolveRef(railUid);
                    if (rail) {
                        const idx = (cmd.id !== undefined) ? cmd.id : seq.Path;
                        const p = this.craft.addPath(rail, cmd.args[1], cmd.args[2], cmd.args[3]);
                        idMap.Path[idx] = p;
                        addToRegistry('Path', p, idx);
                        seq.Path++;
                    }
                    break;
                }
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
