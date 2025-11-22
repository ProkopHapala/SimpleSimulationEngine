class SpaceCraftEngine {
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
                    window.logToUI(msg.payload);
                    break;
                case 'ERROR':
                    window.logToUI(`ERROR: ${msg.payload}`);
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
        window.logger.info("Engine reset.");
    }

    processCommands(cmds) {
        console.log(`Processing ${cmds.length} commands...`);

        // 1. Populate Abstract SpaceCraft
        // We need a temporary map to resolve Shadow IDs from worker to real objects
        const idMap = {
            Node: [],
            Girder: [],
            Rope: []
        };

        for (const cmd of cmds) {
            switch (cmd.method) {
                case 'Node': {
                    const n = this.craft.addNode(cmd.args[0]);
                    idMap.Node.push(n); // Shadow ID corresponds to index
                    break;
                }
                case 'Girder': {
                    // Args: [id1, id2, matName]
                    // Worker sends Shadow IDs (indices)
                    const n1 = idMap.Node[cmd.args[0]];
                    const n2 = idMap.Node[cmd.args[1]];
                    if (n1 && n2) {
                        this.craft.addGirder(n1, n2, cmd.args[2]);
                    } else {
                        console.warn("Invalid node IDs in Girder command", cmd.args);
                    }
                    break;
                }
                case 'Material':
                    // TODO: Store materials
                    break;
            }
        }

        // 2. Generate Concrete Mesh
        BuildCraft_truss(this.mesh, this.craft);

        if (this.verbosity >= 1) {
            window.logToUI(`Generated Mesh: ${this.mesh.verts.length / 3} verts, ${this.mesh.edges.length / 4} edges.`);
        }

        // 3. Notify Renderer to update
        if (window.renderer) {
            window.renderer.updateGeometry(this.mesh);
        }
    }
}
