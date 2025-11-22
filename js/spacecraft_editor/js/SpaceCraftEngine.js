class SpaceCraftEngine {
    constructor() {
        // Configuration
        this.maxNodes = 4096; // Texture size N x 1
        this.verbosity = 1;

        // Data Arrays (CPU Side)
        // Nodes: [x, y, z, type]
        this.nodes = new Float32Array(this.maxNodes * 4);
        this.nodeCount = 0;

        // Girders: [nodeA, nodeB, type, 0]
        this.girders = []; // Using JS array for now for ease of push, can optimize to Int32Array later

        // GPU Data Texture
        this.dataTexture = null; // Created by Renderer

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
                    window.logToUI(`ERROR: ${msg.payload} `);
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
        this.nodeCount = 0;
        this.girders = [];
        this.nodes.fill(0);
        window.logToUI("Engine reset.");
    }

    processCommands(cmds) {
        console.log(`Processing ${cmds.length} commands...`);
        for (const cmd of cmds) {
            switch (cmd.method) {
                case 'Node':
                    this.createNode(cmd.args[0]);
                    break;
                case 'Girder':
                    this.createGirder(cmd.args[0], cmd.args[1], cmd.args[2]);
                    break;
                case 'Material':
                    // TODO: Store materials
                    break;
            }
        }

        // Notify Renderer to update
        if (window.renderer) {
            window.renderer.updateGeometry();
        }
    }

    createNode(pos) {
        if (this.nodeCount >= this.maxNodes) {
            console.warn("Max nodes reached!");
            return;
        }
        const i = this.nodeCount * 4;
        this.nodes[i] = pos[0];
        this.nodes[i + 1] = pos[1];
        this.nodes[i + 2] = pos[2];
        this.nodes[i + 3] = 0; // Default type

        if (this.verbosity >= 2) {
            window.logToUI(`Node[${this.nodeCount}]: ${pos[0].toFixed(2)}, ${pos[1].toFixed(2)}, ${pos[2].toFixed(2)} `);
        }

        this.nodeCount++;
    }

    createGirder(n1, n2, typeName) {
        // Simple mapping for now
        let type = 0;
        if (typeName === "Steel") type = 1;

        if (this.verbosity >= 2) {
            window.logToUI(`Girder: ${n1} -> ${n2} (${typeName})`);
        }

        this.girders.push({
            n1: n1,
            n2: n2,
            type: type
        });
    }
}
