
// SpaceCraftWorker.js

// Internal State for Shadow IDs
const _counters = {
    Node: 0,
    Girder: 0,
    Rope: 0,
    Plate: 0,
    Slider: 0,
    Ring: 0
};

const _commands = [];

// Global UID counter for all objects
let _uid = 0;
const newUid = () => _uid++;

// The API exposed to the user script
const api = {
    Material: function (name, density, Spull, Spush) {
        _commands.push({ method: 'Material', args: [name, density, Spull, Spush] });
    },
    Node: function (pos, size = 1) {
        const uid = newUid(); _counters.Node++;
        _commands.push({ method: 'Node', args: [pos, size], id: uid });
        return uid;
    },
    Girder: function (n1, n2, nseg = 1, matName) {
        const uid = newUid(); _counters.Girder++;
        _commands.push({ method: 'Girder', args: [n1, n2, nseg, matName], id: uid });
        return uid;
    },
    Rope: function (n1, n2, thick, matName) {
        const uid = newUid(); _counters.Rope++;
        _commands.push({ method: 'Rope', args: [n1, n2, thick, matName], id: uid });
        return uid;
    },
    Radiator: function (b1, b2, span1 = [0, 1], span2 = [0, 1], matName, nx = 2, ny = 2, nz = 1, upA = true, upB = true, sideOff = 0, weldDist = 0) {
        const uid = newUid(); _counters.Plate++;
        _commands.push({ method: 'Plate', args: [b1, b2, span1, span2, matName, 'Radiator', nx, ny, nz, upA, upB, sideOff, weldDist], id: uid });
        return uid;
    },
    Shield: function (b1, b2, span1 = [0, 1], span2 = [0, 1], matName, nx = 2, ny = 2, nz = 1, upA = true, upB = true, sideOff = 0, weldDist = 0) {
        const uid = newUid(); _counters.Plate++;
        _commands.push({ method: 'Plate', args: [b1, b2, span1, span2, matName, 'Shield', nx, ny, nz, upA, upB, sideOff, weldDist], id: uid });
        return uid;
    },
    // Args: railUid, slidingUid, calong, matName, slidingVertId, side, methodFlag
    Slider: function (railUid, slidingUid = null, calong = 0.5, matName = null, slidingVertId = -1, side = 0, methodFlag = true) {
        const uid = newUid(); _counters.Slider++;
        _commands.push({ method: 'Slider', args: [railUid, slidingUid, calong, matName, slidingVertId, side, methodFlag], id: uid });
        console.log(`[Worker] api.Slider: rail=${railUid} sliding=${slidingUid} calong=${calong} side=${side}`);
        return uid;
    },
    Ring: function (pos, dir, up = null, R, nseg, wh, matName, st, phase = 0.0) {
        const uid = newUid(); _counters.Ring++;
        _commands.push({ method: 'Ring', args: [pos, dir, up, R, nseg, wh, matName, st, phase], id: uid });
        console.log(`[Worker] api.Ring uid=${uid} phase=${phase}`);
        return uid;
    },
    Ring3P: function (p1, p2, p3, nseg, wh, matName, st, phase = 0.0) {
        const uid = newUid(); _counters.Ring++;
        _commands.push({ method: 'Ring3P', args: [p1, p2, p3, nseg, wh, matName, st, phase], id: uid });
        return uid;
    }
};

// Helper for vectors
const vec = (x, y, z) => [x, y, z];

// Message Handler
onmessage = function (e) {
    const msg = e.data;
    if (msg.type === 'RUN') {
        try {
            // Reset counters
            _uid = 0;
            _counters.Node = 0;
            _counters.Girder = 0;
            _counters.Rope = 0;
            _counters.Plate = 0;
            _counters.Slider = 0;
            _counters.Ring = 0;
            _commands.length = 0;

            // Execute User Code
            // We use a Function constructor to isolate scope somewhat
            const userFunc = new Function('api', 'vec', 'self', 'window', 'fetch',
                '"use strict";' + msg.code
            );

            userFunc(api, vec, undefined, undefined, undefined);

            // Send commands back to main thread
            postMessage({ type: 'CMD_BATCH', cmds: _commands });
            postMessage({ type: 'LOG', payload: `Script executed. Generated ${_counters.Node} nodes, ${_counters.Girder} girders.` });

        } catch (err) {
            postMessage({ type: 'ERROR', payload: err.toString() });
        }
    }
};
