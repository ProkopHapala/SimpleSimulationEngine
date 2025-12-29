
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

// The API exposed to the user script
const api = {
    Material: function (name, density, Spull, Spush) {
        _commands.push({ method: 'Material', args: [name, density, Spull, Spush] });
    },
    Node: function (pos, size = 1) {
        const id = _counters.Node++;
        _commands.push({ method: 'Node', args: [pos, size], id });
        return id;
    },
    Girder: function (n1, n2, nseg = 1, matName) {
        const id = _counters.Girder++;
        _commands.push({ method: 'Girder', args: [n1, n2, nseg, matName], id });
        return id;
    },
    Rope: function (n1, n2, thick, matName) {
        const id = _counters.Rope++;
        _commands.push({ method: 'Rope', args: [n1, n2, thick, matName], id });
        return id;
    },
    Radiator: function (b1, b2, span1 = [0, 1], span2 = [0, 1], matName, nx = 2, ny = 2, nz = 1, upA = true, upB = true, sideOff = 0, weldDist = 0) {
        const id = _counters.Plate++;
        _commands.push({ method: 'Plate', args: [b1, b2, span1, span2, matName, 'Radiator', nx, ny, nz, upA, upB, sideOff, weldDist], id });
        return id;
    },
    Shield: function (b1, b2, span1 = [0, 1], span2 = [0, 1], matName, nx = 2, ny = 2, nz = 1, upA = true, upB = true, sideOff = 0, weldDist = 0) {
        const id = _counters.Plate++;
        _commands.push({ method: 'Plate', args: [b1, b2, span1, span2, matName, 'Shield', nx, ny, nz, upA, upB, sideOff, weldDist], id });
        return id;
    },
    Slider: function (rail, calong, matName, sliding = null, slidingVertId = -1, side = 0, methodFlag = true) {
        const id = _counters.Slider++;
        _commands.push({ method: 'Slider', args: [rail, calong, matName, sliding, slidingVertId, side, methodFlag], id });
        return id;
    },
    Ring: function (pos, dir, up, R, nseg, wh, matName, st) {
        const id = _counters.Ring++;
        _commands.push({ method: 'Ring', args: [pos, dir, up, R, nseg, wh, matName, st], id });
        return id;
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
            _counters.Node = 0;
            _counters.Girder = 0;
            _counters.Rope = 0;
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
