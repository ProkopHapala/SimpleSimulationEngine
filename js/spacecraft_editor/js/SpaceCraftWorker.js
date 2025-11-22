
// SpaceCraftWorker.js

// Internal State for Shadow IDs
const _counters = {
    Node: 0,
    Girder: 0,
    Rope: 0
};

const _commands = [];

// The API exposed to the user script
const api = {
    Material: function (name, density, Spull, Spush) {
        _commands.push({ method: 'Material', args: [name, density, Spull, Spush] });
    },
    Node: function (pos) {
        const id = _counters.Node++;
        _commands.push({ method: 'Node', args: [pos] });
        return id;
    },
    Girder: function (n1, n2, matName) {
        const id = _counters.Girder++;
        _commands.push({ method: 'Girder', args: [n1, n2, matName] });
        return id;
    },
    Rope: function (n1, n2, thick, matName) {
        const id = _counters.Rope++;
        _commands.push({ method: 'Rope', args: [n1, n2, thick, matName] });
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
