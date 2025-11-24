const fs = require('fs');
const path = require('path');

// Mock global logger
global.logger = {
    info: (msg) => console.log(`[INFO] ${msg}`),
    debug: (msg) => console.log(`[DEBUG] ${msg}`),
    warn: (msg) => console.warn(`[WARN] ${msg}`),
    error: (msg) => console.error(`[ERROR] ${msg}`),
};

// Mock window (still needed for some other checks potentially, but logger is now global)
global.window = {
    logger: global.logger
};

// Load dependencies
// Adjust paths as necessary based on where this script is located relative to others
const { Vec3 } = require('../common_js/Vec3.js');
global.Vec3 = Vec3;

const { MeshBuilder } = require('./js/MeshBuilder.js');
const { WireFlags, extendMeshBuilder } = require('./js/MeshesUV.js');

// Extend MeshBuilder with MeshesUV functions
extendMeshBuilder(MeshBuilder);

// Extend MeshBuilder with block/skeleton generators
const { extendMeshBuilderWithGenerators } = require('./js/MeshGenerators.js');
extendMeshBuilderWithGenerators(MeshBuilder);

global.MeshBuilder = MeshBuilder;
global.WireFlags = WireFlags;

module.exports = {
    saveObj: function saveObj(filename, mesh) {
        // Validate mesh before saving
        const validation = mesh.validateMesh();

        if (!validation.valid) {
            logger.error(`Mesh validation failed for ${filename}:`);
            for (const err of validation.errors) {
                logger.error(`  - ${err}`);
            }
            if (validation.nanCount > validation.errors.length) {
                logger.error(`  - ... and ${validation.nanCount - validation.errors.length} more vertices with NaN`);
            }
            if (validation.infCount > validation.errors.length) {
                logger.error(`  - ... and ${validation.infCount - validation.errors.length} more vertices with Inf`);
            }
        }

        const objString = mesh.toObjString();
        const bytes = Buffer.from(objString, 'utf8').length;
        fs.writeFileSync(filename, objString);
        logger.info(`Saved ${filename} (${bytes} bytes)`);
    }
};
