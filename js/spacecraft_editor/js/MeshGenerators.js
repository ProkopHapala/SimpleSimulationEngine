import { Vec3 } from '../../common_js/Vec3.js';

// Simplified Skeleton and Block mesh generators
// These are simplified versions focused on mesh generation only

const MeshGenerator = {

    /**
     * Generate a skeleton (network of girders connecting nodes)
     * Simplified version that generates girder meshes between position pairs
     * 
     * @param {Array} nodes - Array of {pos: Vec3, size: number} 
     * @param {Array} connections - Array of {x, y} pairs (node indices)
     * @param {Object} params - {nseg, width, stickTypes}
     */
    skeleton(nodes, connections, params = {}) {
        const {
            nseg = 3,           // Segments along girder
            width = 0.5,        // Girder width
            stickTypes = { x: 0, y: 0, z: 0, w: 0 }
        } = params;

        logger.info(`skeleton: ${nodes.length} nodes, ${connections.length} connections`);

        // Generate girders between connected nodes  
        for (const conn of connections) {
            const n0 = nodes[conn.x];
            const n1 = nodes[conn.y];

            // Calculate up vector perpendicular to girder direction
            const dir = new Vec3().setSub(n1.pos, n0.pos);
            const len = dir.normalize();

            let up = new Vec3(0, 0, 1);
            if (Math.abs(dir.z) > 0.9) {
                up = new Vec3(1, 0, 0);
            }

            this.girder1(n0.pos, n1.pos, up, nseg, width, stickTypes, false);
        }

        logger.info(`skeleton: generated ${connections.length} girders`);
    },

    /**
     * Generate a simple block mesh (subdivided box)
     * 
     * @param {Vec3} pos - Center position
     * @param {Vec3} size - Box dimensions {x, y, z}
     * @param {Vec3} subdivs - Subdivisions per face {x, y, z}
     * @param {boolean} bWireframe - Generate wireframe edges
     */
    block(pos, size, subdivs = { x: 2, y: 2, z: 2 }, bWireframe = true) {
        const hx = size.x * 0.5;
        const hy = size.y * 0.5;
        const hz = size.z * 0.5;

        const nx = subdivs.x;
        const ny = subdivs.y;
        const nz = subdivs.z;

        logger.info(`block: size=(${size.x},${size.y},${size.z}), subdivs=(${nx},${ny},${nz})`);

        // Generate subdivided faces
        // We'll create vertices on a grid for each face and connect them

        const verts = [];

        // Helper to generate a grid of vertices on a face
        const makeFaceGrid = (corner, du, dv, nu, nv) => {
            const faceVerts = [];
            for (let iv = 0; iv <= nv; iv++) {
                const row = [];
                for (let iu = 0; iu <= nu; iu++) {
                    const p = new Vec3().setAdd(corner, new Vec3().setMul(du, new Vec3(iu / nu, iu / nu, iu / nu)));
                    p.addMul(dv, iv / nv);
                    const idx = this.vert(p);
                    row.push(idx);
                }
                faceVerts.push(row);
            }
            return faceVerts;
        };

        // Front face (z+)
        const front = makeFaceGrid(
            new Vec3(pos.x - hx, pos.y - hy, pos.z + hz),
            new Vec3(size.x, 0, 0),
            new Vec3(0, size.y, 0),
            nx, ny
        );

        // Back face (z-)
        const back = makeFaceGrid(
            new Vec3(pos.x + hx, pos.y - hy, pos.z - hz),
            new Vec3(-size.x, 0, 0),
            new Vec3(0, size.y, 0),
            nx, ny
        );

        // Connect with edges if wireframe
        if (bWireframe) {
            // Edges along each face
            const connectFaceGrid = (grid) => {
                for (let iv = 0; iv < grid.length; iv++) {
                    for (let iu = 0; iu < grid[iv].length; iu++) {
                        const curr = grid[iv][iu];
                        // Connect to right neighbor
                        if (iu < grid[iv].length - 1) {
                            this.edge(curr, grid[iv][iu + 1]);
                        }
                        // Connect to top neighbor
                        if (iv < grid.length - 1) {
                            this.edge(curr, grid[iv + 1][iu]);
                        }
                    }
                }
            };

            connectFaceGrid(front);
            connectFaceGrid(back);

            // Connect front to back
            for (let iv = 0; iv <= ny; iv++) {
                for (let iu = 0; iu <= nx; iu++) {
                    this.edge(front[iv][iu], back[iv][iu]);
                }
            }
        }

        logger.info(`block: generated subdivided block`);
    },

    /**
     * Generate a truss wheel akin to Mesh::Builder2::wheel:
     * - 4 vertices per segment (outer/inner at two successive angles with axial offset)
     * - longitudinal edges (stickTypes.x), rim spokes (stickTypes.y), diagonals (stickTypes.z), cross-bracing (stickTypes.w)
     */
    wheel(pos, p1, ax, nseg, wh = { x: 0.5, y: 0.5 }, stickTypes = { x: 1, y: 1, z: 1, w: 1 }) {
        const dir = new Vec3().setSub(p1, pos);
        const r = dir.norm();
        if (r <= 0) {
            logger.error('wheel: invalid radius');
            return -1;
        }
        const z = new Vec3(ax.x, ax.y, ax.z).normalize();
        // build orthonormal basis
        const tmp = Math.abs(z.x) < 0.9 ? new Vec3(1, 0, 0) : new Vec3(0, 1, 0);
        const x = new Vec3().setCross(tmp, z).normalize();
        const y = new Vec3().setCross(z, x).normalize();

        const dnp = 4; // verts per segment
        const i00 = this.verts.length / 3;

        const angStep = (2 * Math.PI) / nseg;
        let ang = 0;
        const cos = Math.cos;
        const sin = Math.sin;

        const pushVertAt = (a, radial, axial) => {
            const ca = cos(a), sa = sin(a);
            const Rv = new Vec3(
                x.x * (radial * ca) + y.x * (radial * sa),
                x.y * (radial * ca) + y.y * (radial * sa),
                x.z * (radial * ca) + y.z * (radial * sa)
            );
            return this.vert(new Vec3(
                pos.x + Rv.x + z.x * axial,
                pos.y + Rv.y + z.y * axial,
                pos.z + Rv.z + z.z * axial
            ));
        };

        for (let i = 0; i < nseg; i++, ang += angStep) {
            const angNext = ang + angStep;
            // four verts: outer/inner at ang, outer/inner at angNext with axial offsets
            const v0 = pushVertAt(ang,     r + wh.x, +0.0);
            const v1 = pushVertAt(ang,     r - wh.x, +0.0);
            const v2 = pushVertAt(angNext, r + wh.x, +wh.y);
            const v3 = pushVertAt(angNext, r - wh.x, -wh.y);

            // rim spokes
            this.edge(v0, v1, stickTypes.y);
            this.edge(v2, v3, stickTypes.y);
            // diagonals
            this.edge(v0, v2, stickTypes.z);
            this.edge(v0, v3, stickTypes.z);
            this.edge(v1, v2, stickTypes.z);
            this.edge(v1, v3, stickTypes.z);

            // connect to next segment or wrap
            const baseNext = (i === nseg - 1) ? i00 : i00 + (i + 1) * dnp;
            this.edge(v2, baseNext + 0, stickTypes.w);
            this.edge(v2, baseNext + 1, stickTypes.w);
            this.edge(v3, baseNext + 0, stickTypes.w);
            this.edge(v3, baseNext + 1, stickTypes.w);
            this.edge(v0, baseNext + 0, stickTypes.x);
            this.edge(v1, baseNext + 1, stickTypes.x);
            this.edge(v2, baseNext + 2, stickTypes.x);
            this.edge(v3, baseNext + 3, stickTypes.x);
        }
        logger.info(`wheel: generated truss wheel nseg=${nseg} r=${r}`);
        return i00;
    }
};

// Function to extend MeshBuilder with these methods
export function extendMeshBuilderWithGenerators(MeshBuilderClass) {
    Object.assign(MeshBuilderClass.prototype, MeshGenerator);
}

// Optional browser global for legacy usage
if (typeof window !== 'undefined') {
    window.extendMeshBuilderWithGenerators = extendMeshBuilderWithGenerators;

    // Auto-extend if MeshBuilder is already present
    if (window.MeshBuilder) {
        extendMeshBuilderWithGenerators(window.MeshBuilder);
    }
}

