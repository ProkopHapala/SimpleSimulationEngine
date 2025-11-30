import { logger } from './Logger.js';

export const WireFlags = {
    NONE: 0,
    CLOSED_CIRCUM: 1 << 0,  // 1
    CAPPED_CENTER: 1 << 1,  // 2
    RADIAL_EDGES: 1 << 2,   // 4
    AZIMUTHAL_EDGES: 1 << 3,// 8
    DIAGONAL1_EDGES: 1 << 4,// 16
    DIAGONAL2_EDGES: 1 << 5,// 32
    FIRST_RING: 1 << 6,     // 64
    LAST_RING: 1 << 7,      // 128
    ALTERNATE_DIAG: 1 << 8  // 256
};

WireFlags.BASIC_GRID = WireFlags.AZIMUTHAL_EDGES | WireFlags.RADIAL_EDGES;
WireFlags.FULL_GRID = WireFlags.BASIC_GRID | WireFlags.DIAGONAL1_EDGES | WireFlags.DIAGONAL2_EDGES;
WireFlags.DEFAULT_WIRE = WireFlags.BASIC_GRID | WireFlags.CLOSED_CIRCUM | WireFlags.FIRST_RING | WireFlags.LAST_RING;
WireFlags.STAR = WireFlags.CLOSED_CIRCUM | WireFlags.RADIAL_EDGES | WireFlags.FIRST_RING;
WireFlags.TRIMESH = WireFlags.BASIC_GRID | WireFlags.DIAGONAL1_EDGES | WireFlags.ALTERNATE_DIAG | WireFlags.FIRST_RING | WireFlags.LAST_RING | WireFlags.CLOSED_CIRCUM;

export const MeshesUV = {

    // --- UV Function Generators ---

    UVFunc2wire_new(n, UVmin, UVmax, voff, func, wire_flags = WireFlags.DEFAULT_WIRE) {
        const duv = { x: (UVmax.x - UVmin.x) / n.x, y: (UVmax.y - UVmin.y) / n.y };

        // Unpack flags
        const bPeriodicB = (wire_flags & WireFlags.CLOSED_CIRCUM) !== 0;
        const bHasCenter = (wire_flags & WireFlags.CAPPED_CENTER) !== 0;
        const bAzimEdges = (wire_flags & WireFlags.AZIMUTHAL_EDGES) !== 0;
        const bRadialEdges = (wire_flags & WireFlags.RADIAL_EDGES) !== 0;
        const bDiagEdges1 = (wire_flags & WireFlags.DIAGONAL1_EDGES) !== 0;
        const bDiagEdges2 = (wire_flags & WireFlags.DIAGONAL2_EDGES) !== 0;
        const bFirstRing = (wire_flags & WireFlags.FIRST_RING) !== 0;
        const bLastRing = (wire_flags & WireFlags.LAST_RING) !== 0;
        const bAlternateDiag = (wire_flags & WireFlags.ALTERNATE_DIAG) !== 0;

        const startIndex = this.verts.length;
        let centerIndex = -1;

        // Create center vertex if needed
        if (bHasCenter) {
            const p_center = func({ x: 0.0, y: 0.0 });
            centerIndex = this.vert(p_center);
        }

        // Create grid vertices
        const nCols = bPeriodicB ? n.y : n.y + 1;
        for (let ia = (bHasCenter ? 1 : 0); ia <= n.x; ia++) {
            let uv = { x: UVmin.x + duv.x * ia, y: UVmin.y + voff * duv.y * ia };
            for (let ib = 0; ib < nCols; ib++) {
                const p = func(uv);
                this.vert(p);
                uv.y += duv.y;
            }
        }

        if (typeof logger !== 'undefined') logger.debug(`UVFunc2wire_new: n.x=${n.x} n.y=${n.y} | flags=${wire_flags}`);

        // Create edges
        for (let ir = 0; ir <= n.x - (bHasCenter ? 1 : 0); ir++) {
            const rowStart = startIndex + (bHasCenter ? 1 : 0) + ir * nCols;

            // Connect to center for first ring
            if (bHasCenter && ir === 0) {
                for (let ib = 0; ib < nCols; ib++) {
                    this.edge(centerIndex, rowStart + ib);
                }
            }

            // Create edges within current row
            let bDoAzumith = bAzimEdges;
            if (ir === 0) bDoAzumith = bFirstRing;
            else if (ir === n.x - 1) bDoAzumith = bLastRing;

            if (bDoAzumith) {
                for (let ib = 0; ib < nCols - 1; ib++) {
                    this.edge(rowStart + ib, rowStart + ib + 1);
                }
                if (bPeriodicB) {
                    this.edge(rowStart + nCols - 1, rowStart);
                }
            }

            // Create edges between rows
            if (ir > 0 && bRadialEdges) {
                const prevRow = startIndex + (bHasCenter ? 1 : 0) + (ir - 1) * nCols;
                for (let ib = 0; ib < nCols; ib++) {
                    this.edge(prevRow + ib, rowStart + ib);
                }
            }

            // Create diagonal edges
            if (ir > 0 && (bDiagEdges1 || bDiagEdges2)) {
                const prevRow = startIndex + (bHasCenter ? 1 : 0) + (ir - 1) * nCols;
                let nc = nCols - 1;
                if (bPeriodicB) nc++;

                for (let ib = 0; ib < nc; ib++) {
                    let i0 = 0, i1 = 1;
                    if (bAlternateDiag && ((ir ^ ib) & 1)) { i0 = 1; i1 = 0; }

                    let idx_prev_0 = prevRow + ib + i1;
                    let idx_curr_0 = rowStart + ib + i0;

                    let idx_prev_1 = prevRow + ib + i0;
                    let idx_curr_1 = rowStart + ib + i1;

                    // Wrap indices if needed
                    if (bPeriodicB && ib === nCols - 1) {
                        if (i1 === 1) {
                            idx_prev_0 = prevRow + 0;
                            idx_curr_1 = rowStart + 0;
                        }
                        if (i0 === 1) {
                            idx_curr_0 = rowStart + 0;
                            idx_prev_1 = prevRow + 0;
                        }
                    }

                    if (bDiagEdges1) this.edge(idx_prev_0, idx_curr_0);
                    if (bDiagEdges2) this.edge(idx_prev_1, idx_curr_1);
                }
            }
        }
    },

    // Parabolic dish using ParametricQuadPatch as a UV generator. We first
    // build a strip over [0,1]x[0,1] in (u,v) using ParametricQuadPatch, then
    // reinterpret (x,y) as (u,v) and map to a paraboloid via ParabolaUVfunc.
    // - nTop:    segment count on the OUTER ring (radius R2)
    // - nBottom: segment count on the INNER ring (radius R1)
    // - nRows:   number of radial rows from inner to outer radius
    // - R1:      inner radius
    // - R2:      outer radius of the paraboloid (used for curvature)
    // - L:       axial height at radius R2
ParametricParabolaPatch(nTop, nBottom, nRows, R1, R2, L) {
    const iv0 = this.verts.length;

    const p00 = { x: 0.0, y: 0.0, z: 0.0 };
    const p01 = { x: 0.0, y: 1.0, z: 0.0 };
    const p10 = { x: 1.0, y: 0.0, z: 0.0 };
    const p11 = { x: 1.0, y: 1.0, z: 0.0 };
    // NOTE: We want nBottom on the INNER ring and nTop on the OUTER ring.
    // ParametricQuadPatch uses nTop for the FIRST row and nBottom for the LAST.
    // Therefore we swap the order here so that row 0 (inner) gets nBottom and
    // row nRows-1 (outer) gets nTop.
    this.ParametricQuadPatch(nBottom, nTop, nRows, p00, p01, p10, p11);

    const iv1 = this.verts.length;

    // DEBUG: basic info
    console.log("ParametricParabolaPatch: header", {
        iv0, iv1,
        nTop, nBottom, nRows,
        R1, R2, L,
        addedVerts: iv1 - iv0
    });

    if (iv1 <= iv0) {
        console.warn("ParametricParabolaPatch: no verts added by ParametricQuadPatch");
        return;
    }

    // DEBUG: show first vertex before remap
    const v0_before = { ...this.verts[iv0] };
    console.log("ParametricParabolaPatch: v0 before remap", v0_before);

    // Remap (x,y) stored in vertex position -> (u_rad, u_ang) -> paraboloid using ParabolaUVfunc
    for (let i = iv0; i < iv1; i++) {
        const v = this.verts[i];
        // ParametricQuadPatch stored the strip in v.pos.x (along row), v.pos.y (row index)
        // We interpret v.pos.y as radial coordinate (0..1 inner->outer) and v.pos.x as
        // angular fraction around the ring.
        const u_rad = v.pos.y;              // 0..1 radial fraction (inner -> outer)
        const u_ang = v.pos.x;              // 0..1 angular fraction (around ring)
        const r     = R1 + (R2 - R1) * u_rad;
        const ang   = u_ang * 2.0 * Math.PI;
        const p     = this.ParabolaUVfunc({ x: r, y: ang }, R2, L);
        v.pos.x = p.x;
        v.pos.y = p.y;
        v.pos.z = p.z;
    }

    // DEBUG: show first vertex after remap
    const v0_after = { ...this.verts[iv0] };
    console.log("ParametricParabolaPatch: v0 after remap", v0_after);
},

    // Wrap a ParametricQuadPatch (built over a unit square in x,y) into an
    // annular / disk-like patch using cylindrical coordinates. We treat the
    // generated x,y as (u_radial, u_angle) in [0,1]^2 and remap to:
    //   r   = R_inner + u_radial * (R_outer - R_inner)
    //   phi = phi0 + u_angle * dphi
    //   p   = (r*cos(phi), r*sin(phi), z0)
    // Connectivity (including variable nTop/nBottom) is entirely reused from
    // ParametricQuadPatch; only vertex positions are updated afterwards.
    ParametricAnnulusPatch(nTop, nBottom, nRows, R_inner, R_outer, phi0 = 0.0, dphi = 2.0 * Math.PI, z0 = 0.0) {
        const Rin = R_inner;
        const Rout = R_outer;
        const dR = Rout - Rin;

        // 1) Build a parametric strip over unit square [0,1]x[0,1] in x,y
        const iv0 = this.verts.length;
        const p00 = { x: 0.0, y: 0.0, z: 0.0 };
        const p01 = { x: 0.0, y: 1.0, z: 0.0 };
        const p10 = { x: 1.0, y: 0.0, z: 0.0 };
        const p11 = { x: 1.0, y: 1.0, z: 0.0 };
        this.ParametricQuadPatch(nTop, nBottom, nRows, p00, p01, p10, p11);
        const iv1 = this.verts.length;

        // 2) Remap positions to cylindrical annulus
        for (let i = iv0; i < iv1; i++) {
            const v = this.verts[i];
            const u_rad = v.x;  // 0..1 across thickness (inner->outer)
            const u_ang = v.y;  // 0..1 around the ring
            const r   = Rin + dR * u_rad;
            const phi = phi0 + dphi * u_ang;
            const cs  = Math.cos(phi);
            const sn  = Math.sin(phi);
            v.x = cs * r;
            v.y = sn * r;
            v.z = z0;
        }
    },

    Parabola_Wire_new(n, UVmin, UVmax, R, L, voff, wire_flags = WireFlags.DEFAULT_WIRE) {
        const K = L / (R * R);
        // UVmin/max modification
        const uvMinMod = { x: UVmin.x * R, y: UVmin.y };
        const uvMaxMod = { x: UVmax.x * R, y: UVmax.y };

        const uvfunc = (uv) => {
            // ParabolaUVfunc
            // uv.x is r, uv.y is angle
            const ang = uv.y;
            const csb = { x: Math.cos(ang), y: Math.sin(ang) };
            const r = uv.x;
            const l = r * r * K;
            return new Vec3(csb.x * r, csb.y * r, l);
        };

        this.UVFunc2wire_new(n, uvMinMod, uvMaxMod, voff, uvfunc, wire_flags);
    },

    // --- Slab & Panel Generators ---

    QuadUVfunc(uv, p00, p01, p10, p11) {
        const u = uv.x;
        const v = uv.y;
        const mu = 1 - u;
        const mv = 1 - v;

        // Bilinear interpolation: p00 * (mu * mv) + p01 * (mu * v) + p10 * (u * mv) + p11 * (u * v)
        const res = new Vec3();
        res.setLincomb3(mu * mv, p00, mu * v, p01, u * mv, p10);
        res.addMul(p11, u * v);

        return res;
    },

    UV_slab_verts(n, uv0, duv, idx, func) {
        const iv0 = this.verts.length;
        for (let iy = 0; iy < n.y; ++iy) {
            for (let ix = 0; ix < n.x; ++ix) {
                const uv = { x: uv0.x + ix * duv.x, y: uv0.y + iy * duv.y };
                const p = func(uv);
                const iFlat = iy * n.x + ix;
                idx[iFlat] = this.vert(p);
            }
        }
        return iv0;
    },

    slabEdges(n, idx0, idx1, dirMask, stickTypes) {
        // Unique positive directions (13)
        const DIRS = [
            // axial
            { x: 1, y: 0, z: 0 }, // 0 x
            { x: 0, y: 1, z: 0 }, // 1 y
            { x: 0, y: 0, z: 1 }, // 2 z
            // face diag
            { x: 1, y: 1, z: 0 }, // 3 x+y
            { x: 1, y: -1, z: 0 }, // 4 x-y
            { x: 1, y: 0, z: 1 }, // 5 x+z
            { x: -1, y: 0, z: 1 }, // 6 z-x
            { x: 0, y: 1, z: 1 }, // 7 y+z
            { x: 0, y: 1, z: -1 }, // 8 y-z
            // space diag
            { x: 1, y: 1, z: 1 }, // 9 x+y+z
            { x: 1, y: 1, z: -1 }, // 10 x+y-z
            { x: 1, y: -1, z: 1 }, // 11 x-y+z
            { x: 1, y: -1, z: -1 }  // 12 x-y-z
        ];
        const nz = 2;
        const getVert = (ix, iy, iz) => {
            const i = iy * n.x + ix;
            return (iz === 0) ? idx0[i] : idx1[i];
        };
        const edgeTypeByDir = (d) => {
            const comps = (d.x !== 0 ? 1 : 0) + (d.y !== 0 ? 1 : 0) + (d.z !== 0 ? 1 : 0);
            if (comps === 1) { return d.z ? stickTypes.x : stickTypes.y; }         // axis
            if (comps === 2) { return d.z ? stickTypes.w : stickTypes.z; }         // face diag
            return stickTypes.w;                                              // space diag
        };

        const ie0 = this.edges.length;
        for (let iz = 0; iz < nz; ++iz) {
            for (let iy = 0; iy < n.y; ++iy) {
                for (let ix = 0; ix < n.x; ++ix) {
                    for (let id = 0; id < 13; ++id) {
                        if (!(dirMask & (1 << id))) continue;
                        const d = DIRS[id];
                        const jx = ix + d.x;
                        const jy = iy + d.y;
                        const jz = iz + d.z;
                        if (jx < 0 || jx >= n.x || jy < 0 || jy >= n.y || jz < 0 || jz >= nz) continue; // stay inside slab
                        const a = getVert(ix, iy, iz);
                        const b = getVert(jx, jy, jz);
                        this.edge(a, b, edgeTypeByDir(d));
                    }
                }
            }
        }
        return ie0;
    },

    // Periodic-in-y variant of slabEdges for tubes: y index wraps modulo n.y,
    // while x (axial) and z (layer) stay clamped. This removes duplicate seam
    // vertices in the angular direction and matches TubeSheet's periodic
    // topology when interpreting y as angle.
    slabEdges_wrap(n, idx0, idx1, dirMask, stickTypes) {
        const DIRS = [
            // axial
            { x: 1, y: 0, z: 0 }, // 0 x
            { x: 0, y: 1, z: 0 }, // 1 y
            { x: 0, y: 0, z: 1 }, // 2 z
            // face diag
            { x: 1, y: 1, z: 0 }, // 3 x+y
            { x: 1, y: -1, z: 0 }, // 4 x-y
            { x: 1, y: 0, z: 1 }, // 5 x+z
            { x: -1, y: 0, z: 1 }, // 6 z-x
            { x: 0, y: 1, z: 1 }, // 7 y+z
            { x: 0, y: 1, z: -1 }, // 8 y-z
            // space diag
            { x: 1, y: 1, z: 1 }, // 9 x+y+z
            { x: 1, y: 1, z: -1 }, // 10 x+y-z
            { x: 1, y: -1, z: 1 }, // 11 x-y+z
            { x: 1, y: -1, z: -1 }  // 12 x-y-z
        ];
        const nz = 2;
        const wrapY = (iy) => {
            const r = iy % n.y;
            return r < 0 ? r + n.y : r;
        };
        const getVert = (ix, iy, iz) => {
            const i = iy * n.x + ix;
            return (iz === 0) ? idx0[i] : idx1[i];
        };
        const edgeTypeByDir = (d) => {
            const comps = (d.x !== 0 ? 1 : 0) + (d.y !== 0 ? 1 : 0) + (d.z !== 0 ? 1 : 0);
            if (comps === 1) { return d.z ? stickTypes.x : stickTypes.y; }
            if (comps === 2) { return d.z ? stickTypes.w : stickTypes.z; }
            return stickTypes.w;
        };

        const ie0 = this.edges.length;
        for (let iz = 0; iz < nz; ++iz) {
            for (let iy = 0; iy < n.y; ++iy) {
                for (let ix = 0; ix < n.x; ++ix) {
                    for (let id = 0; id < 13; ++id) {
                        if (!(dirMask & (1 << id))) continue;
                        const d = DIRS[id];
                        const jx = ix + d.x;
                        const jyRaw = iy + d.y;
                        const jz = iz + d.z;
                        if (jx < 0 || jx >= n.x || jz < 0 || jz >= nz) continue; // x,z clamped
                        const jy = wrapY(jyRaw);
                        const a = getVert(ix, iy, iz);
                        const b = getVert(jx, jy, jz);
                        this.edge(a, b, edgeTypeByDir(d));
                    }
                }
            }
        }
        return ie0;
    },

    UV_slab(n, UVmin, UVmax, up, dirMask, stickTypes, func1, func2) {
        const duv = {
            x: (UVmax.x - UVmin.x) / (n.x - 1),
            y: (UVmax.y - UVmin.y) / (n.y - 1)
        };
        const idx0 = new Array(n.x * n.y);
        const idx1 = new Array(n.x * n.y);

        // Note: In C++, UVmin+duv*up.xy() is passed to second call.
        // up is Vec3f, up.xy() is Vec2f.
        // JS: up.x, up.y used for UV offset?
        // C++: UVmin+duv*up.xy() -> UVmin.x + duv.x*up.x, UVmin.y + duv.y*up.y
        const uvOffset = { x: duv.x * up.x, y: duv.y * up.y };
        const uvStart2 = { x: UVmin.x + uvOffset.x, y: UVmin.y + uvOffset.y };

        this.UV_slab_verts(n, UVmin, duv, idx0, func1);
        this.UV_slab_verts(n, uvStart2, duv, idx1, func2);
        this.slabEdges(n, idx0, idx1, dirMask, stickTypes);
    },

    // Periodic variant of UV_slab for tubes: angular direction (y) is treated as
    // periodic, so there are exactly n.y unique azimuthal samples and no
    // duplicate seam vertices. Connectivity wraps in index space instead.
    UV_slab_wrap(n, UVmin, UVmax, up, dirMask, stickTypes, func1, func2) {
        const duv = {
            x: (UVmax.x - UVmin.x) / (n.x - 1),
            y: (UVmax.y - UVmin.y) / n.y
        };
        const idx0 = new Array(n.x * n.y);
        const idx1 = new Array(n.x * n.y);

        const uvOffset = { x: duv.x * up.x, y: duv.y * up.y };
        const uvStart2 = { x: UVmin.x + uvOffset.x, y: UVmin.y + uvOffset.y };

        this.UV_slab_verts(n, UVmin, duv, idx0, func1);
        this.UV_slab_verts(n, uvStart2, duv, idx1, func2);
        this.slabEdges_wrap(n, idx0, idx1, dirMask, stickTypes);
    },

    QuadSlab(n, UVmin, UVmax, p00, p01, p10, p11, up, dirMask, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        // Compute normal: cross(p10-p00, p01-p00)
        const d1 = new Vec3().setSub(p10, p00);
        const d2 = new Vec3().setSub(p01, p00);
        const nor = new Vec3().setCross(d1, d2);
        nor.normalize();

        const uvfunc1 = (uv) => this.QuadUVfunc(uv, p00, p01, p10, p11);
        const uvfunc2 = (uv) => {
            const p = this.QuadUVfunc(uv, p00, p01, p10, p11);
            p.addMul(nor, up.z);
            return p;
        };

        this.UV_slab(n, UVmin, UVmax, up, dirMask, stickTypes, uvfunc1, uvfunc2);
    },

    // Periodic angular variant of SlabTube. Interprets n.y as the number of
    // unique azimuthal directions (polygon sides). Uses UV_slab_wrap and
    // slabEdges_wrap so there is a single seam in index space, without
    // duplicated vertices at angles 0 and 2π.
    SlabTube_wrap(n, UVmin, UVmax, Rs, L, up, dirMask, twist = 0.5, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const dudv = twist * (n.x - 1.0) / n.y;

        const uvfunc1 = (uv) => {
            const uv_y = (uv.y + uv.x * dudv) * 2 * Math.PI;
            return this.ConeUVfunc({ x: uv.x, y: uv_y }, Rs.x, Rs.y, L);
        };
        const uvfunc2 = (uv) => {
            const uv_y = (uv.y + uv.x * dudv) * 2 * Math.PI;
            return this.ConeUVfunc({ x: uv.x, y: uv_y }, Rs.x + up.z, Rs.y + up.z, L);
        };

        this.UV_slab_wrap(n, UVmin, UVmax, up, dirMask, stickTypes, uvfunc1, uvfunc2);
    },

    // Parabolic analogue of SlabTube_wrap: double-layer paraboloid with slab connectivity
    // r in [0, Rs.x], angle in [0, 2π], z = (L / (Rs.x^2)) * r^2 for both layers,
    // with the outer layer shifted by up.z along +z to give thickness.
    ParabolaSlab_wrap(n, UVmin, UVmax, Rs, L, up, dirMask, twist = 0.5, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const dudv = twist * (n.x - 1.0) / n.y;
        const R = Rs.x;
        const K = (R > 0) ? (L / (R * R)) : 0.0;

        const baseParabola = (uv, zOffset) => {
            // uv.x ~ radial parameter in [0,1], uv.y ~ angular parameter in [0,1]
            const r = R * uv.x;
            const ang = (uv.y + uv.x * dudv) * 2.0 * Math.PI;
            const cs = Math.cos(ang);
            const sn = Math.sin(ang);
            const z = K * r * r + zOffset;
            return new Vec3(cs * r, sn * r, z);
        };

        const uvfunc1 = (uv) => baseParabola(uv, 0.0);
        const uvfunc2 = (uv) => baseParabola(uv, up.z);

        this.UV_slab_wrap(n, UVmin, UVmax, up, dirMask, stickTypes, uvfunc1, uvfunc2);
    },

    UV_panel(n, UVmin, UVmax, width, stickTypes, func) {
        const na = n.x;
        const nb = n.y;
        const du = { x: (UVmax.x - UVmin.x) / (na - 1), y: 0.0 };
        const dv = { x: 0.0, y: (UVmax.y - UVmin.y) / (nb - 1) };

        const base0 = this.verts.length;

        // 1) base grid vertices
        const baseIdx = new Array(na * nb);
        for (let ib = 0; ib < nb; ++ib) {
            for (let ia = 0; ia < na; ++ia) {
                const uv = { x: UVmin.x + ia * du.x + ib * dv.x, y: UVmin.y + ia * du.y + ib * dv.y };
                const idx = this.verts.length;
                baseIdx[ib * na + ia] = idx;
                this.vert(func(uv));
            }
        }

        // 2) center layer (raised) vertices
        const cenIdx = new Array((na - 1) * (nb - 1));
        for (let ib = 0; ib < nb - 1; ++ib) {
            for (let ia = 0; ia < na - 1; ++ia) {
                const uv00 = { x: UVmin.x + ia * du.x + ib * dv.x, y: UVmin.y + ia * du.y + ib * dv.y };
                const uv10 = { x: UVmin.x + (ia + 1) * du.x + ib * dv.x, y: UVmin.y + (ia + 1) * du.y + ib * dv.y };
                const uv01 = { x: UVmin.x + ia * du.x + (ib + 1) * dv.x, y: UVmin.y + ia * du.y + (ib + 1) * dv.y };

                const p00 = func(uv00);
                const p10 = func(uv10);
                const p01 = func(uv01);

                // bilinear center
                const uvC = { x: uv00.x + 0.5 * du.x + 0.5 * dv.x, y: uv00.y + 0.5 * du.y + 0.5 * dv.y };
                const pc = func(uvC);

                // normal
                const uvec = new Vec3().setSub(p10, p00);
                const vvec = new Vec3().setSub(p01, p00);
                const n = new Vec3().setCross(uvec, vvec);
                n.normalize();
                pc.addMul(n, width);

                const idx = this.verts.length;
                cenIdx[ib * (na - 1) + ia] = idx;
                this.vert(pc);
            }
        }

        // 3) base grid edges
        for (let ib = 0; ib < nb; ++ib) {
            for (let ia = 0; ia < na; ++ia) {
                const v = baseIdx[ib * na + ia];
                if (ia < na - 1) this.edge(v, baseIdx[ib * na + ia + 1], stickTypes.y);
                if (ib < nb - 1) this.edge(v, baseIdx[(ib + 1) * na + ia], stickTypes.y);
            }
        }

        // 4) connect center vertices
        for (let ib = 0; ib < nb - 1; ++ib) {
            for (let ia = 0; ia < na - 1; ++ia) {
                const c = cenIdx[ib * (na - 1) + ia];
                const v00 = baseIdx[ib * na + ia];
                const v10 = baseIdx[ib * na + ia + 1];
                const v01 = baseIdx[(ib + 1) * na + ia];
                const v11 = baseIdx[(ib + 1) * na + ia + 1];

                this.edge(c, v00, stickTypes.w);
                this.edge(c, v10, stickTypes.w);
                this.edge(c, v01, stickTypes.w);
                this.edge(c, v11, stickTypes.w);

                if (ia < na - 2) {
                    this.edge(c, cenIdx[ib * (na - 1) + ia + 1], stickTypes.z);
                }
                if (ib < nb - 2) {
                    this.edge(c, cenIdx[(ib + 1) * (na - 1) + ia], stickTypes.z);
                }
            }
        }

        return base0;
    },

    QuadPanel(n, UVmin, UVmax, p00, p01, p10, p11, width, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const uvfunc = (uv) => this.QuadUVfunc(uv, p00, p01, p10, p11);
        this.UV_panel(n, UVmin, UVmax, width, stickTypes, uvfunc);
    },

    // --- New Functions for Sheets and Tubes ---

    ConeUVfunc(uv, R1, R2, L) {
        // uv.x = u (0..1 along length), uv.y = v (angle)
        const ang = uv.y;
        const csb = { x: Math.cos(ang), y: Math.sin(ang) };
        const R = (1 - uv.x) * R1 + uv.x * R2;
        return new Vec3(csb.x * R, csb.y * R, L * uv.x);
    },

    // Parabolic analogue of ConeUVfunc used by parabolic generators.
    // uv.x = r (radius), uv.y = angle; z = (L / R^2) * r^2
    ParabolaUVfunc(uv, R, L) {
        const ang = uv.y;
        const csb = { x: Math.cos(ang), y: Math.sin(ang) };
        const r   = uv.x;
        const K   = (R > 0) ? (L / (R * R)) : 0.0;
        const z   = K * r * r;
        return new Vec3(csb.x * r, csb.y * r, z);
    },

    TorusUVfunc(uv, r, R) {
        // uv.x = angle small circle, uv.y = angle large circle
        const csa = { x: Math.cos(uv.x), y: Math.sin(uv.x) };
        const csb = { x: Math.cos(uv.y), y: Math.sin(uv.y) };
        return new Vec3(csb.x * (R + r * csa.x), csb.y * (R + r * csa.x), r * csa.y);
    },

    stickEdges2D(n, idx, dirMask, stickTypes, bPeriodicX = false, bPeriodicY = false) {
        const DIRS = [
            { x: 1, y: 0 }, // 0
            { x: 0, y: 1 }, // 1
            { x: 1, y: 1 }, // 2
            { x: 1, y: -1 } // 3
        ];
        const ie0 = this.edges.length;
        for (let iy = 0; iy < n.y; ++iy) {
            for (let ix = 0; ix < n.x; ++ix) {
                for (let id = 0; id < 4; ++id) {
                    if (!(dirMask & (1 << id))) continue;
                    const d = DIRS[id];
                    let jx = ix + d.x;
                    let jy = iy + d.y;

                    if (bPeriodicX) { jx = (jx + n.x) % n.x; } else if (jx < 0 || jx >= n.x) continue;
                    if (bPeriodicY) { jy = (jy + n.y) % n.y; } else if (jy < 0 || jy >= n.y) continue;

                    const a = idx[iy * n.x + ix];
                    const b = idx[jy * n.x + jx];

                    // Edge type logic: if diagonal (d.x && d.y) use stickTypes.z, else stickTypes.y
                    const type = (d.x !== 0 && d.y !== 0) ? stickTypes.z : stickTypes.y;
                    this.edge(a, b, type);
                }
            }
        }
        return ie0;
    },

    UV_sheet(n, uv0, duv, dirMask, stickTypes, func, bPeriodicX = false, bPeriodicY = false) {
        const idx = new Array(n.x * n.y);
        const iv0 = this.UV_slab_verts(n, uv0, duv, idx, func);
        this.stickEdges2D(n, idx, dirMask, stickTypes, bPeriodicX, bPeriodicY);
        return iv0;
    },

    SlabTube(n, UVmin, UVmax, Rs, L, up, dirMask, twist = 0.5, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const dudv = twist * (n.x - 1.0) / (n.y - 1.0);

        const uvfunc1 = (uv) => {
            const uv_y = (uv.y + uv.x * dudv) * 2 * Math.PI;
            return this.ConeUVfunc({ x: uv.x, y: uv_y }, Rs.x, Rs.y, L);
        };
        const uvfunc2 = (uv) => {
            const uv_y = (uv.y + uv.x * dudv) * 2 * Math.PI;
            return this.ConeUVfunc({ x: uv.x, y: uv_y }, Rs.x + up.z, Rs.y + up.z, L);
        };

        this.UV_slab(n, UVmin, UVmax, up, dirMask, stickTypes, uvfunc1, uvfunc2);
    },

    QuadSheet(n, UVmin, UVmax, p00, p01, p10, p11, dirMask, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const uvfunc = (uv) => this.QuadUVfunc(uv, p00, p01, p10, p11);
        const duv = { x: (UVmax.x - UVmin.x) / (n.x - 1), y: (UVmax.y - UVmin.y) / (n.y - 1) };
        this.UV_sheet(n, UVmin, duv, dirMask, stickTypes, uvfunc, false, false);
    },

    // Parametric quad patch with varying vertex counts per row, following the
    // Python generate_mesh_variations() + triangulate_strip(..., mode='parametric')
    // logic. This works directly in world space, not UV space.
    // - nTop:    vertex count on the first row (between p00 and p10)
    // - nBottom: vertex count on the last row  (between p01 and p11)
    // - nRows:   number of rows (>=2)
    ParametricQuadPatch(nTop, nBottom, nRows, p00, p01, p10, p11) {
        nTop = Math.max(1, nTop | 0);
        nBottom = Math.max(1, nBottom | 0);
        nRows = Math.max(2, nRows | 0);

        const rowCounts = [];
        for (let i = 0; i < nRows; i++) {
            const t = (nRows === 1) ? 0 : i / (nRows - 1);
            const cnt = Math.round(nTop + (nBottom - nTop) * t);
            rowCounts.push(Math.max(1, cnt));
        }

        const allRows = [];

        const lerp = (a, b, t) => ({
            x: a.x * (1 - t) + b.x * t,
            y: a.y * (1 - t) + b.y * t,
            z: a.z * (1 - t) + b.z * t
        });

        // Build rows of vertices between the two edges p00-p01 and p10-p11
        for (let i = 0; i < nRows; i++) {
            const v = (nRows === 1) ? 0 : i / (nRows - 1);
            const left = lerp(p00, p01, v);
            const right = lerp(p10, p11, v);
            const count = rowCounts[i];
            const row = [];

            if (count === 1) {
                row.push(this.vert(left));
            } else {
                for (let j = 0; j < count; j++) {
                    const u = j / (count - 1);
                    const pos = lerp(left, right, u);
                    row.push(this.vert(pos));
                }
            }
            allRows.push(row);
        }

        // Edges along each row
        for (let i = 0; i < nRows; i++) {
            const row = allRows[i];
            for (let j = 0; j < row.length - 1; j++) {
                this.edge(row[j], row[j + 1]);
            }
        }

        // Parametric triangulation between successive rows
        for (let r = 0; r < nRows - 1; r++) {
            const rowA = allRows[r];
            const rowB = allRows[r + 1];
            const countA = rowA.length;
            const countB = rowB.length;

            let ia = 0;
            let ib = 0;

            while (ia < countA - 1 || ib < countB - 1) {
                let choice;

                if (ia === countA - 1) {
                    choice = 'B';
                } else if (ib === countB - 1) {
                    choice = 'A';
                } else {
                    const uNextA = (ia + 1) / (countA - 1);
                    const uNextB = (ib + 1) / (countB - 1);
                    choice = (uNextA <= uNextB) ? 'A' : 'B';
                }

                if (choice === 'A') {
                    // Triangle (A_curr, B_curr, A_next)
                    const aCurr = rowA[ia];
                    const aNext = rowA[ia + 1];
                    const bCurr = rowB[ib];

                    this.edge(aCurr, bCurr);
                    this.edge(bCurr, aNext);
                    this.edge(aCurr, aNext);

                    ia += 1;
                } else {
                    // Triangle (A_curr, B_curr, B_next)
                    const aCurr = rowA[ia];
                    const bCurr = rowB[ib];
                    const bNext = rowB[ib + 1];

                    this.edge(aCurr, bCurr);
                    this.edge(bCurr, bNext);
                    this.edge(aCurr, bNext);

                    ib += 1;
                }
            }
        }
    },

    TubeSheet(n, UVmin, UVmax, Rs, L, dirMask = 0b1011, twist = 0.5, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const dudv = (twist * (n.x - 1.0)) / n.y;
        const uvfunc = (uv) => {
            const uv_y = (uv.y + uv.x * dudv) * 2 * Math.PI;
            const p = this.ConeUVfunc({ x: uv.x, y: uv_y }, Rs.x, Rs.y, L);
            return p;
        };
        // use n.x for periodic wrap and n.y-1 for length spacing?
        // C++: Vec2f duv = { (UVmax.x-UVmin.x)/(n.x-1), (UVmax.y-UVmin.y)/(n.y) }; 
        const duv = { x: (UVmax.x - UVmin.x) / (n.x - 1), y: (UVmax.y - UVmin.y) / n.y };
        this.UV_sheet(n, UVmin, duv, dirMask, stickTypes, uvfunc, false, true);
    },

    // Parabolic analogue of TubeSheet: same connectivity/topology, but mapped to a paraboloid
    // r in [0, Rs.x], angle in [0, 2π], z = (L / (Rs.x^2)) * r^2
    ParabolaSheet(n, UVmin, UVmax, Rs, L, dirMask = 0b1011, twist = 0.5, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const dudv = (twist * (n.x - 1.0)) / n.y;
        const R = Rs.x;
        const K = (R > 0) ? (L / (R * R)) : 0.0;
        const uvfunc = (uv) => {
            // uv.x ∈ [0,1] -> radius r ∈ [0,R]
            const r = R * uv.x;
            const ang = (uv.y + uv.x * dudv) * 2.0 * Math.PI;
            const cs = Math.cos(ang);
            const sn = Math.sin(ang);
            const z = K * r * r;
            return new Vec3(cs * r, sn * r, z);
        };
        const duv = { x: (UVmax.x - UVmin.x) / (n.x - 1), y: (UVmax.y - UVmin.y) / n.y };
        this.UV_sheet(n, UVmin, duv, dirMask, stickTypes, uvfunc, false, true);
    },

    // Variant of TubeSheet where n.x is treated as azimuthal (ring) direction.
    // Consecutive ix indices (0,1,2,...) run around the ring; iy runs along axis.
    TubeSheet_swapped(n, UVmin, UVmax, Rs, L, dirMask = 0b1011, twist = 0.5, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const nx = n.x | 0;
        const ny = n.y | 0;
        if (nx <= 0 || ny <= 0) return;
        // Interpret: uv.x ~ angle, uv.y ~ axial parameter u (0..1)
        const dudv = (twist * (ny - 1.0)) / nx;
        const uvfunc = (uv) => {
            // Map uv.y (0..1) to axial u, uv.x (0..1) to angle.
            const u  = (UVmax.y === UVmin.y) ? uv.y : (uv.y - UVmin.y) / (UVmax.y - UVmin.y);
            const ang = (uv.x + u * dudv) * 2.0 * Math.PI;
            return this.ConeUVfunc({ x: u, y: ang }, Rs.x, Rs.y, L);
        };
        const duv = {
            x: (UVmax.x - UVmin.x) / nx,      // step around ring
            y: (UVmax.y - UVmin.y) / (ny - 1) // step along axis
        };
        // Periodic in X (around ring), non-periodic in Y (axis).
        this.UV_sheet(n, UVmin, duv, dirMask, stickTypes, uvfunc, true, false);
    },

    // Post-process a regular grid (e.g. TubeSheet) to add "skip" edges.
    // - n: {x: nAlongX, y: nAlongY}
    // - iv0: starting vertex index of this grid in this.verts
    // - xskip: stride along X between connected vertices (can be negative)
    // - yskip: stride along Y between connected vertices, wrapped modulo n.y (can be negative)
    // - matID: edge type / material ID to store in edge.z
    addSkipEdges(n, iv0, xskip, yskip, matID = 0) {
        const nx = n.x | 0;
        const ny = n.y | 0;
        if (nx <= 0 || ny <= 0) return;
        const wrap = (k, m) => {
            const r = k % m;
            return r < 0 ? r + m : r;
        };
        if (xskip === 0) return;
        // For a given (ix, iy), connect to (ix + xskip, iy + yskip mod ny)
        for (let ix = 0; ix < nx; ix++) {
            const jx = wrap(ix + xskip, nx);
            if (jx < 0 || jx >= nx) continue; // stay within axial range
            for (let iy = 0; iy < ny; iy++) {
                //const jy = wrap(iy + yskip, ny);
                const jy = iy + yskip;
                if (jy < 0 || jy >= ny) continue; // stay within axial range
                const a = iv0 + iy * nx + ix;
                const b = iv0 + jy * nx + jx;
                this.edge(a, b, matID);
            }
        }
    },

    TorusSheet(n, UVmin, UVmax, Rs, dirMask = 0b1011, twist = 0.5, stickTypes = { x: 0, y: 0, z: 0, w: 0 }) {
        const dudv = (twist * n.y) / n.x;
        const uvfunc = (uv) => {
            const uv_x = (uv.x + uv.y * dudv) * 2 * Math.PI;
            const uv_y = uv.y * 2 * Math.PI;
            return this.TorusUVfunc({ x: uv_x, y: uv_y }, Rs.x, Rs.y);
        };
        // C++: Vec2f duv = { (UVmax.x-UVmin.x)/(n.x), (UVmax.y-UVmin.y)/(n.y) };
        const duv = { x: (UVmax.x - UVmin.x) / n.x, y: (UVmax.y - UVmin.y) / n.y };
        this.UV_sheet(n, UVmin, duv, dirMask, stickTypes, uvfunc, true, true);
    }

};

// Function to mixin these methods into MeshBuilder
export function extendMeshBuilder(MeshBuilderClass) {
    Object.assign(MeshBuilderClass.prototype, MeshesUV);
}

// Also expose on window for legacy global-script users
if (typeof window !== 'undefined') {
    window.WireFlags = WireFlags;
    window.extendMeshBuilderWithUV = extendMeshBuilder;

    // Auto-extend if MeshBuilder is already present
    if (window.MeshBuilder) {
        extendMeshBuilder(window.MeshBuilder);
    }
}
