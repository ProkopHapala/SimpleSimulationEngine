#ifndef DrawUV_h
#define DrawUV_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "MeshBuilder2.h"
#include "geometry/UVfuncs.h"

namespace Mesh {

// enum class WireFlags : uint16_t {
//     NONE            = 0,
//     CLOSED_CIRCUM   = 1 << 0,  // 1: Close circumference (periodic in angular direction)
//     CAPPED_CENTER   = 1 << 1,  // 2: Add central vertex and connect to first ring
//     RADIAL_EDGES    = 1 << 2,  // 4: Add radial edges between rings
//     AZIMUTHAL_EDGES = 1 << 3,  // 8: Add edges along each ring
//     DIAGONAL1_EDGES = 1 << 4,  // 16: Add diagonal edges (top-right to bottom-left)
//     DIAGONAL2_EDGES = 1 << 5,  // 32: Add diagonal edges (top-left to bottom-right)
//     FIRST_RING      = 1 << 6,  // 64: Add edges for first radial loop
//     LAST_RING       = 1 << 7,  // 128: Add edges for last radial loop
//     ALTERNATE_DIAG  = 1 << 8,  // 256: Add alternate edges
    
//     // Common combinations
//     BASIC_GRID      = AZIMUTHAL_EDGES | RADIAL_EDGES,
//     FULL_GRID       = BASIC_GRID | DIAGONAL1_EDGES | DIAGONAL2_EDGES,
//     DEFAULT_WIRE    = BASIC_GRID | CLOSED_CIRCUM | FIRST_RING | LAST_RING,
//     STAR            = CLOSED_CIRCUM | RADIAL_EDGES | FIRST_RING,
//     TRIMESH         = BASIC_GRID | DIAGONAL1_EDGES | ALTERNATE_DIAG | FIRST_RING | LAST_RING | CLOSED_CIRCUM
// };


namespace WireFlags{enum{
    NONE            = 0,
    CLOSED_CIRCUM   = 1 << 0,  // 1: Close circumference (periodic in angular direction)
    CAPPED_CENTER   = 1 << 1,  // 2: Add central vertex and connect to first ring
    RADIAL_EDGES    = 1 << 2,  // 4: Add radial edges between rings
    AZIMUTHAL_EDGES = 1 << 3,  // 8: Add edges along each ring
    DIAGONAL1_EDGES = 1 << 4,  // 16: Add diagonal edges (top-right to bottom-left)
    DIAGONAL2_EDGES = 1 << 5,  // 32: Add diagonal edges (top-left to bottom-right)
    FIRST_RING      = 1 << 6,  // 64: Add edges for first radial loop
    LAST_RING       = 1 << 7,  // 128: Add edges for last radial loop
    ALTERNATE_DIAG  = 1 << 8,  // 256: Add alternate edges
    // Common combinations
    BASIC_GRID      = AZIMUTHAL_EDGES | RADIAL_EDGES,
    FULL_GRID       = BASIC_GRID | DIAGONAL1_EDGES | DIAGONAL2_EDGES,
    DEFAULT_WIRE    = BASIC_GRID | CLOSED_CIRCUM | FIRST_RING | LAST_RING,
    STAR            = CLOSED_CIRCUM | RADIAL_EDGES | FIRST_RING,
    TRIMESH         = BASIC_GRID | DIAGONAL1_EDGES | ALTERNATE_DIAG | FIRST_RING | LAST_RING | CLOSED_CIRCUM
};};

#define _unpack_WireFlags(flags) \
    bool bPeriodicB   = (bool)(flags & WireFlags::CLOSED_CIRCUM);    \
    bool bHasCenter   = (bool)(flags & WireFlags::CAPPED_CENTER);    \
    bool bAzimEdges   = (bool)(flags & WireFlags::AZIMUTHAL_EDGES);  \
    bool bRadialEdges = (bool)(flags & WireFlags::RADIAL_EDGES);     \
    bool bDiagEdges1  = (bool)(flags & WireFlags::DIAGONAL1_EDGES);  \
    bool bDiagEdges2  = (bool)(flags & WireFlags::DIAGONAL2_EDGES);  \
    bool bFirstRing   = (bool)(flags & WireFlags::FIRST_RING);     \
    bool bLastRing    = (bool)(flags & WireFlags::LAST_RING);      \
    bool bAlternateDiag= (bool)(flags & WireFlags::ALTERNATE_DIAG);


//inline WireFlags operator|(WireFlags a, WireFlags b) { return static_cast<WireFlags>(static_cast<uint16_t>(a) | static_cast<uint16_t>(b)); }
//inline WireFlags operator&(WireFlags a, WireFlags b) { return static_cast<WireFlags>(static_cast<uint16_t>(a) & static_cast<uint16_t>(b)); }
//inline WireFlags operator^(WireFlags a, WireFlags b) { return static_cast<WireFlags>(static_cast<uint16_t>(a) ^ static_cast<uint16_t>(b)); }


// template<typename UVfunc> Vec3f getUVFuncNormal(Vec2f uv, float h, UVfunc func) { 
//     Vec2f o; Vec3f nor,da,db; 
//     o=uv; o.a+=h; da.set(func(o)); 
//     o=uv; o.a-=h; da.sub(func(o)); 
//     o=uv; o.b+=h; db.set(func(o)); 
//     o=uv; o.b-=h; db.sub(func(o)); 
//     nor.set_cross(db,da); nor.normalize(); 
//     return nor; 
// }







template<typename UVfunc> void UVFunc2smooth(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func) {
    Vec2f duv = UVmax-UVmin; duv.mul({1.0f/n.a,1.0f/n.b});
    std::vector<int> verts((n.a+1)*(n.b+1)); 
    int iv = 0;
    for(int ia=0; ia<=n.a; ia++) { 
        Vec2f uv = {UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia}; 
        for(int ib=0; ib<=n.b; ib++) { 
            Vec3d p,nor; 
            convert(func(uv), p); 
            convert(getUVFuncNormal(uv,0.01,func), nor); 
            verts[iv] = builder.vert(p, nor, Vec2d{uv.x,uv.y}); 
            //printf( " %i %i iv(%i): %i uv(%f,%f) p(%f,%f,%f)\n", ia,ib, iv, verts[iv], uv.a, uv.b,  p.x,p.y,p.z  );
            printf( " %3i %3i verts[ %3i ]: %3i \n", ia,ib, iv, verts[iv] );
            iv++;
            uv.b += duv.b; 
        }
    }
    for(int ia=0; ia<n.a; ia++) {
        for(int ib=0; ib<n.b; ib++) { 
            int i=ia*(n.b+1)+ib; 
            printf( " %3i %3i ivs{ %3i %3i %3i %3i  } \n", ia,ib, i, i+1, i+n.b+1, i+n.b+2 );
            builder.tri(verts[i  ], verts[i+1    ], verts[i+n.b+1]); 
            builder.tri(verts[i+1], verts[i+n.b+2], verts[i+n.b+1]); 
            /// TODO:  Make Faces - we need some switch if this should ple plotted as a face
            //builder.chunk({builder.tris.size()-2*n.a*n.b, 2*n.a*n.b, -1, (int)Builder2::ChunkType::face}); 
        }
        /// TODO: Make triangle strip
        //builder.chunk({builder.tris.size()-2*n.a*n.b, 2*n.a*n.b, -1, (int)Builder2::ChunkType::face});
    }
}

template<typename UVfunc> void UVFunc2wire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func) {
    Vec2f duv = UVmax-UVmin; duv.mul({1.0f/n.a,1.0f/n.b});
    //std::vector<int> verts((n.a+1)*(n.b+1)); int iv = 0;
    int iv = builder.verts.size();
    for(int ia=0; ia<=n.a; ia++) { 
        Vec2f uv = {UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia};
        for(int ib=0; ib<=n.b; ib++) { 
            Vec3d p; 
            convert(func(uv), p); 
            //verts[iv] = 
            builder.vert(p); 
            if(ia<n.a) builder.edge( iv, iv+n.b+1 ); 
            if(ib<n.b) builder.edge( iv, iv+1     ); 
            iv++; 
            uv.b += duv.b; 
            /// TODO:  Make Faces - we need some switch if this should ple plotted as a face
            //builder.chunk({builder.edges.size()-n.a*(n.b+1)-n.b*(n.a+1), n.a*(n.b+1)+n.b*(n.a+1), -1, (int)Builder2::ChunkType::edgestrip});
        }
        /// TODO: Make triangle strip
        //builder.chunk({builder.edges.size()-n.a*(n.b+1)-n.b*(n.a+1), n.a*(n.b+1)+n.b*(n.a+1), -1, (int)Builder2::ChunkType::edgestrip});
    }
    
}

template<typename UVfunc> void UVFunc2wireExtruded(Builder2& builder, Vec2i n, float thick, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func) {
    Vec2f duv = UVmax-UVmin; duv.mul({1.0f/n.a,1.0f/n.b}); float eps = 0.001;
    std::vector<int> verts((n.a+1)*(n.b+1)*2); int iv = 0;
    for(int ia=0; ia<=n.a; ia++) { 
        Vec2f uv = {UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia};
        for(int ib=0; ib<=n.b; ib++) { 
            Vec3d p,nr; 
            convert(func(uv), p); 
            convert(getUVFuncNormal(uv,eps,func), nr);
            verts[iv*2  ] = builder.vert(p); 
            verts[iv*2+1] = builder.vert(p+nr*thick);
            if(ia>0 && ib>0) { 
                int i=iv*2, i_prev=(iv-1)*2, i_up=(iv-(n.b+1))*2;
                builder.tri(verts[i_prev], verts[i  ], verts[i_prev+1]); 
                builder.tri(verts[i     ], verts[i+1], verts[i_prev+1]);
                builder.tri(verts[i_up  ], verts[i  ], verts[i_up  +1]); 
                builder.tri(verts[i     ], verts[i+1], verts[i_up  +1]);
            }
            iv++; 
            uv.b += duv.b;
        }
    }
    //builder.chunk({builder.tris.size()-4*n.a*n.b, 4*n.a*n.b, -1, (int)Builder2::ChunkType::face});
}

// New flexible wireframe functions with center/edge closure control
template<typename UVfunc> void UVFunc2wire_new(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float voff, UVfunc func, int wire_flags=WireFlags::DEFAULT_WIRE ) {
    Vec2f duv = UVmax-UVmin; duv.mul({1.0f/n.a,1.0f/n.b});
    _unpack_WireFlags(wire_flags)
    int  startIndex = builder.verts.size();
    int  centerIndex = -1;

    // Create center vertex if needed
    if(bHasCenter) {
        Vec3d p_center; convert(func({0.0,0.0}), p_center);
        centerIndex = builder.vert(p_center);
    }
    
    // Create grid vertices
    int nCols = bPeriodicB ? n.b : n.b+1;
    for(int ia=(bHasCenter?1:0); ia<=n.a; ia++) {
        Vec2f uv = {UVmin.a+duv.a*ia, UVmin.b+voff*duv.b*ia};
        for(int ib=0; ib<nCols; ib++) {
            Vec3d p; convert(func(uv), p);
            builder.vert(p);
            uv.b += duv.b;
        }
    }
    printf("UVFunc2wire_new(): n.a=%i n.b=%i | bAlternateDiag=%i bDiagEdges1=%i bDiagEdges2=%i bPeriodicB=%i bHasCenter=%i bAzimEdges=%i bRadialEdges=%i bFirstRing=%i bLastRing=%i\n", n.a, n.b, bAlternateDiag, bDiagEdges1, bDiagEdges2, bPeriodicB, bHasCenter, bAzimEdges, bRadialEdges, bFirstRing, bLastRing);
    // Create edges
    for(int ir=0; ir<=n.a-(bHasCenter?1:0); ir++) {
        int rowStart = startIndex + (bHasCenter?1:0) + ir * nCols;
        
        // Connect to center for first ring
        if(bHasCenter && ir==0) {
            for(int ib=0; ib<nCols; ib++) { builder.edge(centerIndex, rowStart + ib); }
        }
        
        // Create edges within current row
        bool bDoAzumith = bAzimEdges;
        if     ( ir==0     ) bDoAzumith = bFirstRing;
        else if( ir==n.a-1 ) bDoAzumith = bLastRing;
        //printf("UVFunc2wire_new()[ir=%i] bDoAzumith=%i | bAzimEdges=%i bFirstRing=%i bLastRing=%i\n", ir, bDoAzumith, bAzimEdges, bFirstRing, bLastRing);
        if(bDoAzumith) {
            //printf("UVFunc2wire_new()[ir=%i] bDoAzumith=%i \n", ir, bDoAzumith);
            for(int ib=0; ib<nCols-1; ib++) { builder.edge(rowStart + ib, rowStart + ib+1); }
            if (bPeriodicB)                 { builder.edge(rowStart + nCols-1, rowStart);    }
        }
        // Create edges between rows
        if(ir > 0 && bRadialEdges) {
            int prevRow = startIndex + (bHasCenter?1:0) + (ir-1) * nCols;
            for(int ib=0; ib<nCols; ib++) { builder.edge(prevRow + ib, rowStart + ib); }
        }
        // Create diagonal edges 
        if(ir > 0 && (bDiagEdges1||bDiagEdges2)) {
            int prevRow = startIndex + (bHasCenter?1:0) + (ir-1) * nCols;

            int nc=nCols-1;
            if(bPeriodicB) nc++;
            for(int ib=0; ib<nc; ib++) {
                int i0=0,i1=1;
                if( bAlternateDiag && ( (bool)((ir^ib)&1) ) ) { i0=1; i1=0; }
                if(ib>=(nCols-1)){ int i0_=i0; i0=i1*(nCols-1)-ib; i1=i0_*(nCols-1)-ib; }
                if(bDiagEdges1) builder.edge(prevRow + ib+i1, rowStart + ib+i0);
                if(bDiagEdges2) builder.edge(prevRow + ib+i0, rowStart + ib+i1);
            }
            // Periodic boundary connections
            // if(bPeriodicB) {
            //     int i0=0,i1 = 1;
            //     if( bAlternateDiag && ( (bool)( (ir^(nCols-1))&1) ) ) { i0 = 1; i1 = 0; }
            //     // // Apply the same alternating logic to boundaries
            //     // if(bAlternateDiag) {
            //     //     bool swap = (ir ^ (nCols-1)) & 1;
            //     //     if(swap) { i0 = 1; i1 = 0; }
            //     // }
            //     // // Create boundary connections using the same pattern
            //     if(bDiagEdges1) builder.edge(prevRow + i0*(nCols-1), rowStart + i1*(nCols-1));
            //     if(bDiagEdges2) builder.edge(prevRow + i1*(nCols-1), rowStart + i0*(nCols-1));
            // }
        }
    }
}


/**
//  * Creates a panel of truss elements between four corner points. Adds the points and edges to the Truss object.
//  * // TODO: make also triangular panel
//  * 
//  * @param p00 The first corner point.
//  * @param p01 The second corner point.
//  * @param p10 The third corner point.
//  * @param p11 The fourth corner point.
//  * @param n The number of subdivisions along each side of the panel.
//  * @param width The width of the truss elements.
//  */
// int Builder2::panel( Vec3d p00, Vec3d p01, Vec3d p10, Vec3d p11, Vec2i n, double width, Quat4i stickTypes ){
//     // ToDo: ad p00,p01,p10,p11 etc. - maybe we should rather specify indexes of existing verts rather than positions of new verts ?
//     //printf( "Mesh::panel() n(%i,%i) w=%g p00(%g,%g,%g) p01(%g,%g,%g) p10(%g,%g,%g) p11(%g,%g,%g) \n", n.x,n.y, p00.x,p00.y,p00.z, p01.x,p01.y,p01.z, p10.x,p10.y,p10.z, p11.x,p11.y,p11.z );
//     //int kind_long   = 0;
//     //int kind_perp   = 1;
//     //int kind_zigIn  = 2;
//     //int kind_zigOut = 3;
//     Vec2d step = {1.0/n.a,1.0/n.b};
//     int di = 2*n.a-1;
//     //int ibloc = block();   // this is better to call manually from outside
//     int i0  = verts.size();
//     //int i00 = verts.size();
//     for (int ib=0; ib<n.b; ib++){
//         double db,mb;
//         db = ib*step.b;     mb=1-db;
//         Vec3d p0  = p00*mb + p10*db;
//         Vec3d p1  = p01*mb + p11*db;
//         db += 0.5*step.b; mb=1-db;
//         Vec3d p0_ = p00*mb + p10*db;
//         Vec3d p1_ = p01*mb + p11*db;
//         for (int ia=0; ia<n.a; ia++){
//             double da,ma;
//             da = ia*step.a; ma = 1-da;
//             Vec3d p   = p0 *ma + p1 *da;
//             //points.push_back( p              );
//             vert( p );
//             int bi = i0+di; if( ib==n.b-2 )bi-=ia;
//             int dia = 2;    if( ib==n.b-1 )dia=1;
//             if (ia<(n.a-1)) edge( i0,i0+dia,stickTypes.y );
//             if (ib<(n.b-1)) edge( i0,bi    ,stickTypes.y );
//             if( (ia<(n.a-1))&&(ib<(n.b-1)) ){ // diagonal
//                 Vec3d p_  = p0_*ma + p1_*da;
//                 da += 0.5*step.a; ma=1-da;
//                 Vec3d p__ = p0_*ma + p1_*da;
//                 Vec3d up; up.set_cross( p_-p, p__-p ); up.normalize();
//                 //points.push_back( p__ + up*width );
//                 vert( p__ + up*width );
//                 if( ia<(n.a-2) ) edge( i0+1,i0+1+dia,stickTypes.z );
//                 if( ib<(n.b-2) ) edge( i0+1,bi+1    ,stickTypes.z );
//                 edge( i0+1,i0     ,stickTypes.w );
//                 edge( i0+1,i0+dia ,stickTypes.w );
//                 edge( i0+1,bi     ,stickTypes.w );
//                 if( ib==n.b-2 )dia=1;
//                 edge( i0+1,bi+dia ,stickTypes.w );
//                 i0++;
//             }
//             i0++;
//         }
//     }
//     return i0;
// }


template<typename UVfunc> 
int UV_panel( Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, double width, Quat4i stickTypes, UVfunc func ){
    // === Build double-layer truss panel mapped by a generic UV function ===
    // number of vertices along U (a) and V (b)
    const int na = n.x;
    const int nb = n.y;

    // step in UV space between base-grid vertices
    Vec2f du = { (UVmax.x-UVmin.x)/(na-1), 0.0f };
    Vec2f dv = { 0.0f, (UVmax.y-UVmin.y)/(nb-1) };

    const int base0 = builder.verts.size();

    // --- 1) base grid vertices ---
    std::vector<int> baseIdx(na*nb);
    for(int ib=0; ib<nb; ++ib){
        for(int ia=0; ia<na; ++ia){
            Vec2f uv = { UVmin.x + ia*du.x, UVmin.y + ib*dv.y };
            int idx  = builder.verts.size();
            baseIdx[ib*na+ia] = idx;
            builder.vert( (Vec3d)func(uv) );
        }
    }

    // --- 2) center layer (raised) vertices ---
    std::vector<int> cenIdx( (na-1)*(nb-1) );
    for(int ib=0; ib<nb-1; ++ib){
        for(int ia=0; ia<na-1; ++ia){
            // Fetch the four corner positions of the cell
            Vec2f uv00 = { UVmin.x + ia    *du.x, UVmin.y + ib    *dv.y };
            Vec2f uv10 = { UVmin.x + (ia+1)*du.x, UVmin.y + ib      *dv.y };
            Vec2f uv01 = { UVmin.x + ia      *du.x, UVmin.y + (ib+1)*dv.y };
            // sample positions
            Vec3d p00 = (Vec3d)func(uv00);
            Vec3d p10 = (Vec3d)func(uv10);
            Vec3d p01 = (Vec3d)func(uv01);
            // bilinear center (approx)
            Vec2f uvC = { uv00.x + 0.5f*du.x, uv00.y + 0.5f*dv.y };
            Vec3d  pc = (Vec3d)func( uvC );
            // estimate normal from the two edge vectors in the cell
            Vec3d uvec = p10 - p00;
            Vec3d vvec = p01 - p00;
            Vec3d  n; n.set_cross( uvec, vvec ); n.normalize();
            pc.add_mul( n, width );
            int idx = builder.verts.size();
            cenIdx[ ib*(na-1) + ia ] = idx;
            builder.vert( pc );
        }
    }

    // --- 3) base grid edges (longitudinal & transversal) ---
    for(int ib=0; ib<nb; ++ib){
        for(int ia=0; ia<na; ++ia){
            int v = baseIdx[ib*na+ia];
            if( ia<na-1 ) builder.edge( v, baseIdx[ib*na+ia+1], stickTypes.y );
            if( ib<nb-1 ) builder.edge( v, baseIdx[(ib+1)*na+ia], stickTypes.y );
        }
    }

    // --- 4) connect center vertices to base vertices & diagonal stiffeners ---
    for(int ib=0; ib<nb-1; ++ib){
        for(int ia=0; ia<na-1; ++ia){
            int c  = cenIdx[ ib*(na-1)+ia ];
            int v00 = baseIdx[ ib   *na + ia   ];
            int v10 = baseIdx[ ib   *na + ia+1 ];
            int v01 = baseIdx[(ib+1)*na + ia   ];
            int v11 = baseIdx[(ib+1)*na + ia+1 ];

            // spokes from center to the four corners
            builder.edge( c, v00, stickTypes.w );
            builder.edge( c, v10, stickTypes.w );
            builder.edge( c, v01, stickTypes.w );
            builder.edge( c, v11, stickTypes.w );

            // diagonal connections between centers for zig-zag stiffening
            if( ia<na-2 ){                // horizontal diag to next center
                builder.edge( c, cenIdx[ ib*(na-1)+ia+1 ], stickTypes.z );
            }
            if( ib<nb-2 ){                // vertical diag to next row center
                builder.edge( c, cenIdx[ (ib+1)*(na-1)+ia ], stickTypes.z );
            }
        }
    }

    return base0;
}


/*
UV_slab
axial[3]:
(1,0,0)
(0,1,0)
(0,0,1)
face-diagonal[6]:
(1,1,0)
(1,0,1)
(0,1,1)
( 1,-1, 0)
( 0, 1,-1)
(-1, 0, 1)
space-diagonal[8]:
( 1, 1, 1)
( 1, 1,-1)
( 1,-1, 1)
( 1,-1,-1)
(-1, 1, 1)
(-1, 1,-1)
(-1,-1, 1)
(-1,-1,-1)
*/



// === General cubic-grid based slab builder ===
// dirMask : bitmask selecting which edge directions to include using following order of directions
// 0:(1,0,0)  1:(0,1,0)  2:(0,0,1)
// 3:(1,1,0)  4:(1,0,1)  5:(0,1,1)
// 6:(1,-1,0) 7:(0,1,-1) 8:(-1,0,1)
// 9:(1,1,1) 10:(1,1,-1) 11:(1,-1,1) 12:(1,-1,-1)
// stickTypes : mapping of edge "kind" to Builder2 edge-type slot
//    x : vertical (z-axis)
//    y : axial planar (x or y)
//    z : planar diagonals (x-y plane)
//    w : any edge touching both layers (face or space diagonals)
template<typename UVfunc>
int UV_slab( Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec3f up, int dirMask, Quat4i stickTypes, UVfunc func ){
    const int na = n.x;
    const int nb = n.y;
    const int nz = 2;                 // just bottom & top layer

    // UV spacing
    Vec2f du = { (UVmax.x-UVmin.x)/(na-1), 0.0f };
    Vec2f dv = { 0.0f, (UVmax.y-UVmin.y)/(nb-1) };

    const int base0 = builder.verts.size();

    // Allocate vertex index grids for both layers
    std::vector<int> idx0(na*nb);
    std::vector<int> idx1(na*nb);

    // Build vertices
    for(int iy=0; iy<nb; ++iy){
        for(int ix=0; ix<na; ++ix){
            Vec2f uv = { UVmin.x + ix*du.x, UVmin.y + iy*dv.y };
            Vec3d p  = (Vec3d)func( uv );

            Vec3d u  = (Vec3d)func( uv + du )-p; //u.normalize();
            Vec3d v  = (Vec3d)func( uv + dv )-p; //v.normalize();
            Vec3d nor; nor.set_cross(u,v); nor.normalize();

            int iFlat = iy*na + ix;
            idx0[iFlat] = builder.vert( p );                       // bottom
            //idx1[iFlat] = builder.vert( Vec3d{ p.x, p.y, p.z + width } ); // top (straight up)

            Vec3d p1 = p + nor*up.z + u*up.x + v*up.y;
            idx1[iFlat] = builder.vert( p1 ); // top (straight up)
        }
    }

    // Unique positive directions (13)
    static const Vec3i DIRS[13] = {
        // axial
        { 1, 0, 0 }, // 0 x
        { 0, 1, 0 }, // 1 y
        { 0, 0, 1 }, // 2 z
        // face diag
        { 1, 1, 0 }, // 3 x+y
        { 1,-1, 0 }, // 4 x-y
        { 1, 0, 1 }, // 5 x+z
        {-1, 0, 1 }, // 6 z-x
        { 0, 1, 1 }, // 7 y+z
        { 0, 1,-1 }, // 8 y-z
        // space diag
        { 1, 1, 1 }, // 9 x+y+z
        { 1, 1,-1 }, // 10 x+y-z
        { 1,-1, 1 }, // 11 x-y+z
        { 1,-1,-1 }  // 12 x-y-z
    };

    auto getVert = [&](int ix,int iy,int iz){ int i = iy*na + ix; return (iz==0)? idx0[i] : idx1[i]; };

    auto edgeTypeByDir = [&](const Vec3i& d){
        int comps = (d.x!=0) + (d.y!=0) + (d.z!=0);
        if(comps==1){ return d.z ? stickTypes.x : stickTypes.y; }         // axis
        if(comps==2){ return d.z ? stickTypes.w : stickTypes.z; }         // face diag
        return stickTypes.w;                                              // space diag
    };

    // Build edges according to mask
    for(int iz=0; iz<nz; ++iz){
        for(int iy=0; iy<nb; ++iy){
            for(int ix=0; ix<na; ++ix){
                for(int id=0; id<13; ++id){
                    if(!(dirMask & (1<<id))) continue;
                    const Vec3i& d = DIRS[id];
                    int jx = ix + d.x;
                    int jy = iy + d.y;
                    int jz = iz + d.z;
                    if(jx<0||jx>=na||jy<0||jy>=nb||jz<0||jz>=nz) continue; // stay inside slab
                    int a = getVert(ix,iy,iz);
                    int b = getVert(jx,jy,jz);
                    builder.edge( a, b, edgeTypeByDir(d) );
                }
            }
        }
    }

    return base0;
}

void QuadSlab(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec3f p00, Vec3f p01, Vec3f p10, Vec3f p11, Vec3f up, int dirMask, Quat4i stickTypes ) {
    auto uvfunc = [&](Vec2f uv){ return QuadUVfunc(uv,p00,p01,p10,p11); };
    UV_slab(builder,n,UVmin,UVmax,up,dirMask,stickTypes,uvfunc);
}

void QuadPanel(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec3f p00, Vec3f p01, Vec3f p10, Vec3f p11, float width, Quat4i stickTypes ) {
    auto uvfunc = [&](Vec2f uv){ return QuadUVfunc(uv,p00,p01,p10,p11); };
    UV_panel(builder,n,UVmin,UVmax,width,stickTypes,uvfunc);
}

void Tube(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec2f Rs, float L, float width, Quat4i stickTypes ) {
    auto uvfunc = [&](Vec2f uv){ return ConeUVfunc(uv,Rs.a,Rs.b,L); };
    UV_panel(builder,n,UVmin,UVmax,width,stickTypes,uvfunc);
}

void Parabola_Wire_new(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, int wire_flags=WireFlags::DEFAULT_WIRE) {
    float K = L/(R*R);
    UVmin.a *= R; UVmax.a *= R;
    auto uvfunc = [K](Vec2f uv) { return ParabolaUVfunc(uv, K); };
    UVFunc2wire_new(builder, n, UVmin, UVmax, voff, uvfunc, wire_flags);
}

void Cone2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Sphere2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Torus2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Teardrop2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void NACASegment2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float *coefs1, float *coefs2, float L, float voff, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return NACA4digitUVfunc(uv,coefs1,coefs2,L);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void HarmonicTube2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, bool wire) { 
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Parabola2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, bool wire) { 
    float K = L/(R*R); 
    UVmin.a*=R; UVmax.a*=R; 
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);}; 
    if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
    else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Hyperbola2Mesh(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, bool wire) {
    if(r>0){ 
        float K = R/L; 
        UVmin.a*=L; UVmax.a*=L; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaLUVfunc(uv,r,K);}; 
        if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
        else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
    }else{ 
        r=-r; 
        float K = L/R; 
        UVmin.a*=R; UVmax.a*=R; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaRUVfunc(uv,r,K);}; 
        if(wire) UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
        else     UVFunc2smooth(builder, n,UVmin,UVmax,voff,uvfunc); 
    }
}

void Cone_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return ConeUVfunc(uv,R1,R2,L);}; 
    UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Sphere_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return SphereUVfunc(uv,R);}; 
    UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Torus_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return TorusUVfunc(uv,r,R);}; 
    UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Teardrop_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return TeardropUVfunc(uv,R1,R2,L);}; 
    UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void NACASegment_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float *coefs1, float *coefs2, float L, float voff, float thick) { 
    auto uvfunc = [&](Vec2f uv){return NACA4digitUVfunc(uv,coefs1,coefs2,L);}; 
    UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void HarmonicTube_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R1, float R2, float L, float voff, float freq, float amp, float thick) { 
    auto uvfunc = [&](Vec2f uv){return HarmonicTubeUVfunc(uv,R1,R2,L,freq,amp);}; 
    UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Parabola_Wire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff) { 
    float K = L/(R*R); UVmin.a*=R; UVmax.a*=R; 
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);}; 
    UVFunc2wire(builder, n,UVmin,UVmax,voff,uvfunc); 
}

void Parabola_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float R, float L, float voff, float thick) { 
    float K = L/(R*R); UVmin.a*=R; UVmax.a*=R; 
    auto uvfunc = [&](Vec2f uv){return ParabolaUVfunc(uv,K);}; 
    UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
}

void Hyperbola_ExtrudedWire(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, float r, float R, float L, float voff, float thick) {
    if(r>0){ 
        float K = R/L; UVmin.a*=L; UVmax.a*=L; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaLUVfunc(uv,r,K);}; 
        UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
    }else{ 
        r=-r; float K = L/R; UVmin.a*=R; UVmax.a*=R; 
        auto uvfunc = [&](Vec2f uv){return HyperbolaRUVfunc(uv,r,K);}; 
        UVFunc2wireExtruded(builder, n,thick,UVmin,UVmax,voff,uvfunc); 
    }
}

void ConeFan(Builder2& builder, int n, float r, const Vec3f& base, const Vec3f& tip) {
    Vec3f a,b,c,c_hat; c.set_sub(tip,base); c_hat.set_mul(c,1/c.norm()); c_hat.getSomeOrtho(a,b); a.normalize(); b.normalize();
    float alfa=2*M_PI/n; Vec2f rot; rot.set(1.0f,0.0f); Vec2f drot; drot.set(cos(alfa),sin(alfa));
    Vec3f q=c; q.add_mul(a,-r); float pnab=c_hat.dot(q)/q.norm(), pnc=sqrt(1-pnab*pnab);
    int i0=builder.verts.size(); builder.vert(Vec3d{tip.x,tip.y,tip.z});
    for(int i=0; i<=n; i++) {
        Vec3f p; p.set(rot.x*a.x+rot.y*b.x, rot.x*a.y+rot.y*b.y, rot.x*a.z+rot.y*b.z);
        Vec3f pn; pn.set(pnab*p.x+pnc*c_hat.x, pnab*p.y+pnc*c_hat.y, pnab*p.z+pnc*c_hat.z);
        builder.vert(Vec3d{base.x+r*p.x,base.y+r*p.y,base.z+r*p.z}, Vec3d{pn.x,pn.y,pn.z});
        if(i>0) builder.tri(i0,i0+i,i0+i+1);
        rot.mul_cmplx(drot);
    }
    builder.chunk({builder.tris.size()-n, n, -1, (int)Builder2::ChunkType::face});
}

void CylinderStrip(Builder2& builder, int n, float r1, float r2, const Vec3f& base, const Vec3f& tip) {
    Vec3f a,b,c,c_hat; c.set_sub(tip,base); c_hat.set_mul(c,1/c.norm()); c_hat.getSomeOrtho(a,b); a.normalize(); b.normalize();
    float alfa=2*M_PI/n; Vec2f rot; rot.set(1.0f,0.0f); Vec2f drot; drot.set(cos(alfa),sin(alfa));
    Vec3f q=c; q.add_mul(a,-(r1-r2)); float pnab=c_hat.dot(q)/q.norm(), pnc=sqrt(1-pnab*pnab);
    int i0=builder.verts.size();
    for(int i=0; i<=n; i++) {
        Vec3f p; p.set(rot.x*a.x+rot.y*b.x, rot.x*a.y+rot.y*b.y, rot.x*a.z+rot.y*b.z);
        Vec3f pn; pn.set(pnab*p.x+pnc*c_hat.x, pnab*p.y+pnc*c_hat.y, pnab*p.z+pnc*c_hat.z);
        builder.vert(Vec3d{base.x+r1*p.x,base.y+r1*p.y,base.z+r1*p.z}, Vec3d{pn.x,pn.y,pn.z});
        builder.vert(Vec3d{tip.x+r2*p.x,tip.y+r2*p.y,tip.z+r2*p.z}, Vec3d{pn.x,pn.y,pn.z});
        if(i>0) { builder.tri(i0+(i-1)*2,i0+i*2,i0+(i-1)*2+1); builder.tri(i0+i*2,i0+i*2+1,i0+(i-1)*2+1); }
        rot.mul_cmplx(drot);
    }
    builder.chunk({builder.tris.size()-2*n, 2*n, -1, (int)Builder2::ChunkType::face});
}

void CylinderStrip_wire(Builder2& builder, int n, float r1, float r2, const Vec3f& base, const Vec3f& tip) {
    Vec3f a,b,c,c_hat; c.set_sub(tip,base); c_hat.set_mul(c,1/c.norm()); c_hat.getSomeOrtho(a,b); a.normalize(); b.normalize();
    float alfa=2*M_PI/n; Vec2f rot; rot.set(1.0f,0.0f); Vec2f drot; drot.set(cos(alfa),sin(alfa));
    int i0=builder.verts.size();
    for(int i=0; i<n; i++) {
        Vec3f p; p.set(rot.x*a.x+rot.y*b.x, rot.x*a.y+rot.y*b.y, rot.x*a.z+rot.y*b.z);
        builder.vert(Vec3d{base.x+r1*p.x,base.y+r1*p.y,base.z+r1*p.z});
        builder.vert(Vec3d{tip.x+r2*p.x,tip.y+r2*p.y,tip.z+r2*p.z});
        if(i>0) { builder.edge(i0+(i-1)*2,i0+i*2); builder.edge(i0+(i-1)*2+1,i0+i*2+1); }
        builder.edge(i0+i*2,i0+i*2+1);
        rot.mul_cmplx(drot);
    }
    builder.edge(i0+(n-1)*2,i0); builder.edge(i0+(n-1)*2+1,i0+1);
    builder.chunk({builder.edges.size()-3*n, 3*n, -1, (int)Builder2::ChunkType::edgestrip});
}

void SphereTriangle_wire(Builder2& builder, int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c) {
    float d=1.0f/n; Vec3f da(a-c),db(b-c); da.mul(d); db.mul(d);
    int i0=builder.verts.size();
    for(int ia=0; ia<=n; ia++) {
        Vec3f p0=c; p0.add_mul(da,ia);
        for(int ib=0; ib<=n-ia; ib++) {
            Vec3f p=p0; p.add_mul(db,ib); p.mul(r/p.norm());
            builder.vert(Vec3d{pos.x+p.x,pos.y+p.y,pos.z+p.z});
            if(ia<n && ib<n-ia) {
                int i1=i0+ia*(2*n-ia+1)/2+ib, i2=i0+(ia+1)*(2*n-ia)/2+ib;
                builder.edge(i1,i2); if(ib<n-ia-1) builder.edge(i2,i2+1);
            }
        }
    }
    builder.chunk({builder.edges.size()-2*n*n, 2*n*n, -1, (int)Builder2::ChunkType::edgestrip});
}

void SphereTriangle(Builder2& builder, int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c) {
    float d=1.0f/n; Vec3f da(a-c),db(b-c); da.mul(d); db.mul(d);
    int i0=builder.verts.size();
    for(int ia=0; ia<=n; ia++) {
        Vec3f p0=c; p0.add_mul(da,ia);
        for(int ib=0; ib<=n-ia; ib++) {
            Vec3f p=p0; p.add_mul(db,ib); p.mul(r/p.norm());
            builder.vert(Vec3d{pos.x+p.x,pos.y+p.y,pos.z+p.z}, Vec3d{p.x/r,p.y/r,p.z/r});
            if(ia<n && ib<n-ia) {
                int i1=i0+ia*(2*n-ia+1)/2+ib, i2=i0+(ia+1)*(2*n-ia)/2+ib;
                builder.tri(i1,i2,i1+1); if(ib<n-ia-1) builder.tri(i2,i2+1,i1+1);
            }
        }
    }
    builder.chunk({builder.tris.size()-n*n, n*n, -1, (int)Builder2::ChunkType::face});
}

void Sphere_oct(Builder2& builder, int n, float r, const Vec3f& pos, bool wire) {
    Vec3f a,b,c; 
    a.set(1.0f,0.0f,0.0f); b.set(0.0f,1.0f,0.0f); c.set(0.0f,0.0f,1.0f);
    Vec3f na=a; na.mul(-1.0f);
    Vec3f nb=b; nb.mul(-1.0f);
    Vec3f nc=c; nc.mul(-1.0f);
    if(wire) {
        SphereTriangle_wire(builder,n,r,pos, a, b, c); SphereTriangle_wire(builder,n,r,pos,na, b, c);
        SphereTriangle_wire(builder,n,r,pos, a,nb, c); SphereTriangle_wire(builder,n,r,pos,na,nb, c);
        SphereTriangle_wire(builder,n,r,pos, a, b,nc); SphereTriangle_wire(builder,n,r,pos,na, b,nc);
        SphereTriangle_wire(builder,n,r,pos, a,nb,nc); SphereTriangle_wire(builder,n,r,pos,na,nb,nc);
    } else {
        SphereTriangle(builder,n,r,pos, a, b, c); SphereTriangle(builder,n,r,pos,na, b, c);
        SphereTriangle(builder,n,r,pos, a,nb, c); SphereTriangle(builder,n,r,pos,na,nb, c);
        SphereTriangle(builder,n,r,pos, a, b,nc); SphereTriangle(builder,n,r,pos,na, b,nc);
        SphereTriangle(builder,n,r,pos, a,nb,nc); SphereTriangle(builder,n,r,pos,na,nb,nc);
    }
}

} // namespace Mesh

#endif
