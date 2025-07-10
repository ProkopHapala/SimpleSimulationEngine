#ifndef DrawUV_h
#define DrawUV_h

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "MeshBuilder2.h"
#include "geometry/UVfuncs.h"


// ToDo: 
// 1) UVfunc should return Vec3d 
// 2) UVfunc should be able to return normal
// 3) UVfunc should have input Vec3d instead of Vec2f, 3-components (U,V,h) is along normal direction to the surface


namespace Mesh {

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
    bool bAlternateDiag= (bool)(flags & WireFlags::ALTERNATE_DIAG)

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
    _unpack_WireFlags(wire_flags);
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




/* --- Original dense UV_sheet kept unchanged above --- */

// New variant: generate single-layer sheet with diagonal clipping (ix+iy > imin && ix+iy < imax)
// Stores vertices row-wise; idx0[iy] contains global index of first vertex in row iy
// Neighbor lookup is O(1) using per-row first index and ix0s[ix]
// NOTE: imin/imax are strict inequalities in the condition above
//       Set imin=-1, imax=1e9 for no clipping (full grid)

template<typename UVfunc>
int UV_sheet_clip_verts( Builder2& builder, Vec2i n, Vec2f uv0, Vec2f duv, std::vector<Vec3i>& idx, UVfunc func, int imin, int imax ){
    const int iv0 = builder.verts.size();
    idx.resize(n.y);
    
    for(int iy=0; iy<n.y; ++iy){
        int ix0 = std::max( 0,     imin+1 - iy );
        int ix1 = std::min( n.x-1, imax-1 - iy );
        if(ix0>ix1){ idx[iy] = {-1,-1,-1}; continue; }
        
        idx[iy] = {builder.verts.size(), ix0, ix1};
        for(int ix=ix0; ix<=ix1; ++ix){
            Vec2f uv = { uv0.x + ix*duv.x, uv0.y + iy*duv.y };
            Vec3d p  = (Vec3d)func( uv );
            builder.vert(p);
        }
    }
    return iv0;
}

int UV_sheet_clip_edges( Builder2& builder, Vec2i n, const std::vector<Vec3i>& idx, int dirMask, Quat4i stickTypes ){
    static const Vec2i DIRS[4] = { {1,0}, {0,1}, {1,1}, {1,-1} };
    const int ie0 = builder.edges.size();
    
    for(int iy=0; iy<n.y; ++iy){
        const Vec3i& row = idx[iy];
        if(row.x<0) continue;
        for(int ix=row.y; ix<=row.z; ++ix){
            int a = row.x + (ix - row.y);
            for(int id=0; id<4; ++id){
                if(!(dirMask & (1<<id))) continue;
                const Vec2i& d = DIRS[id];
                int jx = ix + d.x;
                int jy = iy + d.y;
                if(jy<0 || jy>=n.y) continue;
                const Vec3i& jrow = idx[jy];
                if(jrow.x<0) continue;
                if(jx<jrow.y || jx>jrow.z) continue;
                int b = jrow.x + (jx - jrow.y);
                int edgeType = (d.x && d.y) ? stickTypes.z : stickTypes.y; // diagonal or axial
                builder.edge(a,b,edgeType);
            }
        }
    }
    return ie0;
}

template<typename UVfunc>
int UV_sheet_clip( Builder2& builder, Vec2i n, Vec2f uv0, Vec2f duv, int dirMask, Quat4i stickTypes, UVfunc func, int imin=0, int imax=100 ){
    std::vector<Vec3i> idx;
    int iv0 = UV_sheet_clip_verts(builder,n,uv0,duv,idx,func,imin,imax);
    UV_sheet_clip_edges(builder,n,idx,dirMask,stickTypes);
    return iv0;
}

// Single-layer 2D neighborhood with 4 sticks (axial and diagonal)
int stickEdges2D( Builder2& builder, Vec2i n, int* idx, int dirMask, Quat4i stickTypes ){
    // 4 directions in 2D (axial and diagonal)
    static const Vec2i DIRS[4] = {
        {1,0},  // 0 x
        {0,1},  // 1 y
        {1,1},  // 2 x+y
        {1,-1}  // 3 x-y
    };
    //auto getVert       = [&](int ix,int iy ){ return idx[iy*n.x + ix]; };
    //auto edgeTypeByDir = [&](const Vec2i& d){ return (d.x && d.y) ? stickTypes.z : stickTypes.y; }; // diagonal or axial
    const int ie0 = builder.edges.size();
    for(int iy=0; iy<n.y; ++iy){
        for(int ix=0; ix<n.x; ++ix){
            for(int id=0; id<4; ++id){
                if(!(dirMask & (1<<id))) continue;
                const Vec2i& d = DIRS[id];
                int jx = ix + d.x;
                int jy = iy + d.y;
                if(jx<0||jx>=n.x||jy<0||jy>=n.y) continue;
                //int a = getVert(ix,iy);
                //int b = getVert(jx,jy);
                int a = idx[iy*n.x + ix];
                int b = idx[jy*n.x + jx];
                //builder.edge(a, b, edgeTypeByDir(d));
                builder.edge(a, b, (d.x && d.y) ? stickTypes.z : stickTypes.y  );
            }
        }
    }
    return ie0;
}

template<typename UVfunc>
int UV_sheet( Builder2& builder, Vec2i n, Vec2f uv0, Vec2f duv, int dirMask, Quat4i stickTypes, UVfunc func, int imin=0, int imax=100 ){
    std::vector<int> idx(n.x*n.y);
    int iv0 = UV_slab_verts(builder,n,uv0,duv,idx.data(),func);
    stickEdges2D(builder,n,idx.data(),dirMask,stickTypes);
    return iv0;
}

void QuadSheet(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec3f p00, Vec3f p01, Vec3f p10, Vec3f p11, int dirMask, Quat4i stickTypes, int imin=0, int imax=100 ) {
    auto uvfunc = [&](Vec2f uv){ return QuadUVfunc(uv,p00,p01,p10,p11); };
    Vec2f duv = { (UVmax.x-UVmin.x)/(n.x-1), (UVmax.y-UVmin.y)/(n.y-1) };
    //UV_sheet(builder,n,UVmin,duv,dirMask,stickTypes,uvfunc,imin,imax);
    UV_sheet_clip(builder,n,UVmin,duv,dirMask,stickTypes,uvfunc,imin,imax);
}

void TubeSheet(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec2f Rs, float L, int dirMask, Quat4i stickTypes ) {
    float dudv = 0.5*(n.x-1.)/(n.y-1.);
    auto uvfunc = [&](Vec2f uv){ uv.y+=uv.x*dudv; uv.y*=2*M_PI; return ConeUVfunc(uv,Rs.a,Rs.b,L); };
    Vec2f duv = { (UVmax.x-UVmin.x)/(n.x-1), (UVmax.y-UVmin.y)/(n.y-1) };
    UV_sheet(builder,n,UVmin,duv,dirMask,stickTypes,uvfunc);
} 






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

int slabEdges( Builder2& builder, Vec2i n, int* idx0, int* idx1, int dirMask, Quat4i stickTypes ){
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
    const int nz = 2;  
    auto getVert = [&](int ix,int iy,int iz){ int i = iy*n.x + ix; return (iz==0)? idx0[i] : idx1[i]; };
    auto edgeTypeByDir = [&](const Vec3i& d){
        int comps = (d.x!=0) + (d.y!=0) + (d.z!=0);
        if(comps==1){ return d.z ? stickTypes.x : stickTypes.y; }         // axis
        if(comps==2){ return d.z ? stickTypes.w : stickTypes.z; }         // face diag
        return stickTypes.w;                                              // space diag
    };
    // just bottom & top layer
    // Build edges according to mask
    const int ie0 = builder.edges.size();
    for(int iz=0; iz<nz; ++iz){
        for(int iy=0; iy<n.y; ++iy){
            for(int ix=0; ix<n.x; ++ix){
                for(int id=0; id<13; ++id){
                    if(!(dirMask & (1<<id))) continue;
                    const Vec3i& d = DIRS[id];
                    int jx = ix + d.x;
                    int jy = iy + d.y;
                    int jz = iz + d.z;
                    if(jx<0||jx>=n.x||jy<0||jy>=n.y||jz<0||jz>=nz) continue; // stay inside slab
                    int a = getVert(ix,iy,iz);
                    int b = getVert(jx,jy,jz);
                    builder.edge( a, b, edgeTypeByDir(d) );
                }
            }
        }
    }
    return ie0;
}

template<typename UVfunc>
int UV_slab_verts( Builder2& builder, Vec2i n, Vec2f uv0, Vec2f duv, int* idx, UVfunc func ){
    const int iv0 = builder.verts.size();
    for(int iy=0; iy<n.y; ++iy){
        for(int ix=0; ix<n.x; ++ix){
            Vec2f uv = { uv0.x + ix*duv.x, uv0.y + iy*duv.y };
            Vec3d p  = (Vec3d)func( uv );
            int iFlat = iy*n.x + ix;
            idx[iFlat] = builder.vert( p );                       // bottom
        }
    }
    return iv0;
}

template<typename UVfunc1, typename UVfunc2>
int UV_slab( Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec3f up, int dirMask, Quat4i stickTypes, UVfunc1 func1, UVfunc2 func2 ){
    Vec2f duv = { (UVmax.x-UVmin.x)/(n.x-1),  
                  (UVmax.y-UVmin.y)/(n.y-1) };
    std::vector<int> idx0(n.x*n.y); // perhaps we do not need this, if we have iv0, iv1
    std::vector<int> idx1(n.x*n.y);
    int iv0 = UV_slab_verts(builder,n,UVmin            ,duv,idx0.data(),func1);
    int iv1 = UV_slab_verts(builder,n,UVmin+duv*up.xy(),duv,idx1.data(),func2);
    slabEdges( builder, n, idx0.data(), idx1.data(), dirMask, stickTypes );
    return iv0;
}

void QuadSlab(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec3f p00, Vec3f p01, Vec3f p10, Vec3f p11, Vec3f up, int dirMask, Quat4i stickTypes ) {
    Vec3f nor; nor.set_cross(p10-p00,p01-p00); nor.normalize();
    auto uvfunc1 = [&](Vec2f uv){ return QuadUVfunc(uv,p00,p01,p10,p11); };
    auto uvfunc2 = [&](Vec2f uv){ return QuadUVfunc(uv,p00,p01,p10,p11) + nor*up.z; };
    UV_slab(builder,n,UVmin,UVmax,up,dirMask,stickTypes,uvfunc1,uvfunc2);
}

void SlabTube(Builder2& builder, Vec2i n, Vec2f UVmin, Vec2f UVmax, Vec2f Rs, float L, Vec3f up, int dirMask, Quat4i stickTypes ) {
    //auto uvfunc = [&](Vec2f uv){ Vec2f uv_= {uv.x+uv.y*0.5, uv.y*0.86602540378}; return ConeUVfunc(uv_,Rs.a,Rs.b,L); };
    //auto uvfunc = [&](Vec2f uv){ Vec2f uv_= {uv.x*0.86602540378, uv.y+uv.x*0.5}; return ConeUVfunc(uv_,Rs.a,Rs.b,L); };
    float dudv = 0.5*(n.x-1.)/(n.y-1.);
    auto uvfunc1 = [&](Vec2f uv){ uv.y+=uv.x*dudv; uv.y*=2*M_PI; return ConeUVfunc(uv,Rs.a     ,Rs.b     ,L); };
    auto uvfunc2 = [&](Vec2f uv){ uv.y+=uv.x*dudv; uv.y*=2*M_PI; return ConeUVfunc(uv,Rs.a+up.z,Rs.b+up.z,L); };
    UV_slab(builder,n,UVmin,UVmax,up,dirMask,stickTypes,uvfunc1,uvfunc2);
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
