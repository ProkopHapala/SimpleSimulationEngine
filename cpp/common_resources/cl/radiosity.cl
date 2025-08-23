


// ====================== Radiosity and ray-tracing ======================

float3 getSomeUp( const float3 v ){
	if( v.x<v.y){ return (float3){ -v.y*v.y -v.z*v.z, v.x*v.y           , v.x*v.z  }; }
    else        { return (float3){  v.y*v.x         , -v.z*v.z -v.x*v.x ,  v.y*v.z }; }
}

bool originInTriangle( const float2 a, const float2 b, const float2 c ){
    float   sgn = a.x*(b.y-a.y) - a.y*(b.x-a.x);
	if( 0 > sgn*( b.x*(c.y-b.y) - b.y*(c.x-b.x) ) ) return false;
	if( 0 > sgn*( c.x*(a.y-c.y) - c.y*(a.x-c.x) ) ) return false;
    //printf("passed\n");
    return true;
}

bool rayInTriangle( const float3 a_, const float3 b_, const float3 c_, const float3 hX, const float3 hY ){
	float2 a = (float2){ dot(hX,a_), dot(hY,a_) };
    float2 b = (float2){ dot(hX,b_), dot(hY,b_) };
    float2 c = (float2){ dot(hX,c_), dot(hY,c_) };
    return originInTriangle( a, b, c );
}

float sdOriginLine(  const float2 p1, const float2 d ){
    const float2 T   = (float2){ d.y, -d.x }; // normal
    const float  det = dot(p1,T);
    //if(det < 0.0f){ return det; } // Optimization: If we are inside, we don't care how deep
    return det/sqrt( dot(d,d) );
}

float sdOriginTriangle( const float2 a, const float2 b, const float2 c ){
    float r;
    float2 ba = b-a;
    float2 cb = c-b;
    float sgn = ( dot(cb, (float2){-ba.y,ba.x} )>0.f )? 1.f : -1.f;
    float rab = sdOriginLine( a, ba  )*sgn; r=rab;
    float rbc = sdOriginLine( b, cb  )*sgn; r=min(r,rbc);
    float rca = sdOriginLine( c, a-c )*sgn; r=min(r,rca);
    return r;
}

float sdRayTrinagle( const float3 a_, const float3 b_, const float3 c_, const float3 hX, const float3 hY ){
    return sdOriginTriangle( 
        (float2){ dot(hX,a_), dot(hY,a_) },
        (float2){ dot(hX,b_), dot(hY,b_) },
        (float2){ dot(hX,c_), dot(hY,c_) } 
    );
}


float diffusionLightCoupling(float3 h, float r, float4 s1, float4 s2 ){
    //float3  d = elj.pos - eli.pos;
    //double r = d.normalize();
    //float r2 = dot(d,d);
    float r2 = r*r;
    float coupling = dot(h,s1.xyz)*dot(h,s2.xyz)/r2;
    //if( fabs(coupling) < couplingTrashold ){ M[i*n+j]=0.0; continue; };  // to make the matrix more sparse
    coupling /= ( r2 + s1.w + s2.w );
    //coupling *=  eli.area * elj.area;
    coupling *= 2* s1.w/(4*M_PI);
}

/*

IDEA:
We should run Occlusion Ray Tracing in two stages 
using Groups of points (i.e. (1) chunks of points bouned within a sphere) and groups of obstacles (i.e. chunks of triangles bounded within larger triangle )
1) 1st pass check approximative visibility between groups of points, occluded by groups of obstacles (i.e. triangles)
   - For each pair of point-groups i,j the approximative check can result in one of three cases:
      1. No occlusion      -  all obstacles are in save distace (R_safe) from the ray connecting the tow point groups, M[i,j]=1.0
      2. full occlusion    - some obstacle fully occludes the ray between the two point groups, M[i,j]=0.0 
      3. partial occlusion - the cloeses obstacles are in such distace from the ray that approximate bounding-volues (Spheres, trinagle) are not enough to determine the occlusion of individual points.
        - We store 0.0< M[i,j]<1.0
        - We store the indexes of all obstacles that are closer than save-distance (R_safe) from the ray
2) 2nd pass - We process only the partial occlusions (i.e. 0.0< M[i,j]<1.0) and check the exact occlusion of individual points within point-group pairs and tringles within the obstacle-groups.
   - we only iterate over the obstacle groups listed in the 1st pass
*/



// Kernel: makeRadiosityCouplings
// Computes diffuse coupling between element centers with triangle occluders
// (binary blocking via projection test). Used as a simple radiosity prototype.
__kernel void  makeRadiosityCouplings(
    const int4 ns,                    // {npoints, ntris, _, _}
    __global const float4*  points,   // [npoints]{x,y,z,face_id}
    __global const float4*  faces,    // [npoints]{hx,hy,hz,area}
    __global       float4*  obstacles,// [ntris*3]{A,B,C vertices as float4, w=face_id}
    //__global       int*     faces_surf,     // surface index for each face
    //__global       float4*  obstacle_surt,  // surface index for each obstacle
    __global       float*  M              // [npoints*npoints] coupling matrix
){
    __local  float4 LOC[32*3]; // local copy of triangle points for obstacles
    //__local  int    LIs[32]; // local copy of indexes of surfaces for obstacles
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);

    float4 p1   = points[iG];
    float4 S1  = faces  [iG];
    //int    is1 = faces_surf[iG];
    int    is1 = (int)(p1.w); // we store surface index as float to save memory access    
    //for(int j=0; j<iG; j++){
    for(int i=0; i<ns.x; i++){
        float4 p2  = points[i];
        float4 S2  = faces [i];
        //int    is2 = faces_surf[j];
        int    is2 = (int)(p2.w); 

        // ---- calculate coupling coefficient depending on distance and angle between normals and connecting vector
        float3  d = p2.xyz - p1.xyz;
        float r = length(d);
        d /= r;
        float  coupling = diffusionLightCoupling( d, r*r, S1, S2 );

        // ---- calculate local coordinate system for the ray
        float3 hX = getSomeUp( d);
        float3 hY = cross( d, hX );

        // ---- loop over obstacles - accelerate by loading them to local memory
        //double occlusion = getOcclusion( eli.pos, d, sqrt(r2), eli.isurf, elj.isurf );
        float occlusion = 0.0f;
        for (int i0=0; i0<ns.y; i0+= nL ){ 
            int i3 = (i0+iL)*3;
            LOC[iL*3  ] = obstacles[i3  ];
            LOC[iL*3+1] = obstacles[i3+1];
            LOC[iL*3+2] = obstacles[i3+2];
            barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL] 
            for (int j=0; j<nL; j++){
                int   ip = (int)(LOC[j  ].w*1000.0f);
                float3 a = LOC[j  ].xyz;
                float3 b = LOC[j+1].xyz;
                float3 c = LOC[j+2].xyz;
                if( (ip==is1) || (ip==is2) ) continue; // skip self
                if( rayInTriangle( a,b,c, hX, hY ) ){
                    occlusion = 1.0f;
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
        }

        // ---- store coupling coefficient
        coupling*=(1-occlusion);
        coupling = fabs( coupling );
        M[iG*ns.x+i] = coupling;
        //M[i*ns.x+iG] = coupling;
        //nvalid++;
    }
}

// -----------------------------------------------------------------------------
// Kernel: occlusion_matrix
// Brute-force occlusion between element centers using triangle obstacles.
// One work-item per source point i; tests rays to all j against all triangles.
// Keeps as a simple reference implementation for debugging and validation.
// Inputs:
//  - points: [npoints] float4 (xyz + w=face_id)
//  - tris:   [ntris*3] float4 sequence of A,B,C (xyz) with w=face_id (from face)
// Output:
//  - occ:    [npoints*npoints] float, 1.0 if any triangle blocks i->j, else 0.0
__kernel void occlusion_matrix(
    const int npoints,                 // number of element centers
    const int ntris,                   // number of triangles
    __global const float4* points,     // [npoints]{x,y,z,face_id}
    __global const float4* tris,       // [ntris*3]{A,B,C vertices as float4, w=face_id}
    __global       float*  occ         // [npoints*npoints] 1 if blocked else 0
){
    const int i = get_global_id(0);
    if(i >= npoints) return;

    const float4 p1 = points[i];
    const int    s1 = (int)(p1.w + 0.5f);

    for(int j=0; j<npoints; ++j){
        if(j==i){ occ[i*npoints + j] = 0.0f; continue; }
        const float4 p2 = points[j];
        const int    s2 = (int)(p2.w + 0.5f);

        float3 d = p2.xyz - p1.xyz;
        float  r = length(d);
        if(r <= 1e-12f){ occ[i*npoints + j] = 0.0f; continue; }
        d *= 1.0f/r;
        float3 hX = getSomeUp(d);
        float3 hY = cross(d,hX);

        float  hit = 0.0f;
        for(int k=0; k<ntris; ++k){
            const int k3 = k*3;
            const float4 A4 = tris[k3  ];
            const float4 B4 = tris[k3+1];
            const float4 C4 = tris[k3+2];
            const int ts = (int)(A4.w + 0.5f);
            if((ts==s1) || (ts==s2)) continue; // skip triangles from either endpoint's face
            // Shift triangle by p1 so we test along ray from origin
            float3 A = A4.xyz - p1.xyz;
            float3 B = B4.xyz - p1.xyz;
            float3 C = C4.xyz - p1.xyz;
            // Quick param-range rejection: entirely before start or beyond end
            float tA = dot(d, A);
            float tB = dot(d, B);
            float tC = dot(d, C);
            if( (tA<=0.0f && tB<=0.0f && tC<=0.0f) ) continue;
            if( (tA>=r     && tB>=r     && tC>=r    ) ) continue;
            if( rayInTriangle(A,B,C,hX,hY) ){ hit = 1.0f; break; }
        }
        occ[i*npoints + j] = hit;
    }
}

float2 coneTriangleDistance( 
    float4 p1, // {x,y,z,R} 1st point of the cone
    float4 p2, // {x,y,z,R} 2nd point of the cone
    float3 a, //  {x,y,z} 1st point of the triangle
    float3 b, //  {x,y,z} 2nd point of the triangle
    float3 c  //  {x,y,z} 3rd point of the triangle
){
    float3 ax  = p2.xyz - p1.xyz;     // axis  of the cone
    float3 lax = length(ax);          // length of the axis squared
    ax *= 1/lax;
    // NOTE: we will use as origin p0 which is close to cog of the triangle for better numerical stability
    float3 cog  = (a+b+c)/3.0f;   // center of the triangle
    float3 p0   = ax*dot(ax,cog); // projection of cog to the axis
    // lets define coordinate system alligned with the cone
    float3 dab = b-a; dab -= ax*dot(dab,ax); // projection in plane perpendicular to the axis
    float3 dac = c-a; dac -= ax*dot(dac,ax); // projection in plane perpendicular to the axis
    float lab2 = dot(dab,dab); // squared length of the projection of edge ab to the plane
    float lac2 = dot(dac,dac); // squared length of the projection of edge ac to the plane
    // take the loger edge as the up-axis
    float3 up;
    if(lab2>lac2){
        up = dab/sqrt(lab2);
    }else{
        up = dac/sqrt(lac2);
    }

} 



// NOTE: Occlusion matrix can be used also to calculate visibility between larger groups of points with finite radius
#define nOccMax 8
// Kernel: makeOcclusionMatrix
// Broad-phase per-pair occluder listing using distance-to-ray bands.
// For each macro pair (i,j):
//  - If any obstacle closer than R_full -> full occlusion (bPartial=false)
//  - If all obstacles farther than R_safe -> no occlusion
//  - Else collect triangle indices of potential occluders up to nOccMax
__kernel void  makeOcclusionMatrix(
    const int4 ns,                    // {npoints, nTris, _, _}
    __global const float4*  points    ,// [npoints] {x,y,z,face_id}
    __global       float4*  obstacles ,// [nTris*3] {A,B,C as float4, w=face_id}
    __global       float*   distMin   ,// [npoints*npoints] min distance of any obstacle to ray (i->j)
    __global       int*     occluders ,// [npoints*npoints*nOccMax] triangle indices (or -1) per (i,j)
    __global       int*     noccs     ,// [npoints*npoints] number of collected occluders per (i,j) (clamped)
    float2 Rrange                      // {R_full,R_safe}
){
    __local  float4 LOC[32*3]; // local copy of triangle points for obstacles
    //__local  int    LIs[32]; // local copy of indexes of surfaces for obstacles
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    
    float4 p1  = points[iG];
    int    is1 = (int)(p1.w); 
    for(int i=0; i<ns.x; i++){
        float4 p2  = points[i];
        int    is2 = (int)(p2.w); 

        float3 d = p2.xyz - p1.xyz;
        float  r = length(d);
        if(r <= 1e-12f){
            // Degenerate pair (same point); nothing to do
            noccs  [iG*ns.x+i] = 0;
            distMin[iG*ns.x+i] = 1e30f;
            continue;
        }
        d /= r;
        
        // ---- calculate local coordinate system for the ray
        float3 hX = getSomeUp( d);
        float3 hY = cross( d, hX );

        // ---- loop over obstacles - accelerate by loading them to local memory
        //double occlusion = getOcclusion( eli.pos, d, sqrt(r2), eli.isurf, elj.isurf );
        int   nocc     = 0;
        float dist     = 1e30f;
        bool  bPartial = true;
        int   ijocc    = (iG*ns.x + i) * nOccMax;  // base offset into occluders for pair (iG,i)
        // Clear occluder list for this pair using all local threads
        for(int t=iL; t<nOccMax; t+=nL){ occluders[ ijocc + t ] = -1; }
        for (int i0=0; i0<ns.y; i0+= nL ){ 
            int i3 = (i0+iL)*3;
            LOC[iL*3  ] = obstacles[i3  ];
            LOC[iL*3+1] = obstacles[i3+1];
            LOC[iL*3+2] = obstacles[i3+2];
            barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL] 
            for (int j=0; j<nL; j++){
                int   ip = (int)(LOC[j].w);
                float3 a = LOC[j  ].xyz;
                float3 b = LOC[j+1].xyz;
                float3 c = LOC[j+2].xyz;
                if( (ip==is1) || (ip==is2) ) continue; // skip same surface 
                float rj = sdRayTrinagle( a, b, c, hX, hY );
                if( rj<Rrange.x ) bPartial = false; // Full occlusion
                if( bPartial ){
                    if( rj<Rrange.y ){
                        if(nocc < nOccMax){ occluders[ ijocc + nocc ] = i0 + j; }
                        nocc++;
                    }
                }
                dist = min(dist, rj);
            }
            barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
        }
        // Store number of collected occluders, clamped to capacity
        noccs  [iG*ns.x+i] = (nocc > nOccMax) ? nOccMax : nocc;
        distMin[iG*ns.x+i] = dist;
        //M[i*ns.x+iG] = coupling;
        //nvalid++;
    }
}

// =============================================================================
// Narrow-phase occlusion between sub-elements of a macro pair (i,j)
// Processes one macro pair per work-group. Sub-elements for macro i and j
// are loaded into local memory (bounded by MAX_SUB_I/J). Occluders are
// provided as group IDs from the broad-phase; we iterate their triangles
// in global memory in tiles (not fully cached here to keep memory bounded).
//
// Assumptions and tunables:
//  - MAX_SUB_I, MAX_SUB_J: max sub-elements per macro element loaded to local.
//  - WG_SIZE: work-group size; we stride over ni*nj pairs by WG_SIZE.
//  - nOccMax: number of occluders per pair from broad-phase (already defined).
//
// NOTE: This is a skeleton focused on structure and memory movement; details like
// exact fractional occlusion accumulation are left as TODO.

#ifndef MAX_SUB_I
#define MAX_SUB_I 32
#endif
#ifndef MAX_SUB_J
#define MAX_SUB_J 32
#endif
#ifndef TILE_TRIS
#define TILE_TRIS 64   // triangles per tile loaded to local (tune by local mem size)
#endif

__kernel void makeNarrowOcclusion(
    const int4   ns,            // {npoints, nOccMax, MAX_SUB_I, MAX_SUB_J}
    __global const int2*   pairsIJ,   // [nPairs]{i,j} macro indices; one per work-group
    __global const float4* subPoints, // [nSubTotal]{x,y,z,face_id} all sub-element centers, grouped by macro
    __global const int2*   subRanges, // [npoints]{offset,count} subPoints range per macro element
    __global const float4* tris,      // [nTris*3]{A,B,C as float4, w=face_id} concatenated triangle vertices
    __global const int2*   triRanges, // [nOccGroups]{offset3,count3} ranges into tris (in vertices, multiple of 3)
    __global const int*    occluders, // [npoints*npoints*nOccMax]{groupId|-1} occluder list per (i,j)
    __global const int*    noccs,     // [npoints*npoints]{count} valid occluders per (i,j), clamped ≤ nOccMax
    __global       float*  occPairs   // [nPairs]{fraction_blocked in 0..1} output per processed pair
){
    const int gid  = get_global_id(0);
    const int lid  = get_local_id (0);
    const int lsz  = get_local_size(0);

    // Each work-group handles exactly one pair
    if(lid==0 && get_num_groups(0)>0){ /* placeholder to silence unused warnings */ }
    const int pairId = get_group_id(0);
    if(pairId >= ns.x*ns.x) return; // conservative bound if caller doesn't compact pairs

    const int2 ij = pairsIJ[pairId];
    const int i = ij.x;
    const int j = ij.y;

    // Resolve sub-element ranges (clamp to MAX_SUB_*)
    const int2 ri = subRanges[i];
    const int2 rj = subRanges[j];
    const int ni = (ri.y > MAX_SUB_I) ? MAX_SUB_I : ri.y;
    const int nj = (rj.y > MAX_SUB_J) ? MAX_SUB_J : rj.y;

    // Local caches for sub-elements of i and j
    __local float4 LI[MAX_SUB_I];
    __local float4 LJ[MAX_SUB_J];
    for(int t=lid; t<ni; t+=lsz){ LI[t] = subPoints[ri.x + t]; }
    for(int t=lid; t<nj; t+=lsz){ LJ[t] = subPoints[rj.x + t]; }
    barrier(CLK_LOCAL_MEM_FENCE);

    // Read occluder list for this pair
    const int occBase = (i*ns.x + j) * ns.y;  // ns.y == nOccMax
    const int nOcc    = noccs[i*ns.x + j];

    // Accumulate blocked sub-pairs count
    int blocked = 0;
    const int nPairs = ni * nj;

    // Local tile for occluder triangles
    __local float4 TRILOC[3*TILE_TRIS];

    // Each thread iterates over its strided set of sub-pairs
    for(int kk = lid; kk < nPairs; kk += lsz){
        const int ik = kk / nj;
        const int jk = kk - ik*nj;
        const float4 P1 = LI[ik];
        const float4 P2 = LJ[jk];

        float3 d = P2.xyz - P1.xyz;
        float  r = length(d);
        if(r <= 1e-12f) continue; // same sub-point
        d *= 1.0f/r;
        float3 hX = getSomeUp(d);
        float3 hY = cross(d,hX);

        int hit = 0;
        // Iterate occluder groups listed by broad-phase
        for(int o=0; o<nOcc && !hit; ++o){
            const int gidOcc = occluders[occBase + o];
            if(gidOcc < 0) continue;
            const int2 tr = triRanges[gidOcc];  // {offset3, count3}

            // Tile over triangles of this occluder
            for(int base = tr.x; base < tr.x + tr.y && !hit; base += 3*TILE_TRIS){
                // cooperative load of up to TILE_TRIS triangles (3 verts each)
                int nverts = tr.x + tr.y - base;
                if(nverts > 3*TILE_TRIS) nverts = 3*TILE_TRIS;
                for(int t = lid; t < nverts; t += lsz){
                    TRILOC[t] = tris[base + t];
                }
                barrier(CLK_LOCAL_MEM_FENCE);

                // scan loaded triangles
                for(int v = 0; v < nverts; v += 3){
                    const float4 A4 = TRILOC[v  ];
                    const float4 B4 = TRILOC[v+1];
                    const float4 C4 = TRILOC[v+2];
                    // Optional: skip same-surface hits
                    const int ts = (int)(A4.w + 0.5f);
                    const int s1 = (int)(P1.w + 0.5f);
                    const int s2 = (int)(P2.w + 0.5f);
                    if((ts==s1) || (ts==s2)) continue;
                    float3 A = A4.xyz - P1.xyz;
                    float3 B = B4.xyz - P1.xyz;
                    float3 C = C4.xyz - P1.xyz;
                    // Quick t-range rejection
                    float tA = dot(d, A);
                    float tB = dot(d, B);
                    float tC = dot(d, C);
                    if( (tA<=0.0f && tB<=0.0f && tC<=0.0f) ) continue;
                    if( (tA>=r     && tB>=r     && tC>=r    ) ) continue;
                    if( rayInTriangle(A,B,C,hX,hY) ){ hit = 1; break; }
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }
        }
        blocked += hit;
    }

    // Reduce blocked across work-group (simple atomic add to global tmp or local reduction)
    // For simplicity, use local reduction to one thread then write fraction
    __local int Lsum;
    if(lid==0) Lsum = 0;
    barrier(CLK_LOCAL_MEM_FENCE);
#if __OPENCL_VERSION__ >= 120
    atomic_add(&Lsum, blocked);
#else
    // Fallback: naive reduce with one-by-one accumulation
    // Each thread writes to a local slot and thread 0 sums
    __local int Lbuf[256]; // assumes WG_SIZE<=256; tune as needed
    if(lid<256) Lbuf[lid] = blocked; // threads beyond 256 ignored
    barrier(CLK_LOCAL_MEM_FENCE);
    if(lid==0){ int s=0; for(int t=0;t<lsz && t<256;t++) s+=Lbuf[t]; Lsum=s; }
#endif
    barrier(CLK_LOCAL_MEM_FENCE);
    if(lid==0){
        const float frac = (nPairs>0) ? ((float)Lsum)/((float)nPairs) : 0.0f;
        occPairs[pairId] = frac;
    }
}