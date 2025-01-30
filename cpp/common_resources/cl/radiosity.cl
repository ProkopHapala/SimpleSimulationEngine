


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



__kernel void  makeRadiosityCouplings(
    const int4 ns, 
    __global const float4*  points,    // x,y,z,mass
    __global const float4*  faces,     // hx,hy,hz, area
    __global       float4*  obstacles, // (x,y,z)*3 - triangles
    //__global       int*     faces_surf,     // surface index for each face
    //__global       float4*  obstacle_surt,  // surface index for each obstacle
    __global       float*  M // coupling matrix
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



// NOTE: Occlusion matrix can be used also to calculate visiibility between larger groups of points with finite radius

__kernel void  makeOcclusionMatrix(
    const int4 ns, 
    __global const float4*  points    ,    // [npoints] {x,y,z,R}    centers of point-groups with radius R
    __global       float4*  obstacles ,    // [tris   ] {ia,ib,ic,?} triangles representing obstacle-groups
    __global       float*   distMin   ,    // [npoints^2        ] minimum distance of any obstacle to ray between point-groups i,j
    __global       int*     occluders ,    // [npoints^2*nOccMax] list of obstacles which may occlude i,j
    float2 Rrange                          // {R_full,R_safe} minimum and maximum distance for partial occlusion
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
        d /= r;
        
        // ---- calculate local coordinate system for the ray
        float3 hX = getSomeUp( d);
        float3 hY = cross( d, hX );

        // ---- loop over obstacles - accelerate by loading them to local memory
        //double occlusion = getOcclusion( eli.pos, d, sqrt(r2), eli.isurf, elj.isurf );
        int   nocc     = 0;
        float dist     = 0.0f;
        bool  bPartial = true;
        int   ijocc    = (iG*ns.x+i)*nocc;
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
                        occluders[ ijocc + nocc ]=i0+j;
                        nocc++;
                    }
                }
                dist = min(dist,rj);
            }
            barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
        }
        distMin[iG*ns.x+i] = dist;
        //M[i*ns.x+iG] = coupling;
        //nvalid++;
    }
}