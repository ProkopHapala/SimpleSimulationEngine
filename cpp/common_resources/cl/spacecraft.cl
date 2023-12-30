
// ====================== Truss Simulation ( Soft Body ) ======================



// ========== Try Device Side Enqueue ==========
// see: https://github.com/ProkopHapala/FireCore/wiki/OpenCL-performance-tips 
//   *  https://bloerg.net/posts/opencl-2-0-is-here/
//   *  https://www.codeproject.com/Articles/867780/GPU-Quicksort-in-OpenCL-Nested-Parallelism-and-Wor#05_converting


// __kernel void  test_enque(
//     //const int4 ns, 
//     //__global const float4*  inp,  
//     //__global       float4*  work,
//     //int n,  
//     __global       float4*  out
// ){
//     const int i = get_global_id(0);
//     out[i] = (float4){ i*1.0f , 0.0f , 0.0f , 1.0f };
    
//     //if(i==0){
//         //ndrange_t range = ndrange_1D( get_global_size(0), get_local_size(0) );
//         ndrange_t range = ndrange_1D( 5, 1 );
//         enqueue_kernel( get_default_queue(), CLK_ENQUEUE_FLAGS_WAIT_KERNEL, range,
//         ^{
//             int i  = get_global_id (0);
//             if(i==0){ printf( "GPU::inside enqueue_kernel block \n" ); }
//             out[i] = (float4){ i*2.0f , 0.0f , 0.0f , 2.0f };
//         });
//     //}
    
// }

// Problem - we want to avoid asynchronious memory writes, therefore we need to iterate over vertexes rather than over edges

float4 springForce( float3 d, float4 par ){
    float l  = length(d);
    float dl = l - par.x;
    float f;
    if( dl > 0.0f ){ 
        f = par.z * -dl;
    } else {
        f = par.y * dl;
    }
    return (float4){d*f, f*dl*0.5f};
}

__kernel void  evalTrussForce1(
    const int4 ns, 
    __global const float4*  points, // x,y,z,mass
    __global const float4*  vels,   // velocities are used for damping 
    __global       float4*  forces, 
    __global const int*     neighs, // indexes of neighbor points, if neighs[i] == -1 it is not connected
    __global const float4*  params, // l0, kPress, kPull, damping
    float4 accel, // acceleration of the reference frame
    float4 omega, // angular velocity for simulation of rotating reference frame
    float4 rot0   // center of rotation
){
    const int iG = get_global_id(0);

    //if(iG==0){ printf("GPU::evalTrussForce2() \n" ); }
    if(iG>=ns.x) return;
    float4 p = points[iG];
    float4 v = vels  [iG];
    float4 f =(float4){0.0f,0.0f,0.0f,0.0f};
    
    for(int ij=0; ij<ns.y; ij++){
        int  j  = ns.y*iG + ij;
        int  ja = neighs[j];
        if(ja==-1) break;
        //f += springForce( points[b.x].xyz - p.xyz, bparams[b.y] );
        float3 d = points[ja].xyz - p.xyz;
        float l  = length(d.xyz);
        float dl = l - params[j].x;
        float k = 1e+6f;
        // f.f.add_mul( d, (k*(li-params[j].x)/li) );
        f.xyz += d.xyz * (k*dl/l);
        f.w   += dl*dl*0.5f;
    }
    
    //if(iG==0){ printf("accel(%g,%g,%g,%g) omega(%g,%g,%g,%g) rot0(%g,%g,%g,%g) \n" , accel.x,accel.y,accel.z,accel.w, omega.x,omega.y,omega.z,omega.w, rot0.x,rot0.y,rot0.z,rot0.w ); }
    // acceleration of the reference frame
    f.xyz += accel.xyz*p.w;
    // centrifugal force
    float3 d = p.xyz - rot0.xyz;
    f.xyz += (d - omega.xyz * dot(omega.xyz,d)) * omega.w*omega.w*p.w; 
    //coriolis force (depends on velocity) => we cannot calculate it here
    
    forces[iG] = f; // we may need to do += in future 
}

__kernel void  evalTrussForce2(
    const int4 ns, 
    __global const float4*  points,    // x,y,z,mass
    __global const float4*  vels,      // velocities are used for damping 
    __global       float4*  forces, 
    //__global const int*   neighs,  // indexes of neighbor points, if neighs[i] == -1 it is not connected
    __global const int2*    neighBs,   // indexes of neighbor (point,bond), if neighs[i].x == -1 it is not connected
    __global const float4*  bparams,   // l0, kPress, kPull, damping
    float4 accel, // acceleration of the reference frame
    float4 omega, // angular velocity for simulation of rotating reference frame
    float4 rot0   // center of rotation
){
    const int iG = get_global_id(0);

    //if(iG==0){ printf("GPU::evalTrussForce2() \n" ); }

    if(iG>=ns.x) return;
    float4 p = points[iG];
    float4 v = vels  [iG];
    float4 f =(float4){0.0f,0.0f,0.0f,0.0f};
    
    for(int ij=0; ij<ns.y; ij++){
        int  j = ns.y*iG + ij;
        int2 b = neighBs[j];
        if(b.x==-1) break;
        //f += springForce( points[b.x].xyz - p.xyz, bparams[b.y] );

        float3 d = points[b.x].xyz - p.xyz;
        float l  = length(d.xyz);
        float dl = l - bparams[b.y].x;
        float k = 1e+6f;
        // f.f.add_mul( d, (k*(li-params[j].x)/li) );
        f.xyz += d.xyz * (k*dl/l);
        f.w   += dl*dl*0.5f;
    }
    //if(iG==0){ printf("accel(%g,%g,%g,%g) omega(%g,%g,%g,%g) rot0(%g,%g,%g,%g) \n" , accel.x,accel.y,accel.z,accel.w, omega.x,omega.y,omega.z,omega.w, rot0.x,rot0.y,rot0.z,rot0.w ); }
    //if(iG==1074){ printf("p[%i](%g,%g,%g|%g) v(%g,%g,%g|%g) \n", iG, p.x,p.y,p.z,p.w,  v.x,v.y,v.z,v.w ); }
    // acceleration of the reference frame
    f.xyz += accel.xyz*p.w;
    // centrifugal force
    float3 d = p.xyz - rot0.xyz;
    f.xyz += (d - omega.xyz * dot(omega.xyz,d)) * omega.w*omega.w*p.w; 
    //coriolis force (depends on velocity) => we cannot calculate it here
    forces[iG] = f; // we may need to do += in future 
}

__kernel void  evalTrussBondForce(
    const int4 ns, 
    __global const float4*  points,    // x,y,z,mass
    __global       float4*  bforces,   // bond forces
    __global const int2*    bonds,     // indexes of neighbor (point,bond), if neighs[i].x == -1 it is not connected
    __global const float4*  bparams    // l0, kPress, kPull, damping
){
    const int iG = get_global_id(0);
    //if(iG==0){ printf("GPU::evalTrussBondForce(nBonds=%i,nNeigh=%i)\n", ns.x, ns.y ); }
    int2   b   = bonds  [iG];
    float4 par = bparams[iG];
    float3 d   = points [b.y].xyz - points[b.x].xyz;
    float l    = length(d);
    float dl   = l - par.x;
    float k    = 1e+6f;
    float4 f = (float4){
        d*(k*dl/l), // force 
        dl*dl*0.5f  // potential energy
    };
    bforces[iG] = f; // we may need to do += in future 
    //bforces[iG] = (float4)( iG ,0.0,1.0,2.0);
}

__kernel void  assembleAndMove(
    const int4 ns, 
    float4 MDpars,
    __global       float4* points,    // x,y,z,mass
    __global       float4* velocities, 
    __global       float4* forces,
    //__global const float4* forces,
    __global const int*    neighB2s,  // index of bond for each neighbor, if neighs[i]==0 it is not connected, if negative it is opposite direction
    __global const float4* bforces,
    float4 accel, // acceleration of the reference frame
    float4 omega, // angular velocity for simulation of rotating reference frame
    float4 rot0   // center of rotation
){
    
    const int iG = get_global_id(0);
    if(iG>=ns.x) return;
    //if(iG==0){ printf("GPU::assembleAndMove(nPoint=%i,nNeigh=%i)\n", ns.x, ns.y ); }

    // ------ Assemble bond forces
    //float4 f = forces    [iG];
    float4 f = (float4){0.0f,0.0f,0.0f,0.0f};
    for(int ij=0; ij<ns.y; ij++){
        int  j  = ns.y*iG + ij;
        int  ib = neighB2s[j];
        //if(iG==0){ float4 fb = bforces[ ib-1]; printf("GPU::bond[iG=%i,j=%i,ib=%i]  \n", iG, ij, ib ); }
        if(ib==0) break;
        if(ib>0){
            f += bforces[ ib-1];
        }else{
            f -= bforces[-ib-1];
        }
        //if(iG==0){ float4 fb = bforces[ ib-1]; printf("GPU::bond[iG=%i,j=%i,b=%i] fb(%g,%g,%g|%g) \n", iG, ij, ib, fb.x,fb.y,fb.z,fb.w ); }
    }
    forces[iG] = f;   // DEBUG

    // ------ apply local (pointwise) forces
    float4       p = points    [iG];
    float4       v = velocities[iG];
    f.xyz += accel.xyz*p.w;
    // centrifugal force
    float3 d = p.xyz - rot0.xyz;
    f.xyz += (d - omega.xyz * dot(omega.xyz,d)) * omega.w*omega.w*p.w; 
    // //coriolis force (depends on velocity) => we cannot calculate it here

    // ------ Move (Leap-Frog)
    //if(iG==0){ printf("GPU::move() MDpars(%g,%g,%g,%g)\n" , MDpars.x,MDpars.y,MDpars.z,MDpars.w ); }
    
    v.xyz *= MDpars.y;
    //v.xyz *= 0.95f;
    v.xyz += f.xyz*MDpars.x/p.w;
    p.xyz += v.xyz*MDpars.x;
    // ToDo: something like FIRE ?
    velocities[iG] = v;
    points    [iG] = p;
}

// ====================== Magnetic Interactions ( Amber / Boist-Sawart / Lorenz ) ======================

// https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law
// https://en.wikipedia.org/wiki/Amp%C3%A8re%27s_force_law#General_case
// https://en.wikipedia.org/wiki/Lorentz_force#Force_on_a_current-carrying_wire

//float3 forceLorentz( float3 Bj, float3 Ij ){ return  cross( Ij, B ); }
float3 fieldBiotSavart( float3 d, float3 Ij ){ float r2 = dot(d,d); return cross( Ij, d ) / ( r2* sqrt(r2) ); }
float3 fieldAmper     ( float3 d, float3 Ii, float3 Ij ){ float r2 = dot(d,d); return cross( Ii, cross( Ij, d ) ) / ( r2* sqrt(r2) ); }

__kernel void magnetism_On2(
    const int4 ns,
    __global const float4*  points,    // x,y,z,mass
    __global const float4*  currents,  // Ix,Iy,Iz, Q
    __global       float4*  forces
){
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    __local  float8 LOC[32];
    float4 p   = points  [iG];
    float3 B =(float3){0.0f,0.0f,0.0f};
    for (int i0=0; i0<ns.x; i0+= nL ){ 
        LOC[iL] = (float8){ points[i0 + iL], currents[i0 + iL] }; 
        barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL] 
        for (int j=0; j<nL; j++){
            B += fieldBiotSavart( LOC[j].xyz-p.xyz, LOC[j].hi.xyz );
        }
        barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
    }
    forces[iG] = (float4){ cross( currents[iG].xyz, B ), 0.0f }; // Lorentz force
}

// ====================== Movement ======================

__kernel void  move(
    const int4 ns, 
    float4 MDpars,
    __global       float4* points,    // x,y,z,mass
    __global       float4* velocities, 
    __global const float4* forces
){
    const int iG = get_global_id(0);
    if(iG>=ns.x) return;
    // ------ Move (Leap-Frog)
    float4       p = points    [iG];
    float4       v = velocities[iG];
    const float4 f = forces    [iG];
    //if(iG==0){ printf("GPU::move() MDpars(%g,%g,%g,%g)\n" , MDpars.x,MDpars.y,MDpars.z,MDpars.w ); }
    v.xyz *= MDpars.y;
    v.xyz += f.xyz*MDpars.x/p.w;
    p.xyz += v.xyz*MDpars.x;
    // ToDo: something like FIRE ?
    velocities[iG] = v;
    points    [iG] = p;
}


// ====================== Radiosity and ray-tracing ======================

float3 getSomeUp( float3 v ){
	if( v.x<v.y){ return (float3){ -v.y*v.y -v.z*v.z, v.x*v.y           , v.x*v.z  }; }
    else        { return (float3){  v.y*v.x         , -v.z*v.z -v.x*v.x ,  v.y*v.z }; }
}

bool rayInTriangle( float3 a_, float3 b_, float3 c_, float3 hX, float3 hY ){
	float2 a = (float2){ dot(hX,a_), dot(hY,a_) };
    float2 b = (float2){ dot(hX,b_), dot(hY,b_) };
    float2 c = (float2){ dot(hX,c_), dot(hY,c_) };
    float   sgn = a.x*(b.y-a.y) - a.y*(b.x-a.x);
	if( 0 > sgn*( b.x*(c.y-b.y) - b.y*(c.x-b.x) ) ) return false;
	if( 0 > sgn*( c.x*(a.y-c.y) - c.y*(a.x-c.x) ) ) return false;
    //printf("passed\n");
    return true;
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
float getOcclusion( float3 ray0, float3 hRay, float tmax, int ip1, int ip2 ){
    if( triangleObstacles.size()==0 ) return 0.0;
    // Check occlusion - TODO can be made better
    Vec3d hX,hY;
    float3 hX = getSomeUp( hRay);
    float3 hY = cross( hRay, hX );
    for( int i=0; i<triangleObstacles.size(); i++ ){
        if( (i==ip1) || (i==ip2) ) continue; // skip self
        float3 a = obstacles[i  ];
        float3 b = obstacles[i+1];
        float3 c = obstacles[i+2];
        if( rayInTriangle( a,b,c, hX, hY ) ) return 1.0;
    }
    return 0.0;
}
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
    int    is1 = (int)(p1.w * 1000.0f); // we store surface index as float to save memory access    
    //for(int j=0; j<iG; j++){
    for(int i=0; i<ns.x; i++){
        float4 p2  = points[i];
        float4 S2  = faces [i];
        //int    is2 = faces_surf[j];
        int    is2 = (int)(p2.w * 1000.0f); 

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
        float occlusion = 0.0;
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
                    occlusion = 1.0;
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