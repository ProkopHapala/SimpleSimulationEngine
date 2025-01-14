
// ====================== Truss Simulation ( Soft Body ) ======================



// ========== Try Device Side Enqueue ==========
// see: https://github.com/ProkopHapala/FireCore/wiki/OpenCL-performance-tips 
//   *  https://bloerg.net/posts/opencl-2-0-is-here/
//   *  https://www.codeproject.com/Articles/867780/GPU-Quicksort-in-OpenCL-Nested-Parallelism-and-Wor#05_converting


__kernel void  test_enque(
    //const int4 ns, 
    //__global const float4*  inp,  
    //__global       float4*  work,
    //int n,  
    __global       float4*  out,
    __global       float4*  out2
){
    //const int i = get_global_id(0);
    //out[i] = (float4){ i*1.0f , 0.0f , 1.0f , 1.0f };
    if( 0 == get_global_id(0) ){  // the global_size is 1 anyway, but just to be sure


    ndrange_t range = ndrange_1D( 100, 1 );
    
    // this kernell just clears the output array
    //enqueue_kernel( get_default_queue(), CLK_ENQUEUE_FLAGS_WAIT_KERNEL, range,
    enqueue_kernel( get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, range,
    ^{
        const int iG  = get_global_id (0);
        //if(iG==0){ printf( "GPU::inside enqueue_kernel.clear \n" ); }
        out[iG] = (float4){ 0.0f, 0.0f , 0.0f , 2.0f };
        if(iG==10) out[iG].x = 1.0f;
        //barrier(CLK_GLOBAL_MEM_FENCE);

    });
    //barrier(CLK_GLOBAL_MEM_FENCE);
    for(int itr=0; itr<100; itr++){   // we submit 5 iterations of the same kernel      
        //enqueue_kernel( get_default_queue(), CLK_ENQUEUE_FLAGS_WAIT_KERNEL, range,
        enqueue_kernel( get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, range,
        ^{
            const int iG = get_global_id (0);
            const int nG = get_global_size(0);
            //if(iG==0){ printf( "GPU::inside enqueue_kernel.flip(%i) \n", itr ); }
            //out[iG] = out[iG] + (float4){ 1.0f, iG*1.f, itr , 0.0f };

            //out2[iG].x = out[iG].x + 1.0f;

            // blur
            if ((iG>0)||(iG<nG)) { out2[iG].x = (out[iG-1].x + out[iG].x + out[iG+1].x)/3.0f; }

            //pascal triangle
            // if ((iG==0) || (iG==nG-1)) {
            //    out2[iG].x = 1.0f;
            // }else{ 
            //    out2[iG].x = out[iG].x + out[iG-1].x;
            // }
            //barrier(CLK_GLOBAL_MEM_FENCE);
        });

        //enqueue_kernel( get_default_queue(), CLK_ENQUEUE_FLAGS_WAIT_KERNEL, range,
        enqueue_kernel( get_default_queue(), CLK_ENQUEUE_FLAGS_NO_WAIT, range,
        ^{
            const int iG = get_global_id (0);
            const int nG = get_global_size(0);
            //if(iG==0){ printf( "GPU::inside enqueue_kernel.flop(%i) \n", itr ); }
            //out[iG] = out[iG] + (float4){ 1.0f, iG*1.f, itr , 0.0f };

            //out2[iG].x = out[iG].x + out[iG-1].x;

            //out[iG].x = out2[iG].x + 1.0f;
            // blur
            if ((iG>0)||(iG<nG)) { out[iG].x = (out2[iG-1].x + out2[iG].x + out2[iG+1].x)/3.0f; }

            //pascal triangle
            // if ((iG==0) || (iG==nG)) {
            //    out[iG].x = 1.0f;
            // }else{ 
            //    out[iG].x = out2[iG].x + out2[iG-1].x;
            // }
            //barrier(CLK_GLOBAL_MEM_FENCE);
        });

    }
   
    }
}

__kernel void  test_blur(
    __global       float4*  inp,
    __global       float4*  out
){
    const int iG = get_global_id(0);
    const int nG = get_global_size(0);
    if ((iG>0)||(iG<nG)) { out[iG].x = (inp[iG-1].x + inp[iG].x + inp[iG+1].x)/3.0f; }
}    


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



__kernel void getTrussForces( 
    int npoint,                     // 1 number of points
    int nmax_neigh,                 // 2 max number of neighbors
    __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    __global       float4*  forces, // 4 [npoint] x,y,z,E
    __global const int*     neighs, // 5 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global const float4*  params, // 6 [npoint,nmax_neigh] {l0, kPress, kPull, damping} 
    float inv_dt2                   // 7 1/dt^2, controls scale of inertial term
){
    const int iG     =  get_global_id(0);
    if(iG>=npoint) return;
    const int j0     = iG * nmax_neigh;
    const float4 pi  = ps[iG];                       
    float4       fi  = (float4){0.0f,0.0f,0.0f,0.0f}; 
    for( int jj = 0; jj < nmax_neigh; jj++ ){
        const int j  = j0 + jj;
        const int jG = neighs[j];
        if( jG == -1 ) break;
        const float4 par = params[j];
        const float3 pj  = ps[jG].xyz;
        const float3 dij = pi.xyz - pj;
        const float  l   = length(dij);
        const float  dl  = l - par.x; 
        const float  k   = (dl>0.f)?par.y:par.z;
        const float  fr  = dl*k;
        fi.xyz          += dij  * fr/l  ;
        fi.w            += 0.5f * fr*dl ;
    }
    forces[iG] = fi;
}


__kernel void updatePD_RHS( 
    int npoint,                     // 1 number of points
    int nmax_neigh,                 // 2 max number of neighbors
    __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    __global       float4*  bvec,   // 4 [npoint] x,y,z,Aii
    __global const int*     neighs, // 5 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global const float4*  params, // 6 [npoint,nmax_neigh] {l0, kPress, kPull, damping} 
    float inv_dt2                   // 7 1/dt^2, controls scale of inertial term
){

    //printf( "GPU::updateJacobi_mix() \n" );
    const int iG     =  get_global_id(0);
    if(iG>=npoint) return;
    const float4 pi  = ps[iG];                       
     
    // inertial term
    const float Aii0 = pi.w*inv_dt2;                   // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
    float4      bi   = (float4){pi.xyz * Aii0, Aii0};  // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )
    
    const int    j0  = iG * nmax_neigh;
    for( int jj = 0; jj < nmax_neigh; jj++ ){
        const int j   = j0 + jj;
        const int jG  = neighs[j];
        if( jG == -1 ) break;
        const float4 par = params[j];
        const float3 pj  = ps[jG].xyz;
        const float3 dij = pi.xyz - pj;
        const float  l   = length(dij);
        //const float k = (l<par.x)?par.y:par.z;
        const float k    = par.z;
        bi.xyz          += dij * (k*par.x/l) ;   // RHS      bi  = \sum_j K_{ij} d_{ij} (+ inertial term M_i/dt^2 p'_i )
        bi.w            +=        k          ;   // diagonal Aii =  \sum_j k_{ij}       (+ inertial term M_i/dt^2      )
    }
    bvec[iG] = bi;
}



__kernel void updateJacobi_lin( 
    int npoint,                     // 1 number of points
    int nmax_neigh,                 // 2 max number of neighbors
    __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    __global       float4*  ps_out, // 4 [npoint] x,y,z,mass
    __global       float4*  bvec,   // 5 [npoint] x,y,z,Aii      force
    __global const int*     neighs, // 6 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global const float*   kngs    // 7 [npoint,nmax_neigh] stiffness 
){
    //printf( "GPU::updateJacobi_lin() \n" );
    const int iG     =  get_global_id(0);
    if(iG>=npoint) return;

    //if(iG==0){    for(int i=0; i<npoint; i++){ printf("GPU::updateJacobi_lin(%3i) ps(%10.3e,%10.3e,%10.3e|%10.3e) bi(%10.3e,%10.3e,%10.3e|%10.3e)\n" , i, ps[i].x,ps[i].y,ps[i].z,ps[i].w, bvec[i].x,bvec[i].y,bvec[i].z,bvec[i].w );  } }

    const int j0     = iG * nmax_neigh;
    float3 sum_j = (float3){0.0f,0.0f,0.0f};                         
    for( int jj = 0; jj < nmax_neigh; jj++ ){
        const int j = j0 + jj;
        const int jG = neighs[j];
        if( jG == -1 ) break;
        sum_j       += ps[jG].xyz * kngs[j];   // kPull 
        //if(iG==0){ printf("GPU::updateJacobi_lin(%3i) j: %3i k: %10.3e pj(%10.3e,%10.3e,%10.3e)\n" , iG, kngs[j], ps[jG].x,ps[jG].y,ps[jG].z ); }
    }
    const float4 bi     =  bvec[iG];

    //if(iG==0){  printf("GPU::updateJacobi_lin(%3i) Aii(%10.3e) bi(%10.3e,%10.3e,%10.3e) sum_j(%10.3e,%10.3e,%10.3e)\n" , iG, bi.w, bi.x,bi.y,bi.z, sum_j.x,sum_j.y,sum_j.z );  }

    const float3 pi_new = (bi.xyz + sum_j)/bi.w;
    ps_out[iG] = (float4)(pi_new, ps[iG].w );
}



__kernel void updateJacobi_mix( 
    int npoint,                     // 1 number of points
    int nmax_neigh,                 // 2 max number of neighbors
    __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    __global       float4*  ps_out, // 4 [npoint] x,y,z,mass
    __global       float4*  dps,    // 4 [npoint] x,y,z,mass momentum
    __global const int*     neighs, // 5 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global const float4*  params, // 6 [npoint,nmax_neigh] {l0, kPress, kPull, damping} 
    float inv_dt2,                  // 7 1/dt^2, controls scale of inertial term
    float2 bmix                     // 8 mixing parameter for momentum
){

    //printf( "GPU::updateJacobi_mix() \n" );
    const int iG     =  get_global_id(0);
    if(iG>=npoint) return;

    const int j0     = iG * nmax_neigh;
    const float4 pi  = ps[iG];                       
     
    // inertial term
    const float mi = pi.w*inv_dt2;
    float3 pi_new  = pi.xyz * mi; // p_i    = \sum_j (K_{ij} p_j + d_{ij} )  + M_i/dt^2 p'_i
    float  Aii     =          mi; // A_{ii} = \sum_j  K_{ij}                 + M_i/dt^2

    for( int jj = 0; jj < nmax_neigh; jj++ ){
        const int j  = j0 + jj;
        const int jG = neighs[j];
        if( jG == -1 ) break;
        const float4 par = params[j];
        const float3 pj  = ps[jG].xyz;
        const float3 dij = pi.xyz - pj;
        const float  l   = length(dij);
        const float3 pij = pj + dij * (par.x/l);  // p'_{ij} is ideal position of particle i which satisfy the constraint between i and j
        const float k = (l<par.x)?par.y:par.z;
        pi_new += pij * k;   // 
        Aii    +=       k;   // diagonal of the Projective-Dynamics matrix, \sum_j k_{ij}       (+ inertial term M_i/dt^2)
    }

    pi_new /= Aii;

    // momentum mixing
    pi_new      = pi_new * bmix.x   +   dps[iG].xyz * bmix.y;
    dps[iG].xyz = pi_new - pi.xyz;

    ps_out[iG] = (float4)(pi_new, pi.w);

}

__kernel void updateJacobi_neighs( 
    int npoint,                     // 1 number of points
    int nmax_neigh,                 // 2 max number of neighbors
    __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    __global       float4*  ps_out, // 4 [npoint] x,y,z,mass
    __global const int*     neighs, // 5 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    __global const float4*  params, // 6 [npoint,nmax_neigh] {l0, kPress, kPull, damping} 
    float inv_dt2                   // 7 1/dt^2, controls scale of inertial term
){

    const int iG     =  get_global_id(0);
    if(iG>=npoint) return;

    // if(iG==0){
    //     printf( "updateJacobi_neighs() npoint(%i) nmax_neigh(%i) inv_dt2(%g) \n",  npoint, nmax_neigh, inv_dt2 );
    //     for(int i=0; i<npoint; i++){
    //         printf("GPU::updateJacobi_neighs(%3i) ps(%10.3e,%10.3e,%10.3e,%10.3e) " , i, ps[i].x,ps[i].y,ps[i].z,ps[i].w );
    //         for(int j=0; j<nmax_neigh; j++){
    //             int jG = neighs[i*nmax_neigh+j];
    //             float4 par = params[i*nmax_neigh+j];
    //             if( jG == -1 ) break;
    //             printf(" ngs[%i](%3i,%10.3f,%10.3f)" ,  j, jG, par.x, par.y );
    //         }
    //         printf("\n");
    //     }
    // }

    const int j0     = iG * nmax_neigh;
    const float4 pi  = ps[iG];                       
     
    // inertial term
    const float mi = pi.w*inv_dt2;
    float3 pi_new  = pi.xyz * mi; // p_i    = \sum_j (K_{ij} p_j + d_{ij} )  + M_i/dt^2 p'_i
    float  Aii     =          mi; // A_{ii} = \sum_j  K_{ij}                 + M_i/dt^2

    //float3 bi     = pi.xyz * mi;             // DEBUG 
    //float3 sum_j  = (float3){0.0f,0.0f,0.0f};  // DEBUG

    for( int jj = 0; jj < nmax_neigh; jj++ ){
        const int j  = j0 + jj;
        const int jG = neighs[j];
        if( jG == -1 ) break;
        const float4 par = params[j];
        const float3 pj  = ps[jG].xyz;
        const float3 dij = pi.xyz - pj;
        const float  l   = length(dij);
        const float3 pij = pj + dij * (par.x/l);  // p'_{ij} is ideal position of particle i which satisfy the constraint between i and j
        const float k = (l<par.x)?par.y:par.z;
        //printf("GPU::updateJacobi_neighs(iG=%3i,j=%i) jG(%3i) l(%10.3f) l0(%10.3f) k(%10.3f) \n" , iG, j, jG, l, par.x, k );
        pi_new += pij * k;   // 
        Aii    +=       k;   // diagonal of the Projective-Dynamics matrix, \sum_j k_{ij}       (+ inertial term M_i/dt^2)

        //bi     += dij * (k*par.x/l);   // DEBUG
        //sum_j  += pj  *  k;            // DEBUG
        //if(iG==0){ printf("GPU::updateJacobi_neighs(%3i) j: %3i k: %10.3e pj(%10.3e,%10.3e,%10.3e)\n" , iG, j, k, pj.x,pj.y,pj.z ); }

    }

    // NOTE 1: pi_new is basically weighted average of predicted positions of particles and due to all the constraints weighted by stiffness (and due to intertial term M_i/dt^2)
    // pi_new   = \sum_j pij * k_{ij} / \sum_j k_{ij}
    // NOTE 2: at the same time, pi_new can be seen as one step of Jacobi iteration for Projective-Dynamics equation Ap=b
    // pi_new   = (bi + sum_j) / A_{ii} 
    // where  bi     = \sum_j K_{ij} d_{ij}    (+ inertial term M_i/dt^2 pi )
    // and    A_{ii} = \sum_j K_{ij}           (+ inertial term M_i/dt^2    )

    //printf("GPU::updateJacobi_neighs(iG=%3i) Aii(%10.3f) pi_new(%10.3f,%10.3f,%10.3f) \n" , iG, Aii, pi_new.x,pi_new.y,pi_new.z );

    pi_new /= Aii;

   // if(iG==0){  printf("GPU::updateJacobi_neighs(%3i) Aii(%10.3e) bi(%10.3e,%10.3e,%10.3e) sum_j(%10.3e,%10.3e,%10.3e)\n" , iG, Aii, bi.x,bi.y,bi.z, sum_j.x,sum_j.y,sum_j.z );  }
    //pi_new = (bi + sum_j)/Aii;    // DEBUG



    ps_out[iG] = (float4)(pi_new, pi.w);

}

__kernel void PD_perdictor( 
    int npoint,                     // 1 number of points
    __global const float4*  ps,     // 2 [npoint] x,y,z,mass
    __global       float4*  ps_out, // 3 [npoint] x,y,z,mass
    __global const float4*  fs    , // 4 [npoint] x,y,z,E      force 
    __global const float4*  vs    , // 5 [npoint] x,y,z,?      velocity 
    //__global const float4*  dps , // 6 [npoint] x,y,z,?      solver momentum
    float dt                  
){
    const int iG     =  get_global_id(0);
    if(iG>=npoint) return;
    //if(iG==0){    for(int i=0; i<npoint; i++){ printf("GPU::PD_perdictor(%3i) ps(%10.3e,%10.3e,%10.3e,%10.3e) fs(%10.3e,%10.3e,%10.3e,%10.3e) vs(%10.3e,%10.3e,%10.3e,%10.3e) \n" , i, ps[i].x,ps[i].y,ps[i].z,ps[i].w, fs[i].x,fs[i].y,fs[i].z,fs[i].w, vs[i].x,vs[i].y,vs[i].z,vs[i].w );  } }
    const float4 pi  = ps[iG];                       
    float3 v   = vs[iG].xyz + fs[iG].xyz * (dt/pi.w); // Leap-Frog: v_{k+1/2} = v_{k-1/2} + f_k / m   dt
    ps_out[iG] = (float4){ pi.xyz + v*dt, pi.w };     //            p_{k+1  } = p_k       + v_{k+1/2} dt 
}

__kernel void PD_corrector( 
    int npoint,                     // 1 number of points
    __global const float4*  ps_new, // 2 [npoint] x,y,z,mass
    __global       float4*  ps_old, // 3 [npoint] x,y,z,mass
    __global       float4*  vs    , // 4 [npoint] x,y,z,?      velocity
    __global       float4*  dvs   , // 5 [npoint] x,y,z,?      Impulse
    __global       float4*  fs    , // 6 [npoint] x,y,z,?      force
    float dt                  
){
    const int iG     =  get_global_id(0);
    if(iG>=npoint) return;
    //if(iG==0){  for(int i=0; i<npoint; i++){ printf("GPU::PD_corrector(%3i) ps_new(%10.3e,%10.3e,%10.3e|%10.3e) ps_old(%10.3e,%10.3e,%10.3e) vs(%10.3e,%10.3e,%10.3e) dvs(%10.3e,%10.3e,%10.3e) fs(%10.3e,%10.3e,%10.3e)\n" , i, ps_new[i].x,ps_new[i].y,ps_new[i].z,ps_new[i].w, ps_old[i].x,ps_old[i].y,ps_old[i].z, vs[i].x,vs[i].y,vs[i].z, dvs[i].x,dvs[i].y,dvs[i].z, fs[i].x,fs[i].y,fs[i].z ); }}
    //const 
    float4 pi     = ps_new[iG];

    float3       v_new  = (pi.xyz - ps_old[iG].xyz)/dt;        // Leap-Frog: v_{k+1/2}  = ( p_{k+1  } - p_k ) /   dt
    float3       v_new_ = vs[iG].xyz + fs[iG].xyz * (dt*pi.w); // Leap-Frog: v_{k+1/2}' =   v_{k-1/2} + f_k / m   dt
    { pi.z = 0.0f; v_new.z = 0.0f; v_new_.z = 0.0f; } // DEBUG 2D
    //dvs   [iG] += (float4){ (v_new-v_new_), 0.0f };  // we accumulate impulses due to velocity correction which can be used to correct momentum conservation violated by the constraint solver
    vs    [iG]  = (float4){ v_new       , 0.0f };    
    ps_old[iG]  = (float4){ pi.xyz      , pi.w };
    //fs  [iG] = (float4){  (v_new_-v_new)*pi.w/dt    ,0.0f};
    
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

__kernel void  evalTrussBondForce_loc(
    const int4 ns, 
    __global const float4*  points,    // x,y,z,mass
    __global const int*     nPgroup,   // number of points in the group
    __global const int*     ips,       // indexes of points in the group
    __global       float4*  bforces,   // bond forces
    __global const int2*    bonds,     // indexes of neighbor (point,bond), if neighs[i].x == -1 it is not connected
    __global const float4*  bparams    // l0, kPress, kPull, damping
){
    local float4 points_loc[64];
    const int iG = get_global_id(0);
    const int iL = get_local_id(0);
    const int nL = get_local_size(0);

    // --- preload points to local memory
    //     NOTE: we will read less than 2*nL points ( nL is number of bonds in the group (one bond per thread)
    const int npg = nPgroup[ get_group_id(0) ];
    // 1st point
    if( iL < npg ){
        const int ip   = ips[iL];
        points_loc[iL] = points[ip];
    }
    // 2nd point
    const int i2  = nL+iL;
    if( i2 < npg ){
        const int ip   = ips[i2];
        points_loc[i2] = points[ip];
    }
    barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL]

    //if(iG==0){ printf("GPU::evalTrussBondForce(nBonds=%i,nNeigh=%i)\n", ns.x, ns.y ); }
    int2   b   = bonds  [iG];
    float4 par = bparams[iG];
    float3 d   = points_loc[b.y].xyz - points_loc[b.x].xyz;
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

// The purpose of this kernel is towfold:
// 1) improve numerical stability in single-precision by avoiding subtractions of large numbers (e.g. d = points[b.y].xyz - points[b.x].xyz; dl = length(d) - l0; )
// 2) to improve performance avoding sqrt() and division ( this is perhaps negligible considering that the kernel is memory bound )
// see: https://github.com/ProkopHapala/SimpleSimulationEngine/wiki/Truss-Simulation
__kernel void  evalTrussBondForceLinearized(
    const int4 ns, 
    __global const float4*  dpoints,   // x,y,z,mass - (small) displacements of points from initial position (around which we linearize)
    __global       float4*  bforces,   // bond forces
    __global       float4*  bvecs,     // (x,y,z,dl0) - normalized bond vectors, dl0 is displacement from initial length with respect to the rest-length
    __global const int2*    bonds,     // indexes of neighbor (point,bond), if neighs[i].x == -1 it is not connected
    __global const float4*  bparams    // l0, kPress, kPull, damping
){
    const int iG = get_global_id(0);
    //if(iG==0){ printf("GPU::evalTrussBondForce(nBonds=%i,nNeigh=%i)\n", ns.x, ns.y ); }
    const int2   b   = bonds  [iG];
    const float4 par = bparams[iG];
    const float4 hb  = bforces[iG]; 
    const float3 d   = dpoints[b.y].xyz - dpoints[b.x].xyz;  
    const float  dl  = dot(d,hb.xyz) + hb.w;
    const float  k   = 1e+6f;
    float4 f = (float4){
        hb.xyz*(k*dl), // force 
        k*dl*dl*0.5f   // potential energy
    };
    bforces[iG] = f; // we may need to do += in future 
    //bforces[iG] = (float4)( iG ,0.0,1.0,2.0);
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
    float3 d   = points [b.y].xyz - points[b.x].xyz;  // bottleneck is probably here. If we preload points to local memory it can be faster.
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



__kernel void dot_product_loc(
    const int ntot,
    __global const float* x, 
    __global const float* y, 
    __global float* result
    ) {

    // ntot = nG * nW
    // Algorithm:
    // 1. each processor 

    //__local float LOC[64];  // Local memory for partial sums
    __local float LOC[256];   // Local memory for partial sums

    const int iG = get_global_id(0);
    const int iL = get_local_id(0);
    const int iW = get_group_id(0);

    const int nG = get_global_size(0);
    const int nL = get_local_size(0);
    //const int nW = get_group_size(0);

    // Perform the dot product for a chunk of data
    float sum = 0.0f;
    for (int i=iG; i<ntot; i+=nG ) { sum += x[i] * y[i]; }
    LOC[iL] = sum;   // Store the partial sum in local memory
    barrier(CLK_LOCAL_MEM_FENCE);

    // Perform reduction within the work-group
    for (int offset = nL/2; offset>0; offset >>= 1) {
        if ( iL < offset) { LOC[iL] += LOC[iL+offset];   }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Store the result of this work-group in global memory
    if (iL==0) {  result[iW] = LOC[0]; }
}


__kernel void  dot_mat_vec_loc(
    const int4 ns, 
    __global const float*   Amat,    // [n,nNeighMax] sparse Lmat coefs at postions of neighs
    __global const float4*  xvec,    // [n,m]         right-hand-side of linear system y = A*x
    __global       float4*  yvec     // [n,m]         solution        of linear system y = A*x  
){
    __local float4 x_loc[32];
    __local float  A_loc[32];
    const int iG = get_global_id(0);
    const int iL = get_local_id(0);
    const int nL = get_local_size(0);
    if(iG>=ns.x) return;
    float4 sum = (float4){0.0f,0.0f,0.0f,0.0f};
    const int j0 = iG*ns.x;
    for (int i0=0; i0<ns.x; i0+= nL ){ 
        x_loc[iL] = xvec[   i0+iL];
        A_loc[iL] = Amat[j0+i0+iL]; 
        barrier(CLK_LOCAL_MEM_FENCE);  // wait for loading all  pos_shared[iL] 
        for (int j=0; j<nL; j++){
            sum += x_loc[j]*A_loc[j];
        }
        barrier(CLK_LOCAL_MEM_FENCE);  // block writing new     pos_shared[iL] before inner loop finished 
    }    
    yvec[iG] = sum; // we may need to do += in future 
}


__kernel void  dot_mat_vec_sparse(
    const int4 ns, 
    __global const float*   Amat,    // [n,nNeighMax] sparse Lmat coefs at postions of neighs
    __global const int*     neighs,  // [n,nNeighMax] neighbor indexes
    __global const float4*  xvec,    // [n,m]         solution        of linear system A*x=b
    __global       float4*  yvec     // [n,m]         right-hand-side of linear system A*x=b
){
    const int iG = get_global_id(0);
    if(iG>=ns.x) return;
    const int nNeighMax = ns.y;
    float4 sum = (float4){0.0f,0.0f,0.0f,0.0f};
    const int k0 = iG*nNeighMax;
    const int j0 = iG*ns.x;
    for (int k=0; k<nNeighMax; k++){
        const int j = neighs[j0+k];
        if(j<0)break;
        sum += xvec[j]*Amat[j0+j];
    }
    yvec[iG] = sum; // we may need to do += in future 
}

__kernel void  fwd_subs(
    const int4 ns, 
    __global const float*   Lmat,    // [n,nNeighMax] sparse Lmat coefs at postions of neighs
    __global const int*     neighs,  // [n,nNeighMax] neighbor indexes
    __global const float4*  bvec,    // [n,m]         right-hand-side of linear system A*x=b
    __global       float4*  xvec    // [n,m]         solution        of linear system A*x=b
    //__global const float4*  
){
    const int iG = get_global_id(0);
    if(iG>=ns.x) return;
    const int nNeighMax = ns.y;
    float4 sum = bvec[iG];
    const int k0 = iG*nNeighMax;
    const int j0 = iG*ns.x;
    for (int k=0; k<nNeighMax; k++){
        const int j = neighs[j0+k];
        if(j<0)break;
        if (j<iG){ sum += xvec[j]*Lmat[j0+j]; }
    }
    // PARALELIZATION PROBLEM !!!!   depends on previous x-values
    xvec[iG] = sum; // we may need to do += in future 
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
    //for(int ij=0; ij<8; ij++){     // change from 23 to 8 to test if it is significantly faster. It is not. The bottleneck is elsewhere.
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