#ifndef  OCL_Orb_h
#define  OCL_Orb_h

#include "datatypes.h"  
#include "datatypes2cl.h"
#include "OCL.h"
#include "Vec3.h"
//#include "Mat3.h"
#include "OrbSim_f.h"


size_t roundUp( size_t size, size_t alignment ) { 
    return ((size + alignment - 1) / alignment) * alignment;  // this is slower, but works for any alignment
    //return (size + alignment - 1) & ~(alignment - 1);         // this is faster but works only for 2^n
}


//=======================================================================
class OCL_Orb: public OCLsystem, public OrbSim_f { public:
    cl_program program1=0;

    // ===== From OrbSim_f
    //int nPoint=0, nNeighMax=0, nNeighTot=0; 
    //float4* points=0;
    //float4* forces=0;
    //float4* params=0;
    //int*    neighs=0;

    //  ----  OpenCL buffers and textures ids
    int ibuff_points=-1,ibuff_forces=-1,ibuff_vels=-1,ibuff_params=-1,ibuff_neighs=-1,ibuff_neighBs=-1,ibuff_bparams=-1, ibuff_neighB2s=-1, ibuff_bforces=-1, ibuff_bonds=-1;
    int ibuff_ps1=-1, ibuff_ps2=-1, ibuff_dps=-1, ibuff_impuls=-1, ibuff_bvec=-1, ibuff_kngs;

    int ibuff_Amat=-1;
    int ibuff_xvec=-1;
    int ibuff_yvec=-1;

    
    
    // OpenCL buffers and textures ids
    //int itex_FE_Paul=-1;
    // --- OpenCL tasks
    OCLtask* task_evalTrussForce2=0;
    OCLtask* task_evalTrussForce1=0;
    OCLtask* task_move           =0;   
    OCLtask* task_test_enque     =0; 
    OCLtask* task_blur           =0;
    OCLtask* task_assembleAndMove    =0;
    OCLtask* task_evalTrussBondForce =0;

    OCLtask* task_updateJacobi_neighs=0;

    // ------- Parameters
    int4   nDOFs    {0,0,0,0};       // number of DOFs (nPoints,nNeighMax,0,0 )
    Quat4f MDpars   {1e-3,1e-4,0,0}; // MD parameters (dt, damping, cv, cf    )
    Vec2f  bmix     {1.0,0.7};       // momentum-mixing parameters { p_new, dps_k-1 }
   
    int iker_dot_mat_vec_loc    = -1; //newKernel( "dot_mat_vec_loc",    program1 );
    int iker_dot_mat_vec_sparse = -1; //newKernel( "dot_mat_vec_sparse", program1 );


    int iker_getTrussForces      = -1;
    int iker_updateJacobi_neighs = -1;
    int iker_updateJacobi_mix    = -1;
    int iker_PD_perdictor        = -1;
    int iker_PD_corrector        = -1;

    int iker_updatePD_RHS        = -1;
    int iker_updateJacobi_lin    = -1;


    // ====================== Functions

    void makeKrenels_Orb( const char*  cl_src_dir ){
        //printf( "makeKrenels_Orb() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/spacecraft.cl", cl_src_dir );   
        printf( "makeKrenels_Orb(%s) \n", srcpath );  
        buildProgram( srcpath, program1 );
        //newTask( "getNonBond"        ,program1, 2);
        newTask( "evalTrussForce1"     ,program1, 1);
        newTask( "evalTrussForce2"     ,program1, 1);
        newTask( "move"                ,program1, 1);
        newTask( "assembleAndMove"     ,program1, 1);
        newTask( "evalTrussBondForce"  ,program1, 1);
        newTask( "test_enque"          ,program1, 1);
        newTask( "test_blur"           ,program1, 1);


        iker_dot_mat_vec_loc    = newKernel( "dot_mat_vec_loc",    program1 );
        iker_dot_mat_vec_sparse = newKernel( "dot_mat_vec_sparse", program1 );

        iker_getTrussForces      = newKernel( "getTrussForces",      program1 );
        iker_updateJacobi_neighs = newKernel( "updateJacobi_neighs", program1 );
        iker_updateJacobi_mix    = newKernel( "updateJacobi_mix",    program1 );
        iker_PD_perdictor        = newKernel( "PD_perdictor",        program1 );
        iker_PD_corrector        = newKernel( "PD_corrector",        program1 );
        iker_updatePD_RHS        = newKernel( "updatePD_RHS",        program1 );
        iker_updateJacobi_lin    = newKernel( "updateJacobi_lin",    program1 );


        //OCLtask* task = new OCLtask( this, iker, dim, global, local );
        printf( "... makeKrenels_Orb() DONE \n" );
    }

    int initCLBuffsOrb(){
        //int err=0;
        printf( "initAtomsForces() nPoint %i nNeighMax %i \n", nPoint, nNeighMax );
        ibuff_points  = newBuffer( "atoms",    nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); // points state
        ibuff_ps1     = newBuffer( "ps1",      nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); // working copy of points for solver
        ibuff_ps2     = newBuffer( "ps2",      nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); // working copy of points for solver
        ibuff_dps     = newBuffer( "ps2",      nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); // momentum (change) of points in solver from previous iteration
        ibuff_impuls  = newBuffer( "impuls",   nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); // accumulated impulses corrections

        ibuff_bvec    = newBuffer( "bvec",     nBonds   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE );
        
        ibuff_forces  = newBuffer( "forces",   nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); 
        ibuff_vels    = newBuffer( "vels",     nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); 

        ibuff_neighs  = newBuffer( "neighs",   nNeighTot, sizeof(int)   , 0, CL_MEM_READ_ONLY  ); 
        ibuff_neighBs = newBuffer( "bneighBs", nNeighTot, sizeof(int2)  , 0, CL_MEM_READ_ONLY  );
        ibuff_neighB2s= newBuffer( "bneighB2s",nNeighTot, sizeof(int)   , 0, CL_MEM_READ_ONLY  );
        ibuff_params  = newBuffer( "params",   nNeighTot, sizeof(Quat4f), 0, CL_MEM_READ_ONLY  ); 
        ibuff_kngs    = newBuffer( "kngs",     nNeighTot, sizeof(float) , 0, CL_MEM_READ_ONLY  );

        ibuff_bparams = newBuffer( "bparams",  nBonds   , sizeof(Quat4f), 0, CL_MEM_READ_ONLY  );
        ibuff_bforces = newBuffer( "bforces",  nBonds   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE );
        ibuff_bonds   = newBuffer( "bonds",    nBonds   , sizeof(int2)  , 0, CL_MEM_READ_ONLY  );
        return ibuff_points;
    }

    int initCLBuffs_CG(){
        //int err=0;
        printf( "initCLBuffs_CG() nPoint %i nNeighMax %i \n", nPoint, nNeighMax );
        ibuff_Amat = newBuffer( "Amat",  nPoint*nPoint   , sizeof(float), 0, CL_MEM_READ_WRITE );  
        ibuff_xvec = newBuffer( "xvec",   nPoint, sizeof(Quat4f) , 0, CL_MEM_READ_WRITE );      // actually xvec=points 
        ibuff_yvec = newBuffer( "yvec",   nPoint, sizeof(Quat4f) , 0, CL_MEM_READ_WRITE );
        //ibuff_neighs  = newBuffer( "neighs",   nNeighTot, sizeof(int)   , 0, CL_MEM_READ_ONLY  );  // already allocated 
        return ibuff_points;
    }

    void run_ocl( int niter, int upload_mask=0b001, int download_mask=0b001 ){
        //printf( "# ============ OCL_Orb::run_ocl() \n" );
        int err=0;
        MDpars = Quat4f{ dt, 1-damping, cv, cf };
        if(upload_mask&0b001){ err=upload( ibuff_points, points ); OCL_checkError(err, "run_ocl.1"); }
        if(upload_mask&0b010){ err=upload( ibuff_vels  , vel    ); OCL_checkError(err, "run_ocl.2"); }
        if(upload_mask&0b100){ err=upload( ibuff_forces, forces ); OCL_checkError(err, "run_ocl.3"); }
        for(int itr=0; itr<niter; itr++){
            err= task_evalTrussForce2 ->enque_raw();     //OCL_checkError(err, "run_ocl.4");
            //err= task_evalTrussForce1 ->enque_raw();   //OCL_checkError(err, "run_ocl.4");
            err= task_move            ->enque_raw();     //OCL_checkError(err, "run_ocl.4");
            //err= task_evalTrussBondForce->enque_raw();  OCL_checkError(err, "run_ocl.4");
            //err= task_assembleAndMove   ->enque_raw();  OCL_checkError(err, "run_ocl.4");
        }
        if(download_mask&0b001){err=download( ibuff_points, points ); OCL_checkError(err, "run_ocl.5"); }
        if(download_mask&0b010){err=download( ibuff_vels  , vel    ); OCL_checkError(err, "run_ocl.6"); }
        if(download_mask&0b100){err=download( ibuff_forces, forces ); OCL_checkError(err, "run_ocl.7"); }
        OCL_checkError(err, "run_ocl");
        finishRaw();
        //for(int i=0; i<4; i++){ printf( "forces[%i] (%g,%g,%g) \n", i, forces[i].x, forces[i].y, forces[i].z ); }
        //exit(0);
    }

    // __kernel void getTrussForces( 
    //     int npoint,                     // 1 number of points
    //     int nmax_neigh,                 // 2 max number of neighbors
    //     __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    //     __global       float4*  forces, // 3 [npoint] x,y,z,mass
    //     __global const int*     neighs, // 5 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    //     __global const float4*  params, // 6 [npoint,nmax_neigh] {l0, kPress, kPull, damping} 
    //     float inv_dt2                   // 7 1/dt^2, controls scale of inertial term
    // ){
    void run_getTrussForces( int ibuff_ps, int ibuff_fs ){
        //printf( "run_getTrussForces() nPoint %i \n", nPoint );
        int err = 0;
        cl_kernel ker = kernels[iker_getTrussForces];
        err |= clSetKernelArg(ker, 0, sizeof(int),    &nPoint);                               OCL_checkError(err, "run_getTrussForces.1");
        err |= clSetKernelArg(ker, 1, sizeof(int),    &nNeighMax);                            OCL_checkError(err, "run_getTrussForces.2");
        err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[ibuff_ps    ].p_gpu));        OCL_checkError(err, "run_getTrussForces.3");
        err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[ibuff_fs    ].p_gpu));        OCL_checkError(err, "run_getTrussForces.4");
        err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_neighs].p_gpu));        OCL_checkError(err, "run_getTrussForces.5");
        err |= clSetKernelArg(ker, 5, sizeof(cl_mem), &(buffers[ibuff_params].p_gpu));        OCL_checkError(err, "run_getTrussForces.6");
        float inv_dt2 = 1.0f / (dt * dt);
        err |= clSetKernelArg(ker, 6, sizeof(float), &inv_dt2);                               OCL_checkError(err, "run_getTrussForces.7");
        size_t local_work_size  = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);  OCL_checkError(err, "run_getTrussForces.8");
        OCL_checkError(err, "run_getTrussForces");
    }


    void run_updatePD_RHS( int ps ) {
        int err = 0;
        cl_kernel ker = kernels[iker_updatePD_RHS];
        err |= clSetKernelArg(ker, 0, sizeof(int),    &nPoint);                               OCL_checkError(err, "run_updatePD_RHS.1");
        err |= clSetKernelArg(ker, 1, sizeof(int),    &nNeighMax);                            OCL_checkError(err, "run_updatePD_RHS.2");
        err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[ps          ].p_gpu));        OCL_checkError(err, "run_updatePD_RHS.3");
        err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[ibuff_bvec  ].p_gpu));        OCL_checkError(err, "run_updatePD_RHS.4");
        err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_neighs].p_gpu));        OCL_checkError(err, "run_updatePD_RHS.5");
        err |= clSetKernelArg(ker, 5, sizeof(cl_mem), &(buffers[ibuff_params].p_gpu));        OCL_checkError(err, "run_updatePD_RHS.6");
        float inv_dt2 = 1.0f / (dt * dt);
        err |= clSetKernelArg(ker, 6, sizeof(float), &inv_dt2);                               OCL_checkError(err, "run_updatePD_RHS.7");
        size_t local_work_size  = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);  OCL_checkError(err, "run_updatePD_RHS.8");
        OCL_checkError(err, "run_updatePD_RHS");
    }

    void run_updateJacobi_lin(int ps_in, int ps_out) {
        int err = 0;
        cl_kernel ker = kernels[iker_updateJacobi_lin];
        err |= clSetKernelArg(ker, 0, sizeof(int),    &nPoint);                               OCL_checkError(err, "run_updateJacobi_lin.1");
        err |= clSetKernelArg(ker, 1, sizeof(int),    &nNeighMax);                            OCL_checkError(err, "run_updateJacobi_lin.2");
        err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[      ps_in ].p_gpu));        OCL_checkError(err, "run_updateJacobi_lin.3");
        err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[      ps_out].p_gpu));        OCL_checkError(err, "run_updateJacobi_lin.4");
        err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_bvec  ].p_gpu));        OCL_checkError(err, "run_updateJacobi_lin.5");
        err |= clSetKernelArg(ker, 5, sizeof(cl_mem), &(buffers[ibuff_neighs].p_gpu));        OCL_checkError(err, "run_updateJacobi_lin.6");
        err |= clSetKernelArg(ker, 6, sizeof(cl_mem), &(buffers[ibuff_kngs  ].p_gpu));        OCL_checkError(err, "run_updateJacobi_lin.7");
        size_t local_work_size  = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);  OCL_checkError(err, "run_updateJacobi_lin.8");
        OCL_checkError(err, "run_updateJacobi_lin");
    }

    // __kernel void updateJacobi_neighs( 
    //     int npoint,                     // 1 number of points
    //     int nmax_neigh,                 // 2 max number of neighbors
    //     __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    //     __global       float4*  ps_out, // 4 [npoint] x,y,z,mass
    //     __global const int*     neighs, // 5 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    //     __global const float4*  params, // 6 [npoint,nmax_neigh] {l0, kPress, kPull, damping} 
    //     float inv_dt2                   // 7 1/dt^2, controls scale of inertial term
    // ){
    void run_updateJacobi_neighs( int ps_in, int ps_out ){
        //printf( "run_updateJacobi_neighs() nPoint %i \n", nPoint );
        int err = 0;
        cl_kernel ker = kernels[iker_updateJacobi_neighs];
        err |= clSetKernelArg(ker, 0, sizeof(int),    &nPoint);                               OCL_checkError(err, "run_updateJacobi_neighs.1");   
        err |= clSetKernelArg(ker, 1, sizeof(int),    &nNeighMax);                            OCL_checkError(err, "run_updateJacobi_neighs.2");       
        err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[      ps_in ].p_gpu));        OCL_checkError(err, "run_updateJacobi_neighs.3");
        err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[      ps_out].p_gpu));        OCL_checkError(err, "run_updateJacobi_neighs.4");
        err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_neighs].p_gpu));        OCL_checkError(err, "run_updateJacobi_neighs.5");
        err |= clSetKernelArg(ker, 5, sizeof(cl_mem), &(buffers[ibuff_params].p_gpu));        OCL_checkError(err, "run_updateJacobi_neighs.6");
        float inv_dt2 = 1.0f / (dt * dt);
        err |= clSetKernelArg(ker, 6, sizeof(float), &inv_dt2);
        size_t local_work_size  = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);  OCL_checkError(err, "run_updateJacobi_neighs.7");
        OCL_checkError(err, "run_updateJacobi_neighs");
    }

    // __kernel void updateJacobi_mix( 
    //     int npoint,                     // 1 number of points
    //     int nmax_neigh,                 // 2 max number of neighbors
    //     __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    //     __global       float4*  ps_out, // 4 [npoint] x,y,z,mass
    //     __global       float4*  dps,    // 4 [npoint] x,y,z,mass momentum
    //     __global const int*     neighs, // 5 [npoint,nmax_neigh] indexes of neighbor points, if neighs[i] == -1 it is not connected, includes both bonds and collisions
    //     __global const float4*  params, // 6 [npoint,nmax_neigh] {l0, kPress, kPull, damping} 
    //     float inv_dt2,                  // 7 1/dt^2, controls scale of inertial term
    //     float2 bmix                     // 8 mixing parameter for momentum
    // ){
    void run_updateJacobi_mix(int ps_in, int ps_out ) {
        int err = 0;
        cl_kernel ker = kernels[iker_updateJacobi_mix];
        err |= clSetKernelArg(ker, 0, sizeof(int),    &nPoint);                               OCL_checkError(err, "run_updateJacobi_mix.1");
        err |= clSetKernelArg(ker, 1, sizeof(int),    &nNeighMax);                            OCL_checkError(err, "run_updateJacobi_mix.2");
        err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[      ps_in ].p_gpu));        OCL_checkError(err, "run_updateJacobi_mix.3");
        err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[      ps_out].p_gpu));        OCL_checkError(err, "run_updateJacobi_mix.4");
        err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_dps   ].p_gpu));        OCL_checkError(err, "run_updateJacobi_mix.5");
        err |= clSetKernelArg(ker, 5, sizeof(cl_mem), &(buffers[ibuff_neighs].p_gpu));        OCL_checkError(err, "run_updateJacobi_mix.6");
        err |= clSetKernelArg(ker, 6, sizeof(cl_mem), &(buffers[ibuff_params].p_gpu));        OCL_checkError(err, "run_updateJacobi_mix.7");
        float inv_dt2 = 1.0f / (dt * dt);
        err |= clSetKernelArg(ker, 7, sizeof(float), &inv_dt2);                               OCL_checkError(err, "run_updateJacobi_mix.8");
        err |= clSetKernelArg(ker, 8, sizeof(float2), &bmix);                                 OCL_checkError(err, "run_updateJacobi_mix.9");
        size_t local_work_size  = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);  OCL_checkError(err, "run_updateJacobi_mix.10");
        OCL_checkError(err, "run_updateJacobi_mix");
    }

    double run_updateJacobi_smart( int psa, int psb, int itr ) {
        //run_updateJacobi_neighs(psa, psb);
        // I want some function of form     A/( B + itr ) which interpolates between  bmix_end and bmix_start when itr transitions from nitr_end to nitr_start
        const float bmix_end   = 0.75f;
        const float bmix_start = 0.55f;
        const int   nitr_start = 3;
        const int   nitr_end   = 10; 
        // float A = (bmix_end * nitr_end - bmix_start * nitr_start) / (bmix_end - bmix_start);
        // float B = (nitr_end - nitr_start) / (bmix_end - bmix_start) - nitr_start;
        // bmix.y = A / (B + itr);
        // bmix.y = std::max(bmix_end, std::min(bmix_start, bmix.y));
        //printf( "run_updateJacobi_smart() itr: %i bmix: %g \n", itr, bmix.y );
        if(itr<nitr_start){ run_updateJacobi_neighs(psa, psb);  }
        else              { 
            if ( itr > nitr_end ) { bmix.y = bmix_end; }
            else                  { bmix.y = bmix_start + (bmix_end - bmix_start ) * (itr - nitr_start) / (nitr_end - nitr_start); }
            run_updateJacobi_mix(psa, psb); 
        }
        
       return bmix.y;
    }


    // __kernel void PD_perdictor( 
    //     int npoint,                     // 1 number of points
    //     __global const float4*  ps,     // 3 [npoint] x,y,z,mass
    //     __global       float4*  ps_out, // 4 [npoint] x,y,z,mass
    //     __global const float4*  fs    , // 4 [npoint] x,y,z,E      force 
    //     __global const float4*  vs    , // 4 [npoint] x,y,z,?      velocity 
    //     float dt    
    // ){   
    void run_PD_perdictor( int ps ){ 
        //printf( "run_PD_perdictor() nPoint %i \n", nPoint );
        int err = 0;
        cl_kernel ker = kernels[iker_PD_perdictor];
        err |= clSetKernelArg(ker, 0, sizeof(int),    &nPoint);                             OCL_checkError(err, "run_PD_perdictor().1");
        err |= clSetKernelArg(ker, 1, sizeof(cl_mem), &(buffers[ibuff_points].p_gpu));      OCL_checkError(err, "run_PD_perdictor().2");
        err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[      ps    ].p_gpu));      OCL_checkError(err, "run_PD_perdictor().3");
        err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[ibuff_forces].p_gpu));      OCL_checkError(err, "run_PD_perdictor().4");
        err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_vels  ].p_gpu));      OCL_checkError(err, "run_PD_perdictor().5");
        err |= clSetKernelArg(ker, 5, sizeof(float), &dt);                                  OCL_checkError(err, "run_PD_perdictor().6");
        size_t local_work_size  = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL); OCL_checkError(err, "run_PD_perdictor().7");
        OCL_checkError(err, "run_PD_perdictor");
    }         


    // __kernel void PD_corrector( 
    //     int npoint,                     // 1 number of points
    //     __global const float4*  ps_new, // 3 [npoint] x,y,z,mass
    //     __global       float4*  ps_old, // 4 [npoint] x,y,z,mass
    //     __global       float4*  vs    , // 4 [npoint] x,y,z,?      velocity 
    //     float dt                  
    // ){
    void run_PD_corrector( int ps ){ 
        //printf( "run_PD_corrector() nPoint %i \n", nPoint );
        int err = 0;
        cl_kernel ker = kernels[iker_PD_corrector];
        err |= clSetKernelArg(ker, 0, sizeof(int),    &nPoint);                                OCL_checkError(err, "run_PD_corrector().1");
        err |= clSetKernelArg(ker, 1, sizeof(cl_mem), &(buffers[      ps    ].p_gpu));         OCL_checkError(err, "run_PD_corrector().2");
        err |= clSetKernelArg(ker, 2, sizeof(cl_mem), &(buffers[ibuff_points].p_gpu));         OCL_checkError(err, "run_PD_corrector().3");
        err |= clSetKernelArg(ker, 3, sizeof(cl_mem), &(buffers[ibuff_vels  ].p_gpu));         OCL_checkError(err, "run_PD_corrector().4");
        err |= clSetKernelArg(ker, 4, sizeof(cl_mem), &(buffers[ibuff_impuls].p_gpu));         OCL_checkError(err, "run_PD_corrector().5");
        err |= clSetKernelArg(ker, 5, sizeof(cl_mem), &(buffers[ibuff_forces].p_gpu));         OCL_checkError(err, "run_PD_corrector().6");
        err |= clSetKernelArg(ker, 6, sizeof(float), &dt);
        size_t local_work_size  = 32;
        size_t global_work_size = roundUp(nPoint, local_work_size);
        err |= clEnqueueNDRangeKernel(commands, ker, 1, NULL, &global_work_size, &local_work_size, 0, NULL, NULL);   OCL_checkError(err, "run_PD_corrector().7");
        OCL_checkError(err, "run_PD_corrector");
    }   

    void run_projective_dynamics( int nSolverIters, int  ialg ) {
        // 1) predict positions 
        //    - evaluate external forces ( without truss contrains, like gravity, centrifugal forces etc. )
        //    - update velocities v += f/m * dt
        //    - update positions  p_solve = p + v * dt
        // 2) Solve constrains using projective dynamics matrix ( Jacobi or Gauss-Seidel method ) with M/dt^2 intertial term
        //    - for(i=0; i<Nsolv; i++)  p_solve = solver_step( p_solve, constrains) 
        //        - we should use tick-tack method ( swaping ps1 and ps2 buffers as intput and output )
        // 3) correct postion and velocity
        //    - update velocities v = ( p_solve - p ) / dt
        //    - update positions  p = p_solve
        //printf( "OCL_Orb::run_projective_dynamics() \n" );
        int psa = ibuff_ps1, psb = ibuff_ps2;
        run_PD_perdictor(psa);
        switch(ialg){
            case 0: for (int i = 0; i < nSolverIters; i++) {  run_updateJacobi_neighs( psa, psb    ); _swap(psa, psb);  } break;
            case 1: for (int i = 0; i < nSolverIters; i++) {  run_updateJacobi_mix   ( psa, psb    ); _swap(psa, psb);  } break;
            case 2: for (int i = 0; i < nSolverIters; i++) {  run_updateJacobi_smart ( psa, psb, i ); _swap(psa, psb);  } break;
            case 3:{ run_updatePD_RHS( psa );
                    for (int i = 0; i < nSolverIters; i++) {  
                        //run_updateJacobi_neighs( psa, psb );  // DEBUG
                        run_updateJacobi_lin   ( psa, psb ); _swap(psa, psb);  
                        //exit(0);
                    } 
                   } break;
        }
        run_PD_corrector(psa);
    }

    /// this function is used for debugging and performance benchmarks, it runs iterative solver of linear system, and downloads solution and resudual forces at every step
    Vec2f run_SolverConvergence( int nSolverIters, int ialg, bool bPrint=false ) {
        //printf( "OCL_Orb::run_SolverConvergence() ialg: %i bmix( %g, %g ) \n", ialg, bmix.x, bmix.y );
        int err = 0;
        int psa = ibuff_ps1; 
        int psb = ibuff_ps2;
        for(int i=0; i<nPoint; i++){ forces[i] = Quat4fZero; }
        err |= upload( psa, points       );                OCL_checkError(err, "run_SolverConvergence.1"); 
        err |= upload( ibuff_dps, forces );                OCL_checkError(err, "run_SolverConvergence.2");
        Vec2f fe =Vec2fNAN; 
        Vec2f fe1=Vec2fNAN; 
        int i=0;
        for (i=0; i<nSolverIters; i++) { 
            switch(ialg){
                case 0: run_updateJacobi_neighs( psa, psb );    break;
                case 1: run_updateJacobi_mix   ( psa, psb );    break;
                case 2: run_updateJacobi_smart ( psa, psb, i ); break;
                //case 2: run_updateJacobi_neighs ( psa, psb ); break;
            }
            run_getTrussForces( psb, ibuff_forces );  
            err |= download( psb, points          ); OCL_checkError(err, "run_SolverConvergence.5");
            err |= download( ibuff_forces, forces ); OCL_checkError(err, "run_SolverConvergence.6");
            err |= finishRaw();                      OCL_checkError(err, "run_SolverConvergence.7");
            Vec2f fe_=Vec2fZero;
            for(int ia=0; ia<nPoint; ia++){ fe_.x += forces[ia].f.norm2(); fe_.y += forces[ia].e; }
            _swap(psa, psb);
            if(i==0){ fe1=fe_; }
            //if(bPrint){ printf( "step %3i E: %10.2e Enew/Eold: %g \n", i, fe.y, fe_.y/fe.y ); }
            fe=fe_;
        }
        if(bPrint){ printf( "RESULT OCL_Orb::run_SolverConvergence()  ialg: %i bmix( %3.2f , %3.2f ) nstep: %3i Estart: %10.2e Eend: %10.2e Eend/Estart: %g \n", ialg, bmix.x, bmix.y, nSolverIters, fe1.y, fe.y, fe.y/fe1.y ); }
        return (Vec2f)fe;
    }


    static bool isnan_(const Vec3f& v){ return isnan(v.x) || isnan(v.y) || isnan(v.z); } 

    void run_PDcl( int niter, int nSolverIters, int upload_mask=0b001, int download_mask=0b001, int ialg=2 ){
        //printf( "OCL_Orb::run_PDcl() niter=%i nSolverIters=%i \n", niter, nSolverIters );
        //for(int i=0; i<nPoint; i++){  printf( "OCL_Orb::run_PDcl() points[%i] (%10.2e,%10.2e,%10.2e|%10.2e) \n", i, points[i].x, points[i].y, points[i].z, points[i].w );  }
        int err=0;
        MDpars = Quat4f{ dt, 1-damping, cv, cf };
        // if(upload_mask&0b001){ err=upload( ibuff_points, points ); OCL_checkError(err, "run_ocl.1"); }
        // if(upload_mask&0b010){ err=upload( ibuff_vels  , vel    ); OCL_checkError(err, "run_ocl.2"); }
        // if(upload_mask&0b100){ err=upload( ibuff_forces, forces ); OCL_checkError(err, "run_ocl.3"); }
        // for(int i=0; i<nPoint; i++){ impuls[i] = Quat4fZero; }
        // if(upload_mask&0b010){ err=upload( ibuff_impuls, impuls ); OCL_checkError(err, "run_ocl.4"); }

        for(int i=0; i<nPoint; i++){ impuls[i] = Quat4fZero; }
        err=upload( ibuff_impuls, impuls ); OCL_checkError(err, "run_ocl.1"); 
        err=upload( ibuff_points, points ); OCL_checkError(err, "run_ocl.2"); 
        err=upload( ibuff_vels  , vel    ); OCL_checkError(err, "run_ocl.3"); 
        err=upload( ibuff_forces, forces ); OCL_checkError(err, "run_ocl.4"); 


        for(int itr=0; itr<niter; itr++){ 
            //run_projective_dynamics( nSolverIters, ialg ); 
            //run_projective_dynamics( nSolverIters, 0 );   // non-linear Jaccobi solver without mixing
            //run_projective_dynamics( nSolverIters, 1 );   // non-linear Jaccobi solver with mixing
            //run_projective_dynamics( nSolverIters, 2 );   // non-linear Jaccobi solver with smart mixing
            run_projective_dynamics( nSolverIters, 3 );     // linear Jaccobi solver
            //exit(0);
        }
        // if(download_mask&0b001){err=download( ibuff_points, points ); OCL_checkError(err, "run_ocl.5"); }
        // if(download_mask&0b010){err=download( ibuff_vels  , vel    ); OCL_checkError(err, "run_ocl.6"); }
        // if(download_mask&0b100){err=download( ibuff_forces, forces ); OCL_checkError(err, "run_ocl.7"); }
        // if(download_mask&0b010){err=download( ibuff_impuls, impuls ); OCL_checkError(err, "run_ocl.8"); }
        
        err=download( ibuff_points, points ); OCL_checkError(err, "run_ocl.5"); 
        err=download( ibuff_vels  , vel    ); OCL_checkError(err, "run_ocl.6"); 
        err=download( ibuff_forces, forces ); OCL_checkError(err, "run_ocl.7"); 
        err=download( ibuff_impuls, impuls ); OCL_checkError(err, "run_ocl.8"); 

        OCL_checkError(err, "run_ocl");
        finishRaw();
        Vec3f cog   = Vec3fZero;
        Vec3f vcog  = Vec3fZero;
        Vec3f dvcog = Vec3fZero; 
        Vec3f fcog  = Vec3fZero;
        Vec3f torq  = Vec3fZero;
        float mass=0;
        bool bNan=false;
        for(int i=0; i<nPoint; i++){ 
            float mi = points[i].w;
            cog .add_mul ( points[i].f, mi );
        }
        for(int i=0; i<nPoint; i++){ 
            if( isnan_(points[i].f) ){ printf( "OCL_Orb::run_PDcl() NAN points[%i] (%10.2e,%10.2e,%10.2e) \n", i, points[i].f ); bNan=true; }
            if( isnan_(vel[i].f)    ){ printf( "OCL_Orb::run_PDcl() NAN vel   [%i] (%10.2e,%10.2e,%10.2e) \n", i, vel[i].f    ); bNan=true; }
            if( isnan_(forces[i].f) ){ printf( "OCL_Orb::run_PDcl() NAN forces[%i] (%10.2e,%10.2e,%10.2e) \n", i, forces[i].f ); bNan=true; }
            float mi = points[i].w;
            mass += mi;
            vcog .add_mul  ( vel   [i].f, mi );
            dvcog.add_mul  ( impuls[i].f, mi                     );
            torq .add_cross( impuls[i].f *mi , points[i].f - cog );
            fcog .add      ( forces[i].f );
        }
        if(bNan){ printf( "OCL_Orb::run_PDcl() nan detected in points => Exit() \n" ); exit(0); }
        cog .mul( 1.0/mass );
        vcog.mul( 1.0/mass );
        dvcog.mul( 1.0/mass );
        for(int i=0; i<nPoint; i++){ 
            vel[i].f.sub( dvcog );
        }
        printf( "OCL_Orb::run_PDcl() cog(%10.2e,%10.2e,%10.2e) vcog(%10.2e,%10.2e,%10.2e) dvcog(%10.2e,%10.2e,%10.2e) fcog(%10.2e,%10.2e,%10.2e) \n", cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, dvcog.x,dvcog.y,dvcog.z, fcog.x,fcog.y,fcog.z ); 
        //for(int i=0; i<4; i++){ printf( "forces[%i] (%g,%g,%g) \n", i, forces[i].x, forces[i].y, forces[i].z ); }
        //exit(0);
    }


    // OCLtask* setup_updateJacobi_neighs(){
    //     printf("setup_updateJacobi_neighs()\n" );
    //     OCLtask*& task = task_updateJacobi_neighs;
    //     if(task==0) task = getTask("updateJacobi_neighs");
    //     task->global.x = 1; // we really need just one thread to enque list of kernels
    //     task->local.x  = 1;
    //     useKernel( task->ikernel );
    //     // ------- Maybe We do-not need to do this every frame ?
    //     int err=0;
    //     //err |= _useArg   ( nDOFs  );      //1 
    //     err |= useArgBuff( ibuff_ps1    );  //1
    //     err |= useArgBuff( ibuff_ps2    );  //2
    //     err |= useArgBuff( ibuff_neighs );  //3
    //     err |= useArgBuff( ibuff_params );  //4
    //     OCL_checkError(err, "setup_updateJacobi_neighs");
    //     return task;
    // }





    OCLtask* setup_blur(){
        printf("setup_blur()\n" );
        OCLtask*& task = task_blur;
        if(task==0) task = getTask("test_blur");
        task->global.x = 1; // we really need just one thread to enque list of kernels
        task->local.x  = 1;
        useKernel( task->ikernel );
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        //err |= _useArg   ( nDOFs  );     //1 
        err |= useArgBuff( ibuff_forces ); //1
        err |= useArgBuff( ibuff_vels    );//2
        OCL_checkError(err, "setup_blur");
        return task;
    }

    OCLtask* setup_test_enque(){
        printf("setup_test_enque()\n" );
        OCLtask*& task = task_test_enque;
        if(task==0) task = getTask("test_enque");
        task->global.x = 1; // we really need just one thread to enque list of kernels
        task->local.x  = 1;
        useKernel( task->ikernel );
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        //err |= _useArg   ( nDOFs  );     //1 
        err |= useArgBuff( ibuff_forces ); //1
        err |= useArgBuff( ibuff_vels    );//2
        OCL_checkError(err, "setup_test_enque");
        return task;
    }

    void test_enque(){
        int err=0;
        err = finishRaw(); OCL_checkError(err, "test_enque.enque");
        long t0 = getCPUticks();
        err = task_test_enque->enque_raw();      OCL_checkError(err, "test_enque.enque");
        //for(int i=0; i<100; i++){ task_blur->enque_raw(); }
        err = download( ibuff_forces, forces, 100 );   OCL_checkError(err, "test_enque.download");
        err = finishRaw();                        OCL_checkError(err, "test_enque.finish");
        double T = (getCPUticks()-t0)*1e-6; printf( "test_enque: %g ms\n", T );
        for(int i=0; i<100; i++){ printf( "out[%i] (%g,%g,%g,%g) \n", i, forces[i].x, forces[i].y, forces[i].z, forces[i].w ); }
        exit(0);
    }

    OCLtask* setup_dot_mat_vec_loc( int n, int ibuff_A, int ibuff_x, int ibuff_y ){
        printf("setup_assembleAndMove()\n" );
        //OCLtask*& task = setup_dot_mat_vec_loc;
        //if(task==0) task = getTask("assembleAndMove");
        OCLtask* task = new OCLtask();
        task->ikernel = iker_dot_mat_vec_loc;
        task->global.x = nPoint;
        //task->local.x  = 1;
        task->local.x  = 32; 
        task->roundSizes();
        Quat4i ns=Quat4iZero;
        ns.x=n; 
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        useKernel( task->ikernel );
        err |= _useArg   ( ns    );    //1  
        err |= useArgBuff( ibuff_A  ); //2  
        err |= useArgBuff( ibuff_x  ); //3
        err |= useArgBuff( ibuff_y  ); //4
        OCL_checkError(err, "setup_dot_mat_vec_loc"  );
        return task;
    }
// __kernel void  dot_mat_vec_loc(
//     const int4 ns, 
//     __global const float*   Amat,    // [n,nNeighMax] sparse Lmat coefs at postions of neighs
//     __global const float4*  xvec,    // [n,m]         right-hand-side of linear system y = A*x
//     __global       float4*  yvec     // [n,m]         solution        of linear system y = A*x  
// ){


    OCLtask* setup_dot_mat_vec_sparse(  int n, int nneigh, int ibuff_A, int ibuff_neighs, int ibuff_x, int ibuff_y ){
        printf("setup_dot_mat_vec_sparse()\n" );
        OCLtask* task = new OCLtask();
        task->ikernel = iker_dot_mat_vec_sparse;
        task->global.x = n;
        //task->local.x  = 1;
        task->local.x  = 32; 
        task->roundSizes();
        Quat4i ns=Quat4iZero;
        nDOFs.x=nPoint; 
        nDOFs.y=nNeighMax;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        useKernel( task->ikernel );
        err |= _useArg   ( nDOFs  );       //1 
        err |= useArgBuff( ibuff_A      ); //2
        err |= useArgBuff( ibuff_neighs ); //3
        err |= useArgBuff( ibuff_x      ); //4
        err |= useArgBuff( ibuff_y      ); //5
        OCL_checkError(err, "setup_dot_mat_vec_sparse"  );
        return task;
    }
// __kernel void  dot_mat_vec_sparse(
//     const int4 ns, 
//     __global const float*   Amat,    // [n,nNeighMax] sparse Lmat coefs at postions of neighs
//     __global const int*     neighs,  // [n,nNeighMax] neighbor indexes
//     __global const float4*  xvec,    // [n,m]         solution        of linear system A*x=b
//     __global       float4*  yvec     // [n,m]         right-hand-side of linear system A*x=b
// ){


    OCLtask* setup_assembleAndMove(  bool bUploadParams=true  ){
        printf("setup_assembleAndMove()\n" );
        OCLtask*& task = task_assembleAndMove;
        if(task==0) task = getTask("assembleAndMove");
        task->global.x = nPoint;
        //task->local.x  = 1;
        task->local.x  = 16; task->roundSizes();
        if(bUploadParams){
            upload( ibuff_neighB2s, neighB2s  );
            upload( ibuff_bonds,    bonds    );
            upload( ibuff_bparams,  bparams  );
            upload( ibuff_points,   points   );
            upload( ibuff_forces,   forces   ); // to make sure it is initialized
            upload( ibuff_vels,     vel      ); // to make sure it is initialized
        }
        //task->local.x  = 16;
        useKernel( task->ikernel );
        nDOFs.x=nPoint; 
        nDOFs.y=nNeighMax;
        MDpars = Quat4f{ dt, 1-damping, cv, cf };
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs  );        //1 
        err |= _useArg   ( MDpars );        //1 
        err |= useArgBuff( ibuff_points  ); //2  
        err |= useArgBuff( ibuff_vels    ); //3
        err |= useArgBuff( ibuff_forces  ); //4
        err |= useArgBuff( ibuff_neighB2s );//5
        err |= useArgBuff( ibuff_bforces ); //6
        err |= _useArg   ( accel );         //7
        err |= _useArg   ( omega );         //8
        err |= _useArg   ( rot0  );         //9
        OCL_checkError(err, "setup_assembleAndMove"  );
        return task;
    }
    // __kernel void  assembleAndMove(
    //     const int4 ns, 
    //     float4 MDpars,
    //     __global       float4* points,    // x,y,z,mass
    //     __global       float4* velocities, 
    //     __global const float4* forces,
    //     __global const int*    neighB2s,  // index of bond for each neighbor, if neighs[i]==0 it is not connected, if negative it is opposite direction
    //     __global const float4* bforces,
    //     float4 accel, // acceleration of the reference frame
    //     float4 omega, // angular velocity for simulation of rotating reference frame
    //     float4 rot0   // center of rotation
    // ){

    OCLtask* setup_evalTrussBondForce(){
        printf("setup_evalTrussBondForce()\n" );
        OCLtask*& task = task_evalTrussBondForce;
        if(task==0) task = getTask("evalTrussBondForce");
        task->global.x = nBonds;
        //task->local.x  = 1;
        task->local.x  = 16; task->roundSizes();
        //task->global.x = ((nBonds/task->local.x)+1)*task->local.x; // make sure global size is multiple of local size
        //if( (task->global.x%task->local.x!=0)||(task->global.x<nBonds) ){ printf( "ERROR in setup_evalTrussBondForce() task->global.x %i\%%i= nBond=%i \n", task->global.x, task->local.x, task->global.x%task->local.x, nBonds ); exit(0); }
        useKernel( task->ikernel );
        nDOFs.x=nBonds; 
        //nDOFs.y=nNeighMax;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs  );        //1 
        err |= useArgBuff( ibuff_points  ); //2  
        err |= useArgBuff( ibuff_bforces ); //3
        err |= useArgBuff( ibuff_bonds   ); //4
        err |= useArgBuff( ibuff_bparams ); //5
        OCL_checkError(err, "setup_evalTrussBondForce");
        return task;
    }
    // __kernel void  evalTrussBondForce(
    //     const int4 ns, 
    //     __global const float4*  points,    // x,y,z,mass
    //     __global       float4*  bforces,   // bond forces
    //     __global const int2*    bonds,     // indexes of neighbor (point,bond), if neighs[i].x == -1 it is not connected
    //     __global const float4*  bparams    // l0, kPress, kPull, damping
    // ){

    OCLtask* setup_move(){
        printf("setup_move()\n" );
        OCLtask*& task = task_move;
        if(task==0) task = getTask("move");
        task->global.x = nPoint;
        //task->local.x  = 1;
        task->local.x  = 16; task->roundSizes();
        useKernel( task->ikernel );
        nDOFs.x=nPoint; 
        //nDOFs.y=nNeighMax;
        MDpars = Quat4f{ dt, 1-damping, cv, cf };
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs  );       //1 
        err |= _useArg   ( MDpars );       //1 
        err |= useArgBuff( ibuff_points ); //2  
        err |= useArgBuff( ibuff_vels   ); //3
        err |= useArgBuff( ibuff_forces ); //4
        OCL_checkError(err, "setup_move");
        return task;
    }
    // __kernel void  move(
    //     const int4 ns, 
    //     float4 MDpars,
    //     __global       float4* points,    // x,y,z,mass
    //     __global       float4* velocities, 
    //     __global const float4* forces
    // ){

    OCLtask* setup_evalTrussForce2( bool bUploadParams=true ){
        printf("setup_evalTrussForce2() \n" );
        OCLtask*& task = task_evalTrussForce2;
        if(task==0) task = getTask("evalTrussForce2");
        if(bUploadParams){
            upload( ibuff_neighBs, neighBs  );
            upload( ibuff_bparams, bparams  );
            upload( ibuff_points,  points   );
            upload( ibuff_forces,  forces   ); // to make sure it is initialized
            upload( ibuff_vels,    vel      ); // to make sure it is initialized
        }
        task->global.x = nPoint;
        //task->local.x  = 1;
        //task->local.x  = 16;
        task->local.x  = 16; task->roundSizes();
        useKernel( task->ikernel );
        nDOFs.x=nPoint; 
        nDOFs.y=nNeighMax;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );         //1 
        err |= useArgBuff( ibuff_points  ); //2  
        err |= useArgBuff( ibuff_vels    ); //3  
        err |= useArgBuff( ibuff_forces  ); //4
        err |= useArgBuff( ibuff_neighBs ); //5
        err |= useArgBuff( ibuff_bparams ); //6
        err |= _useArg   ( accel );         //7
        err |= _useArg   ( omega );         //8
        err |= _useArg   ( rot0  );         //9 
        OCL_checkError(err, "setup_evalTrussForce2");
        return task;        
    }
    // __kernel void  evalTrussForce(
    //     const int4 ns,                    //1
    //     __global const float4*  points,   //2 x,y,z,mass
    //     __global const float4*  vels,     //3 velocities are used for damping 
    //     __global       float4*  forces,   //4
    //     __global const int*     neighs,   //5 indexes of neighbors, if neighs[i] == -1 it is not connected
    //     __global const float4*  params    //6 l0, kPress, kPull, damping
    //     float4 accel, //7 acceleration of the reference frame
    //     float4 omega, //8 angular velocity for simulation of rotating reference frame
    //     float4 rot0   //9 center of rotation
    // ){

    OCLtask* setup_evalTrussForce1( bool bUploadParams=true ){
        printf("setup_evalTrussForce1() \n" );
        OCLtask*& task = task_evalTrussForce1;
        if(task==0) task = getTask("evalTrussForce1");
        if(bUploadParams){
            upload( ibuff_neighs, neighs  );
            upload( ibuff_params, params  );
            upload( ibuff_points, points  );
            upload( ibuff_forces, forces  ); // to make sure it is initialized
            upload( ibuff_vels,   vel     ); // to make sure it is initialized
        }
        //int nloc = 64;
        task->global.x = nPoint;
        //task->local.x  = 1;
        task->local.x  = 16; task->roundSizes();
        useKernel( task->ikernel );
        nDOFs.x=nPoint; 
        nDOFs.y=nNeighMax;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );        //1 
        err |= useArgBuff( ibuff_points ); //2  
        err |= useArgBuff( ibuff_vels   ); //2  
        err |= useArgBuff( ibuff_forces ); //3
        err |= useArgBuff( ibuff_neighs ); //4
        err |= useArgBuff( ibuff_params ); //5
        err |= _useArg   ( accel );        //7
        err |= _useArg   ( omega );        //8
        err |= _useArg   ( rot0  );        //9 
        OCL_checkError(err, "setup_evalTrussForce");
        return task;        
    }
    // __kernel void  evalTrussForce(
    //     const int4 ns,                    //1
    //     __global const float4*  points,   //2 x,y,z,mass
    //     __global       float4*  forces,   //3
    //     __global const float4*  vels,     // velocities are used for damping 
    //     __global const int*     neighs,   //4 indexes of neighbors, if neighs[i] == -1 it is not connected
    //     __global const float4*  params    //5 l0, kPress, kPull, damping
    //     float4 accel, // acceleration of the reference frame
    //     float4 omega, // angular velocity for simulation of rotating reference frame
    //     float4 rot0   // center of rotation
    // ){

};

#endif
