#ifndef  OCL_Orb_h
#define  OCL_Orb_h

#include "datatypes.h"  
#include "datatypes2cl.h"
#include "OCL.h"
#include "Vec3.h"
//#include "Mat3.h"
#include "OrbSim.h"

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
    int ibuff_points=-1,ibuff_forces=-1,ibuff_vels=-1,ibuff_params=-1,ibuff_neighs=-1,ibuff_neighBs=-1,ibuff_bparams=-1;
    // OpenCL buffers and textures ids
    //int itex_FE_Paul=-1;
    // --- OpenCL tasks
    OCLtask* task_evalTrussForce2=0;
    OCLtask* task_evalTrussForce =0;
    OCLtask* task_move           =0;    

    // ------- Parameters
    int4   nDOFs    {0,0,0,0};       // number of DOFs (nPoints,nNeighMax,0,0 )
    Quat4f MDpars   {1e-3,1e-4,0,0}; // MD parameters (dt, damping, cv, cf    )
   
    // ====================== Functions

    void makeKrenels_Orb( const char*  cl_src_dir ){
        //printf( "makeKrenels_Orb() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/spacecraft.cl", cl_src_dir );   
        printf( "makeKrenels_Orb(%s) \n", srcpath );  
        buildProgram( srcpath, program1 );
        //newTask( "getNonBond"             ,program1, 2);
        newTask( "evalTrussForce"      ,program1, 1);
        newTask( "evalTrussForce2"     ,program1, 1);
        newTask( "move"                ,program1, 1);
        printf( "... makeKrenels_Orb() DONE \n" );
    }

    int initCLBuffsOrb(){
        //int err=0;
        printf( "initAtomsForces() nPoint %i nNeighMax %i \n", nPoint, nNeighMax );
        ibuff_points  = newBuffer( "atoms",    nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); 
        ibuff_forces  = newBuffer( "forces",   nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); 
        ibuff_vels    = newBuffer( "vels",     nPoint   , sizeof(Quat4f), 0, CL_MEM_READ_WRITE ); 
        ibuff_neighs  = newBuffer( "neighs",   nNeighTot, sizeof(int)   , 0, CL_MEM_READ_ONLY  ); 
        ibuff_neighBs = newBuffer( "bneighBs", nNeighTot, sizeof(int2)  , 0, CL_MEM_READ_ONLY  );
        ibuff_bparams = newBuffer( "params",   nBonds   , sizeof(Quat4f), 0, CL_MEM_READ_ONLY  );
        //ibuff_params  = newBuffer( "params", nNeighTot, sizeof(Quat4f), 0, CL_MEM_READ_ONLY  ); 
        return ibuff_points;
    }

    void run_ocl( int niter, int upload_mask=0b001, int download_mask=0b001 ){
        int err=0;
        if(upload_mask&0b001) err=upload( ibuff_points, points );
        if(upload_mask&0b010) err=upload( ibuff_vels  , vel    );
        if(upload_mask&0b100) err=upload( ibuff_forces, forces );
        for(int itr=0; itr<niter; itr++){
            err=task_evalTrussForce2->enque_raw();
            //err=task_move->enque_raw();
        }
        if(download_mask&0b001)err=download( ibuff_points, points );
        if(download_mask&0b010)err=download( ibuff_vels  , vel    );
        if(download_mask&0b100)err=download( ibuff_forces, forces );
        OCL_checkError(err, "run_ocl");
        finishRaw();
    }

    OCLtask* setup_move(){
        printf("setup_move()\n" );
        OCLtask*& task = task_move;
        if(task==0) task = getTask("move");
        task->global.x = nPoint;
        task->local.x  = 1;
        useKernel( task->ikernel );
        nDOFs.x=nPoint; 
        //nDOFs.y=nNeighMax;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs  );       //1 
        err |= _useArg   ( MDpars );       //1 
        err |= useArgBuff( ibuff_points ); //2  
        err |= useArgBuff( ibuff_vels   ); //3
        err |= useArgBuff( ibuff_forces ); //4
        OCL_checkError(err, "setup_move");
        return task;
        // __kernel void  move(
        //     const int4 ns, 
        //     float4 MDpars,
        //     __global       float4* points,    // x,y,z,mass
        //     __global       float4* velocities, 
        //     __global const float4* forces
        // ){
    }

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
        task->local.x  = 1;
        useKernel( task->ikernel );
        nDOFs.x=nPoint; 
        nDOFs.y=nNeighMax;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );         //1 
        err |= useArgBuff( ibuff_points  ); //2  
        err |= useArgBuff( ibuff_forces  ); //3
        err |= useArgBuff( ibuff_neighBs ); //4
        err |= useArgBuff( ibuff_bparams ); //5
        OCL_checkError(err, "setup_evalTrussForce2");
        return task;
        /*
        __kernel void  evalTrussForce(
            const int4 ns,                    //1
            __global const float4*  points,   //2 x,y,z,mass
            __global       float4*  forces,   //3
            __global const int*     neighs,   //4 indexes of neighbors, if neighs[i] == -1 it is not connected
            __global const float4*  params    //5 l0, kPress, kPull, damping
        ){
        */
    }

    /*
    OCLtask* setup_evalTrussForce( int n, OCLtask* task=0){
        printf("setup_evalTrussForce(n=%i) \n", n );
        if(task==0) task = getTask("evalTrussForce");
        //int nloc = 64;
        task->global.x = nPoint;
        task->local.x  = 1;
        useKernel( task->ikernel );
        nDOFs.x=nPoint; 
        nDOFs.y=nNeighMax;
        // ------- Maybe We do-not need to do this every frame ?
        int err=0;
        err |= _useArg   ( nDOFs );        //1 
        err |= useArgBuff( ibuff_points ); //2  
        err |= useArgBuff( ibuff_forces ); //3
        err |= useArgBuff( ibuff_neighs ); //4
        err |= useArgBuff( ibuff_params ); //5
        OCL_checkError(err, "setup_evalTrussForce");
        return task;
        
        __kernel void  evalTrussForce(
            const int4 ns,                    //1
            __global const float4*  points,   //2 x,y,z,mass
            __global       float4*  forces,   //3
            __global const int*     neighs,   //4 indexes of neighbors, if neighs[i] == -1 it is not connected
            __global const float4*  params    //5 l0, kPress, kPull, damping
        ){
        
    }
    */

};

#endif
