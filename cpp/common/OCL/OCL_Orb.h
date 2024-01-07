#ifndef  OCL_Orb_h
#define  OCL_Orb_h

#include "datatypes.h"  
#include "datatypes2cl.h"
#include "OCL.h"
#include "Vec3.h"
//#include "Mat3.h"
#include "OrbSim_f.h"

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
        //newTask( "getNonBond"        ,program1, 2);
        newTask( "evalTrussForce1"     ,program1, 1);
        newTask( "evalTrussForce2"     ,program1, 1);
        newTask( "move"                ,program1, 1);
        newTask( "assembleAndMove"     ,program1, 1);
        newTask( "evalTrussBondForce"  ,program1, 1);
        newTask( "test_enque"          ,program1, 1);
        newTask( "test_blur"           ,program1, 1);
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
        ibuff_neighB2s= newBuffer( "bneighB2s",nNeighTot, sizeof(int)   , 0, CL_MEM_READ_ONLY  );
        ibuff_params  = newBuffer( "params",   nNeighTot, sizeof(Quat4f), 0, CL_MEM_READ_ONLY  ); 
        ibuff_bparams = newBuffer( "bparams",  nBonds   , sizeof(Quat4f) , 0, CL_MEM_READ_ONLY  );
        ibuff_bforces = newBuffer( "bforces",  nBonds   , sizeof(Quat4f) , 0, CL_MEM_READ_WRITE );
        ibuff_bonds   = newBuffer( "bonds",    nBonds   , sizeof(int2)   , 0, CL_MEM_READ_ONLY  );
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
            err= task_evalTrussForce2 ->enque_raw();  //OCL_checkError(err, "run_ocl.4");
            //err= task_evalTrussForce1 ->enque_raw();  //OCL_checkError(err, "run_ocl.4");
            err= task_move            ->enque_raw();  //OCL_checkError(err, "run_ocl.4");
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
