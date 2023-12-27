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

    // dimensions
    //int nPoint=0, nNeighMax=0, nNeighTot=0; 
    
    // OpenCL buffers and textures ids
    int ibuff_points=-1,ibuff_forces=-1,ibuff_params=-1,ibuff_neighs=-1;
    // OpenCL buffers and textures ids
    //int itex_FE_Paul=-1;

    int4 nDOFs    {0,0,0,0};  // number of DOFs (nPoints,nNeighMax,0,0)


    // cpu buffers    - Form OrbSim_f
    //float4* points=0;
    //float4* forces=0;
    //float4* params=0;
    //int*    neighs=0;

    

    // ------- Grid
    
    // ====================== Functions

    void makeKrenels_Orb( const char*  cl_src_dir ){
        printf( "makeKrenels_Orb() \n" );
        char srcpath[1024];
        sprintf( srcpath, "%s/spacecraft.cl", cl_src_dir );     
        buildProgram( srcpath, program1 );
        //newTask( "getNonBond"             ,program1, 2);
        newTask( "evalTrussForce"      ,program1, 1);
        printf( "... makeKrenels_Orb() DONE \n" );
    }

    int initBuffs_Orb( int nPoint_, int nNeighMax_, bool bCPU=true ){
        nPoint=nPoint_;
        nNeighMax=nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        printf( "initAtomsForces() nPoint %i nNeighMax %i \n", nPoint, nNeighMax );
        ibuff_points  = newBuffer( "atoms",  nPoint   , sizeof(float4), 0, CL_MEM_READ_WRITE ); if(bCPU)points = (float4*)malloc( nPoint*sizeof(float4)    );
        ibuff_forces  = newBuffer( "forces", nPoint   , sizeof(float4), 0, CL_MEM_READ_WRITE ); if(bCPU)forces = (float4*)malloc( nPoint*sizeof(float4)    );
        ibuff_params  = newBuffer( "forces", nNeighTot, sizeof(float4), 0, CL_MEM_READ_WRITE ); if(bCPU)params = (float4*)malloc( nNeighTot*sizeof(float4) );
        ibuff_neighs  = newBuffer( "REQs",   nNeighTot, sizeof(int)   , 0, CL_MEM_READ_ONLY  ); if(bCPU)neighs = (int*   )malloc( nNeighTot*sizeof(int)    );
        return ibuff_points;
    }

    OCLtask* setup_evalTrussForce( int n, OCLtask* task=0){
        printf("setup_evalTrussForce(n=%i) \n", n );
        if(task==0) task = getTask("evalTrussForce");
        //int nloc = 1;
        //int nloc = 2;
        //int nloc = 4;
        //int nloc = 8;
        //int nloc = 16;
        int nloc = 32;
        //int nloc = 64;
        task->global.x = nPoint;
        //task->local.x  = nPoint;
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

};

#endif
