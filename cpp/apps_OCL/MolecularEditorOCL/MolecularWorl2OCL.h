#ifndef MolecularWorldOCL_h
#define MolecularWorldOCL_h

#include "OCL.h"
#include "MMFF.h" // we can perhaps make some #ifdef OCL insde MMFF.h
#include "arrayConvert.h"

/*
ToDo:
 - make evaluation of gridFF inside OpenCL kernel
*/

class GridFF_OCL{ public:

    OCLsystem* cl = 0;
    OCLtask *task_FFPLE = 0;
    
    Quat4i nGrid=Quat4iZero,npbc=Quat4iZero;
    Quat4f pos0=Quat4fZero,dA=Quat4fZero,dB=Quat4fZero,dC=Quat4fZero;
    float alpha = 0.0;

    int nAtoms = 0;
    int nGridTot  = 0;
    
    int id_FFPauli    = -1;
    int id_FFLondon   = -1;
    int id_FFelec     = -1;
    int id_gridPoints = -1;
    int id_atoms      = -1;
    
    // ==== Functions

/*
    void init(OCLsystem* cl_, char* fname ){ 
        cl=cl_; 
        int err = cl->buildProgram( fname );                              
        OCL_checkError(err, "cl.buildProgram");
    };
*/

    void init( OCLsystem* cl_, char* fname ){
        cl = cl_; 
        int err = cl->buildProgram( fname );                              
        OCL_checkError(err, "cl.buildProgram");
        
        int id_evalPLE = cl->newKernel("evalPLE"); DEBUG
        task_FFPLE = new OCLtask( cl, id_evalPLE, 1, -1, 32 ); DEBUG
    }
    
    void prepareBuffers( int nAtoms_, int nGrid_ ){
        nAtoms = nAtoms_;
        nGridTot  = nGrid_;
        DEBUG
        
        id_atoms      = cl->newBuffer( "atoms",      nAtoms*8, sizeof(float), NULL, CL_MEM_READ_ONLY  ); DEBUG
        id_gridPoints = cl->newBuffer( "gridPoints", nGridTot*4,  sizeof(float), NULL, CL_MEM_READ_ONLY );     DEBUG 
        //id_FFPauli    = cl->newBuffer( "FFPauli",    nGridTot*4, sizeof(float), NULL, CL_MEM_READ_WRITE );
        //id_FFLondon   = cl->newBuffer( "FFLondon",   nGridTot*4, sizeof(float), NULL, CL_MEM_READ_WRITE );
        //id_FFelec     = cl->newBuffer( "FFelec",     nGridTot*4, sizeof(float), NULL, CL_MEM_READ_WRITE );
        id_FFPauli    = cl->newBuffer( "FFPauli",    nGridTot*4, sizeof(float), NULL, CL_MEM_WRITE_ONLY ); DEBUG
        id_FFLondon   = cl->newBuffer( "FFLondon",   nGridTot*4, sizeof(float), NULL, CL_MEM_WRITE_ONLY ); DEBUG
        id_FFelec     = cl->newBuffer( "FFelec",     nGridTot*4, sizeof(float), NULL, CL_MEM_WRITE_ONLY ); DEBUG
        
        DEBUG
        
        //task_FFPLE->global[0] = nGrid; 
        //task_FFPLE->args = { INTarg(nAtoms), BUFFarg(id_atoms), BUFFarg(id_gridPoints), BUFFarg(id_FFPauli), BUFFarg(id_FFLondon), BUFFarg(id_FFelec) };
        //task_FFPLE->print_arg_list();
                
        DEBUG
    }
    
    
    void setupKernel( GridFF& gridFF ){
    
        task_FFPLE->global[0] = nGridTot; 
        //task_FFPLE->args = { INTarg(nAtoms), BUFFarg(id_atoms), BUFFarg(id_gridPoints), BUFFarg(id_FFPauli), BUFFarg(id_FFLondon), BUFFarg(id_FFelec) };
        //task_FFPLE->print_arg_list();
        //__kernel void evalPLE(
        //    const int   nAtoms,
        //    const uint4 nGrid,
        //    const uint4 npbc,
        //    const float4 pos0,
        //    const float4 dA,
        //    const float4 dB,
        //    const float4 dC,
        //    const float alpha,
        //    __global float8*   atoms,
        //    __global float4*   FFPauli,
        //    __global float4*   FFLondon,
        //    __global float4*   FFelec 
        //){
        
        nGrid.setXYZ( gridFF.grid.n );
        npbc  = {1,1,1,0};
        pos0.setXYZ( (Vec3f)gridFF.grid.pos0    );
        dA  .setXYZ( (Vec3f)gridFF.grid.dCell.a );
        dB  .setXYZ( (Vec3f)gridFF.grid.dCell.b );
        dC  .setXYZ( (Vec3f)gridFF.grid.dCell.c );
        alpha = gridFF.alpha;
        
        task_FFPLE->args = { 
            INTarg(nAtoms), 
            REFarg(nGrid),
            REFarg(npbc),
            REFarg(pos0),
            REFarg(dA),
            REFarg(dB),
            REFarg(dC),
            FLOATarg(alpha),
            BUFFarg(id_atoms), 
            //BUFFarg(id_gridPoints), 
            BUFFarg(id_FFPauli), 
            BUFFarg(id_FFLondon), 
            BUFFarg(id_FFelec)
        };
    }
    
    void uploadAtoms(int n, Vec3d* apos, Vec3d* REQs ){
        if( nAtoms!=n ){ printf("ERROR: GridFF_OCL::uploadAtoms() Wrong Number  of Atoms: n(%i) != nAtoms(%i) \n, ", n, nAtoms ); exit(0); }
        float * buff = new float[n*8];
        Vec3dTofloat8( n, apos, REQs, buff );
        cl->upload( id_atoms, buff );
        delete [] buff;
    }
    
    void uploadGridPos(Vec3i ns, Mat3f dCell ){
        int n = ns.x*ns.y*ns.z;
        float * buff = new float[n*4];
        int i4=0; 
        for(int iz=0; iz<ns.z; iz++){
            for(int iy=0; iy<ns.y; iy++){
                Vec3f p = dCell.b*iy + dCell.c*iz;
                for(int ix=0; ix<ns.x; ix++){
                    buff[i4+0] = p.x;
                    buff[i4+1] = p.y;
                    buff[i4+2] = p.z;
                    i4+=4;           
                    p.add( dCell.a );
                }
            }
        }
        cl->upload( id_gridPoints, buff );
        delete [] buff;
    }
    
    void downloadFF(int n, Vec3d* FFPauli, Vec3d* FFLondon, Vec3d* FFelec ){
        if( nGridTot!=n ){ printf("ERROR: GridFF_OCL::downloadFF() Wrong number of grid points: n(%i) != nGrid(%i) \n, ", n, nGridTot ); exit(0); }
        float * buff = new float[n*4];
        if(FFPauli ){  cl->download( id_FFPauli,  buff ); float4ToVec3d( n, buff, FFPauli  ); printf("FFPauli  downloaded\n"); }
        if(FFLondon){  cl->download( id_FFLondon, buff ); float4ToVec3d( n, buff, FFLondon ); printf("FFLondon downloaded\n"); }
        if(FFelec  ){  cl->download( id_FFelec ,  buff ); float4ToVec3d( n, buff, FFelec   ); printf("FFelec   downloaded\n"); }
        delete [] buff;
    }
    
    void evalGridFFs( GridFF& gridFF, const Vec3i& nPBC ){
        printf( "gridFF.natoms %i \n", gridFF.natoms );
        prepareBuffers( gridFF.natoms, gridFF.grid.getNtot() ); DEBUG
        setupKernel( gridFF );
        uploadAtoms( gridFF.natoms, gridFF.apos, gridFF.aREQs ); DEBUG
        task_FFPLE->enque(); DEBUG
        downloadFF( gridFF.grid.getNtot(), gridFF.FFPauli, gridFF.FFLondon, gridFF.FFelec ); DEBUG
        clFinish(cl->commands); DEBUG
    }
        
};

class MolecularWorldOCL{ public:


    void prepareBuffers( const MMFF& mmff ){
        //natoms=0, nbonds=0, nang=0, ntors=0;
        //Vec2i  * bond2atom = NULL;
        //double * bond_0    = NULL;  // [A]
        //double * bond_k    = NULL;  // [eV/A] ?
        //Vec2i  * ang2bond  = NULL;
        //Vec3i  * ang2atom  = NULL;
        //Vec2d  * ang_0     = NULL; // [1]
        //double * ang_k     = NULL; // [eV/A^2]
        //Vec3i  * tors2bond = NULL;
        //Quat4i * tors2atom = NULL;
        //Vec2d  * tors_0    = NULL; // [1]
        //double * tors_k    = NULL; // [eV/A^2]
    
        //cl.init();
        /*
        cl.newBuffer( "molecule",   nAtoms*8,       sizeof(float), (float*)molecule,  CL_MEM_READ_ONLY );
        cl.newBuffer( "pos",        nMols*8,        sizeof(float), (float*)pos,       CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR );
        cl.newBuffer( "force",      nMols*8,        sizeof(float), (float*)force,     CL_MEM_READ_WRITE );
        cl.newBuffer( "atomsT",     nAtoms*nMols*4, sizeof(float), (float*)atomsT_,   CL_MEM_READ_WRITE );
        */
    }
    
};

#endif
