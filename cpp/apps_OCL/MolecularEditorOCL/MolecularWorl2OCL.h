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

    int nAtoms = 0;
    int nGrid  = 0;
    
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
    
    void prepareBuffers( int nAtoms, int nGrid ){
    
        id_atoms      = cl->newBuffer( "atoms",      nAtoms*8, sizeof(float), NULL, CL_MEM_READ_ONLY  );
        id_gridPoints = cl->newBuffer( "gridPoints", nGrid*4,  sizeof(float), NULL, CL_MEM_READ_ONLY );      
        //id_FFPauli    = cl->newBuffer( "FFPauli",    nGrid*4, sizeof(float), NULL, CL_MEM_READ_WRITE );
        //id_FFLondon   = cl->newBuffer( "FFLondon",   nGrid*4, sizeof(float), NULL, CL_MEM_READ_WRITE );
        //id_FFelec     = cl->newBuffer( "FFelec",     nGrid*4, sizeof(float), NULL, CL_MEM_READ_WRITE );
        id_FFPauli    = cl->newBuffer( "FFPauli",    nGrid*4, sizeof(float), NULL, CL_MEM_WRITE_ONLY );
        id_FFLondon   = cl->newBuffer( "FFLondon",   nGrid*4, sizeof(float), NULL, CL_MEM_WRITE_ONLY );
        id_FFelec     = cl->newBuffer( "FFelec",     nGrid*4, sizeof(float), NULL, CL_MEM_WRITE_ONLY );
        
        task_FFPLE->global[0] = nGrid; 
        task_FFPLE->args = { INTarg(nAtoms), BUFFarg(id_atoms), BUFFarg(id_gridPoints), BUFFarg(id_FFPauli), BUFFarg(id_FFLondon), BUFFarg(id_FFelec) };
        task_FFPLE->print_arg_list();
    }
    
    void uploadAtoms(int n, Vec3d* apos, Vec3d* REQs ){
        if( nAtoms!=n ){ printf("ERROR: GridFF_OCL::uploadAtoms() Wrong Number  of Atoms \n"); exit(0); }
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
        float * buff = new float[n*8];
        
        if(FFPauli ){ 
            cl->download( id_FFPauli, buff );
            float4ToVec3d( n, buff, FFPauli );
        }
        if(FFLondon){
            cl->download( id_FFLondon, buff );
            float4ToVec3d( n, buff, FFLondon );
        }
        if(FFelec  ){ 
            cl->download( id_FFelec , buff );
            float4ToVec3d( n, buff, FFelec );
        }
        delete [] buff;
    }
    
    void evalGridFFs( GridFF& gridFF, Vec3i& nPBC ){
        prepareBuffers( gridFF.natoms, gridFF.grid.getNtot() );
        uploadAtoms( gridFF.natoms, gridFF.apos, gridFF.aREQs );
        task_FFPLE->enque();
        downloadFF( gridFF.grid.getNtot(), gridFF.FFPauli, gridFF.FFLondon, gridFF.FFelec );
        clFinish(cl->commands);
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
