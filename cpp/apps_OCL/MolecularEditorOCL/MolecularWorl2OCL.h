#ifndef MolecularWorldOCL_h
#define MolecularWorldOCL_h

#include <vector>

#include <CL/cl.h>

#include "OCL.h"
#include "MMFF.h" // we can perhaps make some #ifdef OCL insde MMFF.h
#include "Molecule.h"
#include "arrayConvert.h"

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

/*
ToDo:
 - make evaluation of gridFF inside OpenCL kernel
*/

// =======================
//   MolecularWorldOCL
// =======================

void frag2atoms(const Vec3f& pos, const Quat4f& qrot, int n, float8* atom0s, float8* atoms ){
    Mat3f mrot; qrot.toMatrix(mrot);
    //printf( "pos (%g,%g,%g) qrot (%g,%g,%g,%g) \n", pos.x,pos.y,pos.z,   qrot.x,qrot.y,qrot.z,qrot.w );
    for( int j=0; j<n; j++ ){
        Vec3f Mp;
        atoms[j] = atom0s[j];
        //mrot.dot_to_T( *((Vec3f*)(atom0s+j)), Mp );  // this is in  MMFF::frag2atoms
        mrot.dot_to( *((Vec3f*)(atom0s+j)), Mp );
        ((Vec3f*)(atoms+j))->set_add( pos, Mp );
        //printf( "atom %i  xyz(%g,%g,%g) REQ(%g,%g,%g) \n", j ,atoms[j].x, atoms[j].y, atoms[j].z,  atoms[j].hx, atoms[j].hy, atoms[j].hz  ); 
    }
}

class RigidMolecularWorldOCL{ public:

    OCLsystem* cl = 0;
    OCLtask *task_getForceRigidSystemSurfGrid = 0;
    OCLtask *task_getFEtot = 0;
    
    bool bGetFEtot = true;
    
    int nSystems  = 0;
    int nMols     = 0;
    int nMolTypes = 0;
    int nAtomInMolMax  = 0;
    
    //int nAtomsInTypes = 0; // derived 
    int nMolInstances = 0; // derived
   
    Quat4f pos0=Quat4fZero, dA=Quat4fZero, dB=Quat4fZero, dC=Quat4fZero, testPLQ=Quat4fZero;
    float alpha = 0.0;
    
    std::vector<float8> atomsInTypes;
    std::vector<int2>   molTypes;
    
    int2   * mol2atoms    = 0;
    float8 * poses        = 0;
    float8 * fposes       = 0; 
    
    Quat4f * poss      = 0;
    Quat4f * FEs       = 0; 
    
    //cl_mem img_FFPauli;   //  = -1;
    //cl_mem img_FFLondon;  //  = -1;
    //cl_mem img_FFelec;    //  = -1;
    
    int id_FFPauli;   //  = -1;
    int id_FFLondon;  //  = -1;
    int id_FFelec;    //  = -1;
    
    int id_mol2atoms    = -1;
    int id_atomsInTypes = -1;
    int id_poses        = -1;
    int id_fposes       = -1;
    
    int id_poss         = -1;
    int id_FEs          = -1;
    
    // ==== Functions
    
    void init( OCLsystem* cl_, char* fname ){
        cl = cl_; 
        int err = cl->buildProgram( fname );                              
        OCL_checkError(err, "cl.buildProgram");
        
        int id_evalPLE = cl->newKernel("getForceRigidSystemSurfGrid");               DEBUG;
        task_getForceRigidSystemSurfGrid = new OCLtask( cl, id_evalPLE, 1, -1, 32 ); DEBUG;
        
        if( bGetFEtot ){  
            int id_getFEtot   = cl->newKernel("getFEtot");              DEBUG;
            task_getFEtot = new OCLtask( cl, id_getFEtot, 1, -1, 1 );  DEBUG;
        }
    }
    
    void addMolType( const Molecule& molecule ){
        molTypes.push_back( (int2){ atomsInTypes.size(), molecule.natoms} );
        for( int ia=0; ia<molecule.natoms; ia++ ){
            float8 atom={0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0};
            //float8 atom = (float8){ (float4){ (float3)(Vec3f)molecule.pos [ia], 0.0},  (float3)(float4){(Vec3f)molecule.REQs[ia], 0.0}  };
            *((Vec3f*)(((float*)&atom)+0)) = (Vec3f)molecule.pos [ia];
            *((Vec3f*)(((float*)&atom)+4)) = (Vec3f)molecule.REQs[ia]; 
            atom.hy = sqrt(atom.hy); // eps_ii 
            atomsInTypes.push_back( atom );
        }
    }
    
    inline void setMolInstance( int isystem, int imol, int iType ){
        mol2atoms[ nMols*isystem + imol ] = molTypes[iType]; 
    }
    
    void updateMolStats(){
        nAtomInMolMax = 0;
        for(int i=0; i<nMolInstances; i++){
            if( nAtomInMolMax < mol2atoms[i].y ) nAtomInMolMax =  mol2atoms[i].y;
        }    
    }
    
    void prepareBuffers_getFEtot( int nPoss ){ 
        FEs  = new Quat4f[nPoss];
        poss = new Quat4f[nPoss];
        id_poss = cl->newBuffer( "poss",  nPoss, sizeof(Quat4f), NULL, CL_MEM_READ_WRITE  ); DEBUG;
        id_FEs  = cl->newBuffer( "FEs",   nPoss, sizeof(Quat4f), NULL, CL_MEM_READ_WRITE  ); DEBUG;
    };
    
    void prepareBuffers( int nSystems_, int nMols_, Vec3i nGrid, float* FFpauli, float* FFlondon, float* FFelec ){
        int err;
        nMols    = nMols_;
        nSystems = nSystems_;
        nMolInstances = nMols * nSystems;
        mol2atoms    = new int2  [nMolInstances]; 
        poses        = new float8[nMolInstances]; 
        fposes       = new float8[nMolInstances];
                
        id_mol2atoms    = cl->newBuffer( "molTypes",     nMolInstances,       sizeof(int2)  , NULL,                CL_MEM_READ_WRITE ); DEBUG;
        id_atomsInTypes = cl->newBuffer( "atomsInTypes", atomsInTypes.size(), sizeof(float8), atomsInTypes.data(), CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR ); DEBUG;
        
        id_poses        = cl->newBuffer( "poses",        nMolInstances, sizeof(float8), NULL, CL_MEM_READ_WRITE  ); DEBUG;
        id_fposes       = cl->newBuffer( "fposes",       nMolInstances, sizeof(float8), NULL, CL_MEM_READ_WRITE  ); DEBUG;
      
        cl_image_format imgFormat = (cl_image_format){CL_RGBA,CL_FLOAT};
        //img_FFPauli  = clCreateImage3D( cl->context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, &imgFormat, nGrid.x, nGrid.y, nGrid.z, 0, 0, FFpauli,  &err ); OCL_checkError(err, "clCreateImage3D img_FFPauli" );
        //img_FFLondon = clCreateImage3D( cl->context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, &imgFormat, nGrid.x, nGrid.y, nGrid.z, 0, 0, FFlondon, &err ); OCL_checkError(err, "clCreateImage3D img_FFLondon");
        //img_FFelec   = clCreateImage3D( cl->context, CL_MEM_WRITE_ONLY | CL_MEM_COPY_HOST_PTR, &imgFormat, nGrid.x, nGrid.y, nGrid.z, 0, 0, FFelec,   &err ); OCL_checkError(err, "clCreateImage3D img_FFelec"  );
    
        id_FFPauli  = cl->newBufferImage3D( "FFPauli",  nGrid.x, nGrid.y, nGrid.z, sizeof(float), FFpauli,  CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        id_FFLondon = cl->newBufferImage3D( "FFLondon", nGrid.x, nGrid.y, nGrid.z, sizeof(float), FFlondon, CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        id_FFelec   = cl->newBufferImage3D( "FFelec",   nGrid.x, nGrid.y, nGrid.z, sizeof(float), FFelec,   CL_MEM_READ_ONLY  | CL_MEM_COPY_HOST_PTR , {CL_RGBA, CL_FLOAT} );
        
        printf( " id_FFPauli, id_FFLondon, id_FFelec %i %i %i  | sizeof(cl_mem) %i \n ", id_FFPauli, id_FFLondon, id_FFelec , sizeof(cl_mem) );
    
    }
    
    void prepareBuffers( int nSystems_, int nMols_, GridFF& gridFF ){
        prepareBuffers( nSystems_, nMols_, gridFF.grid.n, (float*)gridFF.FFPauli_f, (float*)gridFF.FFLondon_f, (float*)gridFF.FFelec_f );
    }
    
    void upload_mol2atoms(){ cl->upload  ( id_mol2atoms,  mol2atoms ); }
    void upload_poses    (){ cl->upload  ( id_poses,      poses     ); }
    void download_poses  (){ cl->download( id_poses,      poses     ); }
    void download_fposes (){ cl->download( id_fposes,     fposes    ); }
    
    void upload_poss  ( ){ cl->upload  ( id_poss, poss ); }
    void download_FEs ( ){ cl->download( id_FEs,  FEs  ); }
    
    void setupKernel_getFEtot( GridShape& grid, Vec3d testREQ, float alpha_, int nPoss ){        
        //__kernel void getFEtot(
        //    __read_only image3d_t  imgPauli,
        //    __read_only image3d_t  imgLondon,
        //    __read_only image3d_t  imgElec,
        //    __global  float4*  poss,
        //    __global  float4*  FEs,
        //    float4 dinvA,
        //    float4 dinvB,
        //    float4 dinvC,
        //    float4 PLQ
        
        int nLoc = _max( nMols, nAtomInMolMax );
        
        task_getFEtot->dim       = 1;
        task_getFEtot->local [0] = 1;
        task_getFEtot->global[0] = nPoss; 
        //pos0.setXYZ( (Vec3f)grid.pos0    );
        dA.f = (Vec3f)grid.diCell.a*(1.0/grid.n.a);
        dB.f = (Vec3f)grid.diCell.b*(1.0/grid.n.b);
        dC.f = (Vec3f)grid.diCell.c*(1.0/grid.n.c);
        testPLQ.f = (Vec3f)REQ2PLQ( testREQ, alpha );
        
        alpha = alpha_;
        task_getFEtot->args = {             
            BUFFarg(id_FFPauli), 
            BUFFarg(id_FFLondon), 
            BUFFarg(id_FFelec),    
            
            BUFFarg(id_poss),
            BUFFarg(id_FEs),
            REFarg(dA),
            REFarg(dB),
            REFarg(dC),
            REFarg(testPLQ)
        };
    }
    
    
    void setupKernel_getForceRigidSystemSurfGrid( GridShape& grid, float alpha_ ){        
        //__kernel void getForceRigidSystemSurfGrid(
        //    __read_only image3d_t  imgPauli,
        //    __read_only image3d_t  imgLondon,
        //    __read_only image3d_t  imgElec,
        //    // found - previously found molecular configurations
        //    __global  int2*    id_mol2atoms, // molTypes[type.x:type.y]
        //    //__global  int2*  confs,        // pointer to poses ... since we have constant number of molecules, we dont need this
        //    __global  float8*  atomsInTypes, // atoms in molecule types
        //    __global  float8*  poses,        // pos, qrot
        //    float4 dinvA,
        //    float4 dinvB,
        //    float4 dinvC,
        //    int nSystems,
        //    int nMols, // nMols should be approx local size
        //    float alpha
        
        int nLoc = _max( nMols, nAtomInMolMax );
        
        task_getForceRigidSystemSurfGrid->dim       = 1;
        task_getForceRigidSystemSurfGrid->local [0] = nLoc;
        task_getForceRigidSystemSurfGrid->global[0] = nSystems* task_getForceRigidSystemSurfGrid->local[0]; 
        //pos0.setXYZ( (Vec3f)grid.pos0    );
        dA  .setXYZ( (Vec3f)grid.diCell.a*(1.0/grid.n.a) );
        dB  .setXYZ( (Vec3f)grid.diCell.b*(1.0/grid.n.b) );
        dC  .setXYZ( (Vec3f)grid.diCell.c*(1.0/grid.n.c) );
        alpha = alpha_;
        task_getForceRigidSystemSurfGrid->args = { 
            //LBUFFarg(img_FFPauli), 
            //LBUFFarg(img_FFLondon), 
            //LBUFFarg(img_FFelec),
            
            BUFFarg(id_FFPauli), 
            BUFFarg(id_FFLondon), 
            BUFFarg(id_FFelec),    
            
            BUFFarg(id_mol2atoms),
            BUFFarg(id_atomsInTypes),
            BUFFarg(id_poses),
            BUFFarg(id_fposes),
            REFarg(dA),
            REFarg(dB),
            REFarg(dC),
            INTarg(nSystems),
            INTarg(nMols),
            FLOATarg(alpha)
        };
    }
    
    void evalForceGPU(){
       upload_poses   ();
       task_getForceRigidSystemSurfGrid->enque();
       download_fposes();  
       clFinish(cl->commands);      
    };
    
    int system2atoms( int isystem, float8* atoms ){
        int isoff   = isystem * nMols;
        int atom_count = 0;
        float8* atom0s = atomsInTypes.data();
        for(int imol=0; imol<nMols; imol++){
            //float* posi     = (float*)(poses+isoff+imol);
            Quat4f* posi     = (Quat4f*)(poses+isoff+imol);
            
            const int2& m2a = mol2atoms[imol];
            //printf( "isystem %i imol %i m2a (%i,%i) atom_count %i %i \n", isystem, imol, m2a.x, m2a.y, atom_count, atom_count );
            //frag2atoms( *((Vec3f*)(posi)), *((Quat4f*)(posi+4)), m2a.y, atom0s+m2a.x, atoms+atom_count );
            frag2atoms( posi[0].f, posi[1], m2a.y, atom0s+m2a.x, atoms+atom_count );
            atom_count += m2a.y;
        }
        return atom_count;
    }
    
    int evalForceCPU( int isystem, const GridFF& gridFF, float8* atoms, Vec3f* fatoms=0 ){
        printf("evalForceCPU \n");
        int isoff   = isystem * nMols;
        //float8* atom0s = atomsInTypes.data();
        
        Quat4f*  fs    = (Quat4f*)(fposes+isoff);
        Quat4f*  posei = (Quat4f*)( poses+isoff);
        
        float8*  atomim  = atoms;
        Vec3f*   fatomim = fatoms;
        
        const double R2ELEC = 1.0;
        
        for( int imol=0; imol<nMols; imol++ ){
            Vec3f force = Vec3fZero;
            Vec3f torq  = Vec3fZero;
            int natomi = mol2atoms[isoff+imol].y;
            
            Quat4f*  atomi = ((Quat4f*)atomim);
            Vec3f*  fatomi = fatomim;
            for(int ia=0; ia<natomi; ia++){ // atoms of molecule i
                
                const Vec3f& aposi  = atomi[0].f;
                const Vec3f& REQi   = atomi[1].f;
        
                // molecule grid interactions
                //float eps    = sqrt(REQi.y); //  THIS SHOULD BE ALREADY DONE
                
                //float eps   = sqrt(REQi.y);
                //float expar = exp(-alpha*REQi.x);
                //float CP =    eps*expar*expar;
                //float CL = -2*eps*expar;
                
                float expar    = exp(-alpha*REQi.x);
                float CP  =    REQi.y*expar*expar;
                float CL  = -2*REQi.y*expar;
                
                //printf( "CPU fgrid: imol %i ia %i ", imol, ia );
                
                Vec3d fd = Vec3dZero;
                gridFF.addForce( (Vec3d)aposi, (Vec3d){CP,CL,REQi.z}, fd );
                Vec3f f = (Vec3f)fd;
        
                //printf( "CPU fgrid: imol %i ia %i p(%5.5e,%5.5e,%5.5e) PLQ(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", imol, ia, aposi.x, aposi.y, aposi.z, CP,CL,REQi.z, f.x, f.y, f.z );
                //f = (Vec3f){0.1,0.0,0.0};
        
                if(fatoms){
                    fatomi[ia] = f ;
                    //printf( "FMG: imol %i ia %i f(%g,%g,%g) \n", imol, ia, f.x,f.y,f.z );
                    //printf( "FMG: imol %i ia %i f(%g,%g,%g) fatomi(%g,%g,%g) \n", imol, ia, f.x,f.y,f.z, fatomi[ia].x, fatomi[ia].y, fatomi[ia].z );
                }
                        
                force.add(f);
                torq.add_cross(aposi-posei->f, f);
                //printf( "%i tq.f %g tq.r %g r.r %g tq.tq %g \n", ia, torq.dot(f), torq.dot(aposi-posei->f), (aposi-posei->f).norm(), torq.norm() );
                
                atomi+=2;
            } // ia
            
            
            // ==== Molecule - Molecule Interaction
            float8* atomjm = atoms;
            for(int jmol=0; jmol<nMols; jmol++){
                int   natomj = mol2atoms[isoff+jmol].y;
                if (jmol!=imol){  // prevent self-interaction and reading outside buffer
                    
                    //float* atomj = ((float*)(atoms));
                    Quat4f* atomi  = ((Quat4f*)atomim);
                    Vec3f*  fatomi = fatomim;
                    for(int ia=0; ia<natomi; ia++){ // atoms of molecule i
                        //Vec3f aposi  = *(Vec3f*)(atomi  );
                        //Vec3f REQi   = *(Vec3f*)(atomi+4);
                        Vec3f aposi  = atomi[0].f;
                        Vec3f REQi   = atomi[1].f;
                        // molecule-molecule interaction
                        
                        Vec3f f = Vec3fZero;
                        
                        Quat4f* atomj = ((Quat4f*)atomjm);
                        for(int ja=0; ja<natomj; ja++){            // atoms of molecule j
                            //const Vec3f  dp   = *(Vec3f*)(atomj  ) - aposi;
                            //const Vec3f& REQj = *(Vec3f*)(atomj+4);
                            
                            const Vec3f  dp   = atomj[0].f - aposi;
                            const Vec3f& REQj = atomj[1].f;
                            
                            // force
                            float r0    = REQj.x + REQi.x;
                            float eps   = REQj.y * REQi.y; 
                            float cElec = REQj.z * REQi.z * COULOMB_CONST;
                            float r     = sqrt( dp.dot(dp) + +R2SAFE );
                            float expar = exp( alpha*(r-r0));
                            float fr    = eps*2*alpha*( expar*expar - expar ) - cElec/( r*r + R2ELEC );
                            
                            //double fr    = eps*2*alpha*( expar*expar - expar ) + COULOMB_CONST*qq/( r*r + R2ELEC );       
                            //printf( "CPU (%i,%i) (%i,%i) r %g expar %g fr %g kqq %g a %g eps %g \n", imol, ia, jmol, ja, r, expar, fr, cElec, alpha, eps );
                            f.add_mul( dp,  fr/r );
                            //e   += eps*( expar*expar - 2*expar ) + cElec*ir; // Energy
                            //atomj+=8;
                            atomj+=2;
                        } // ja

                        if(fatoms){
                            fatomi[ia].add( f );
                            //printf( "FMM: imol %i ia %i f(%g,%g,%g) \n", imol, ia, f.x,f.y,f.z );
                            //printf( "FMG: imol %i ia %i f(%g,%g,%g) fatomi(%g,%g,%g) \n", imol, ia, f.x,f.y,f.z, fatomi[ia].x, fatomi[ia].y, fatomi[ia].z );
                        }
                                    
                        force.add(f);
                        torq .add_cross(aposi-posei->f, f);
                        //torq .add_cross(f, aposi-posei->f);
                        
                        //atomi+=8;
                        atomi+=2;
                        
                    } // ia
                
                }
                
                atomjm += natomj; 
                
            } // jmol
            
            //float*  fs     = (float*)(fposes+isoff+imol);
            //*((Vec3f*)(fs  )) = force;
            //*((Vec3f*)(fs+4)) = torq;
            
            printf( "CPU imol %i f(%g,%g,%g) tq(%g,%g,%g) \n", imol, force.x, force.y, force.z, torq.x, torq.y, torq.z );
            
            fs[0].f = force;
            fs[1].f = torq;
            //printf( "imol %i force(%g,%g,%g)    torq (%g,%g,%g,%g) \n",  imol, force.x, force.y, force.z,    torq.x, torq.y, torq.z );
            //printf( "imol %i fs[0](%g,%g,%g,%g) fs[1](%g,%g,%g,%g) \n", imol, fs[0].x, fs[0].y, fs[0].z, fs[0].w,  fs[1].x, fs[1].y, fs[1].z, fs[1].w );
            fs+=2;

            //printf( "imol %i (%g,%g,%g,%g|%g,%g,%g,%g) \n", imol, fs[0], fs[1], fs[2], fs[3], fs[4], fs[5], fs[6], fs[7] );
            
            atomim  += natomi;
            fatomim += natomi;
            posei  +=2;
        } // imol
    } // evalForceCPU
    
    void moveSystemGD( int isystem, float dt, float sct, float scr ){
        int isoff   = isystem * nMols;
        Quat4f* ps = (Quat4f*)( poses + isoff); 
        Quat4f* fs = (Quat4f*)(fposes + isoff); 
        for(int imol=0; imol<nMols; imol++){
            ps[0].f.add_mul   ( fs[0].f, dt*sct  );
            ps[1]  .dRot_exact( dt*scr , fs[1].f );
            ps+=2; fs+=2;
        }
    }
       
};


// =======================
//      GridFF_OCL
// =======================

class GridFF_OCL{ public:

    OCLsystem* cl = 0;
    OCLtask *task_FFPLE = 0;
    
    Quat4i nGrid=Quat4iZero,npbc=Quat4iZero;
    Quat4f pos0=Quat4fZero,dA=Quat4fZero,dB=Quat4fZero,dC=Quat4fZero;
    float alpha = 0.0;

    int nAtoms    = 0;
    int nGridTot  = 0;
    
    int id_FFPauli    = -1;
    int id_FFLondon   = -1;
    int id_FFelec     = -1;
    int id_gridPoints = -1;
    int id_atoms      = -1;
    
    // ==== Functions

    void init( OCLsystem* cl_, char* fname ){
        cl = cl_; 
        int err = cl->buildProgram( fname );                              
        OCL_checkError(err, "cl.buildProgram");
        
        int id_evalPLE = cl->newKernel("evalPLE"); DEBUG
        task_FFPLE = new OCLtask( cl, id_evalPLE, 1, -1, 32 ); DEBUG;
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
        
        DEBUG;
        
        //task_FFPLE->global[0] = nGrid; 
        //task_FFPLE->args = { INTarg(nAtoms), BUFFarg(id_atoms), BUFFarg(id_gridPoints), BUFFarg(id_FFPauli), BUFFarg(id_FFLondon), BUFFarg(id_FFelec) };
        //task_FFPLE->print_arg_list();
                
        DEBUG;
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
    
    void downloadFF( GridFF& gridFF ){
        int n = gridFF.grid.getNtot();
        if( nGridTot!=n){ printf("ERROR: GridFF_OCL::downloadFF() Wrong number of grid points: n(%i) != nGrid(%i) \n, ", n, nGridTot ); exit(0); }
        
        gridFF.FFPauli_f  = new Quat4f[nGridTot];
        gridFF.FFLondon_f = new Quat4f[nGridTot];
        gridFF.FFelec_f   = new Quat4f[nGridTot];
        
        bool copyToDouble = true;
        
        if(gridFF.FFPauli_f ){  cl->download( id_FFPauli,  (float*)gridFF.FFPauli_f  ); if(copyToDouble ) float4ToVec3d( n, (float*)gridFF.FFPauli_f,  gridFF.FFPauli  ); printf("FFPauli  downloaded\n"); }
        if(gridFF.FFLondon_f){  cl->download( id_FFLondon, (float*)gridFF.FFLondon_f ); if(copyToDouble ) float4ToVec3d( n, (float*)gridFF.FFLondon_f, gridFF.FFLondon ); printf("FFLondon downloaded\n"); }
        if(gridFF.FFelec_f  ){  cl->download( id_FFelec ,  (float*)gridFF.FFelec_f   ); if(copyToDouble ) float4ToVec3d( n, (float*)gridFF.FFelec_f ,  gridFF.FFelec   ); printf("FFelec   downloaded\n"); }
    }
    
    void evalGridFFs( GridFF& gridFF, const Vec3i& nPBC ){
        printf( "gridFF.natoms %i \n", gridFF.natoms );
        prepareBuffers( gridFF.natoms, gridFF.grid.getNtot() ); DEBUG
        setupKernel( gridFF );
        uploadAtoms( gridFF.natoms, gridFF.apos, gridFF.aREQs ); DEBUG
        task_FFPLE->enque(); DEBUG
        //downloadFF( gridFF.grid.getNtot(), gridFF.FFPauli, gridFF.FFLondon, gridFF.FFelec ); DEBUG;
        downloadFF( gridFF );
        clFinish(cl->commands); DEBUG;
    }
        
};


#endif
