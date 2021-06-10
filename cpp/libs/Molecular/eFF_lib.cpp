
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <vector>
#include <unordered_map>
#include <string>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
//#include "VecN.h"

#include "testUtils.h"

//#include "integration.h"
//#include "AOIntegrals.h"
//#include "AOrotations.h"

Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;
int i_DEBUG=0;

#include "InteractionsGauss.h"
#include "eFF.h"

// ============ Global Variables

EFF ff;

std::unordered_map<std::string, double*>  buffers;
std::unordered_map<std::string, int*>     ibuffers;

extern "C"{
// ========= Grid initialization

double eval(){ return ff.eval(); };

void init_buffers(){
    // atoms (ions)
    buffers.insert( { "pDOFs",   ff.pDOFs } );
    buffers.insert( { "fDOFs",   ff.fDOFs } );
    buffers.insert( { "apos",   (double*)ff.apos   } );
    buffers.insert( { "aforce", (double*)ff.aforce } );
    buffers.insert( { "epos",   (double*)ff.epos   } );
    buffers.insert( { "eforce", (double*)ff.eforce } );
    buffers.insert( { "esize",  (double*)ff.esize   } );
    buffers.insert( { "fsize",  (double*)ff.fsize } );
    //buffers.insert ( { "aQ",            ff.aQ    } );
    //buffers.insert ( { "aAbWs",         (double*)ff.aAbWs } );
    //buffers.insert ( { "eAbWs",         (double*)ff.eAbWs } );
    buffers.insert ( { "aPars",         (double*)ff.aPars } );
    ibuffers.insert( { "espin",         ff.espin  } );
}

bool load_xyz( const char* fname ){ 
    //printf( "load_xyz \n" );
    bool b = ff.loadFromFile_xyz( fname );
    init_buffers();
    return b; 
}

void init( int na, int ne ){
    ff.realloc ( na, ne );
    //ff.autoAbWs( default_aAbWs, default_eAbWs );
    init_buffers();
}

void info(){ ff.info(); }

double* getEnergyPointer(){ return &ff.Ek; }
int*    getDimPointer   (){ return &ff.ne; }

double* getBuff(const char* name){ 
    auto got = buffers.find( name );
    if( got==buffers.end() ){ printf( "buffer[%s] NOT FOUNT !!! \n", name ); return 0;        }
    else                    { return got->second; }
}

void setBuff(const char* name, double* buff){ 
    buffers[name] = buff;
}

int* getIBuff(const char* name){ 
    auto got = ibuffers.find( name );
    if( got == ibuffers.end() ){ return 0;        }
    else                    { return got->second; }
}

void setIBuff(const char* name, int* buff){ 
    ibuffers[name] = buff;
    //auto got = buffers.find( name );
    //if( got==buffers.end() ){ return null;        }
    //else                    { return got->second; }
}

//#define NEWBUFF(name,N)   double* name = new double[N]; buffers.insert( {#name, name} );

void setPauliModel(int i  ){ ff.iPauliModel = i; }
void setKPauli( double KPauli ){
    if( ff.iPauliModel == 0 ){ ff.KPauliOverlap=KPauli; }
    else                  { ff.KPauliKin=KPauli;     }
}

} // extern "C"
