
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "fastmath.h"
#include "Shock1D.h"

ShockSystem1D   simulation;

int nMat;
ShockMaterial * materials; 

// ============ Exported functions

extern "C"{

    void init( char * fname ){  
        
    };

    void update( int n, double dt ){
        for(int i=0; i<n; i++){ simulation.update( dt ); };
    };
    
    //double* getPointer( int i ){    };

}
