
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

    void init_materials( char * fname ){  
        
    };

    void update( int n, double dt ){
        for(int i=0; i<n; i++){ simulation.update( dt ); };
    };
    
    void init_python( int n, double * bounds, double * bforce, double * velocity, int * imats  ){
        simulation.nlayers   = n;
        simulation.bounds    = bounds;
        simulation.bforce    = bforce;
        //simulation.imass     = imass;
        simulation.velocity  = velocity;
        simulation.cells     = new ShockVolume[simulation.nlayers]; 
        for(int i=0; i<n; i++){ simulation.cells[i].material = materials + imats[i]; };
    }
        
    void getArray( int which, int n, double * buff ){
        switch( which ){
            case( 0 ): for(int i=0; i<n; i++ ){ buff[i] = simulation.cells[i].mass; }  break;
            case( 1 ): for(int i=0; i<n; i++ ){ buff[i] = simulation.cells[i].V; }  break;          
            case( 2 ): for(int i=0; i<n; i++ ){ buff[i] = simulation.cells[i].p; }  break;       
        }
    }

    void setArray( int which, int n, double * buff ){
        switch( which ){
            case( 0 ): for(int i=0; i<n; i++ ){ buff[i] = simulation.cells[i].mass; }  break;
            case( 1 ): for(int i=0; i<n; i++ ){ buff[i] = simulation.cells[i].V; }  break;          
            case( 2 ): for(int i=0; i<n; i++ ){ buff[i] = simulation.cells[i].p; }  break;       
        }
    }
    
}
