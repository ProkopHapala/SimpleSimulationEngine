

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

//#include <list>
#include <vector>

//#include <SDL2/SDL.h>
//#include <SDL2/SDL_opengl.h>
//#include "Draw.h"
//#include "Draw3D.h"
//#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"

/*
#include "Multipoles.h"
#include "PotentialFlow.h"
#include "grids3D.h"
#include "MultipoleGrid.h"
*/

//#include "AppSDL2OGL_3D.h"
//#include "testUtils.h"
//#include "SDL_utils.h"
//#include "Plot2D.h"

//#include "MMFF.h"

//#include "RARFF.h"
//#include "RARFF2.h"
#include "RARFFarr.h"

#define R2SAFE  1.0e-8f

// ============ Global Variables

//RARFF2     ff;
RARFF2arr  ff;
//std::list<RigidAtomType>  atomTypes;
std::vector<RigidAtomType*> atomTypes;

extern "C"{

// ========= Grid initialization

int insertAtomType( int nbond, int ihyb, double rbond0, double aMorse, double bMorse, double c6, double R2vdW ){
    RigidAtomType* atyp = new RigidAtomType();
    atomTypes.push_back( atyp );
    //RigidAtomType& atyp = atomTypes.back();
    atyp->nbond  =  nbond;    // number bonds
    atyp->rbond0 =  rbond0;
    atyp->aMorse =  aMorse;
    atyp->bMorse =  bMorse;
    atyp->c6     =  c6;
    atyp->R2vdW  =  R2vdW;
    printf("insertAtomType %i %i  %g %g %g %g %g ", nbond, ihyb, rbond0, aMorse, bMorse, c6, R2vdW );
    switch(ihyb){
        case 0: atyp->bh0s = (Vec3d*)sp1_hs; printf("sp1\n"); break;
        case 1: atyp->bh0s = (Vec3d*)sp2_hs; printf("sp2\n"); break;
        case 2: atyp->bh0s = (Vec3d*)sp3_hs; printf("sp3\n"); break;
    };
    return atomTypes.size()-1;
}

void ralloc(int natom){
    ff.realloc(natom);
    ff.cleanAux();
}

int*    getTypes (){ return (int*)   ff.types;  }
double* getPoss  (){ return (double*)ff.poss;   }
double* getQrots (){ return (double*)ff.qrots;  }
double* getHbonds(){ return (double*)ff.hbonds; }
double* getEbonds(){ return (double*)ff.ebonds; }

void setTypes( int natoms, int* types ){
    for(int i=0; i<natoms; i++){ ff.types[i]=atomTypes[types[i]]; };
}

/*
void setAtoms( int natoms, int* types, double* apos, double* qrots[i]){
    ff.realloc(natoms);
    for(int i=0; i<natoms; i++){
        if(randf()>0.5){ ff.atoms[i].type=&type1;  }else{ ff.atoms[i].type=&type2; }
        ff.atoms[i].pos  = apos [i];
        ff.atoms[i].qrot = qrots[i];
        ff.atoms[i].cleanAux();
    }
}
*/

double relaxNsteps( int nsteps, double F2conf, double dt, double damp ){
    double F2=1.0;
    for(int itr=0; itr<nsteps; itr++){
        ff.cleanAtomForce();
        ff.projectBonds();
        ff.interEF();
        ff.evalTorques();
        F2 = ff.evalF2pos() + ff.evalF2rot(); // some scaling ?
        //printf( "itr %i F2 %g \n", itr, F2 ); 
        if(F2<F2conf) break;
        ff.moveMDdamp(dt, damp);
    }
    return F2;
}

} // extern "C"{


