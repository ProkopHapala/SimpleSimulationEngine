
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "Multipoles.h"
#include "DynamicOpt.h"

#include "radial_splines.h"
#include "AtomTypes.h"
#include "MoleculeType.h"
#include "MolecularWorld.h"

//#include "testUtils.h"

// ============ Global Variables

MolecularWorld world;

bool   converged   = false;
//double fmaxConverg = 0.00001;

// ============ Exported functions

extern "C"{

    void initWorld( char * workdir ){
        char fname[256];
        world.fromDir( workdir, "atomTypes.ini", "molTypes.ini", "instances.ini" );
        //world.fromDir( "inputs/", "atomTypes.ini", "molTypes.ini", "instances.ini" );
        //strcpy(fname,workdir ); strcat(fname,"splines.ini"); world.loadSplines(fname); //exit(0);
        strcpy(fname,workdir ); strcat(fname,"bonds.ini");   world.loadBonds  (fname);

        world.checkBonds( 0.9, 1.2 );
        world.setCutoff ( 6.0 );
        world.makeFF    ( );
        world.optimizer->initOpt( 0.05, 0.15 );
    }

    double relax( int niter, double fmaxConverg ){
        int iter;
        for(iter=0; iter<niter; iter++){
            world.rigidOptStep( );
            printf(" opt step %i fmax %g \n", world.optimizer->stepsDone, world.fmax );
            if( world.fmax < fmaxConverg ){
                converged = true;
                return world.fmax;
            }
        }
        return world.fmax;
    }

    void exportAtoms( char * fname ){
        FILE * fout = fopen(fname, "w");
        char str[256];
        sprintf(str,"# fmax = %g", world.fmax );
        world.exportAtomsXYZ( fout, str );
        fclose(fout);
    }

    void testMultipole( int order,
        int np, double * ps_, double * Qs,
        int nsamples, double * psamples_, double * Eref, double * Eaprox
    ){
        Vec3d * ps       = (Vec3d*)ps_;
        Vec3d * psamples = (Vec3d*)psamples_;
        double coefs[10];
        getMultiPole( {0.0,0.0,0.0}, np, (Vec3d*)ps, Qs, order, coefs );
        for(int i=0;i<10;i++){printf("%f ", coefs[i]);} printf("\n");
        for(int i=0; i<nsamples; i++){
            Vec3d& p = psamples[i];
            Eref  [i] = evalEelectrostatic( p, 1.0, np, (Vec3d*)ps, Qs );
            Eaprox[i] = Emultipole        ( p, order, coefs );
        }
    }

}
