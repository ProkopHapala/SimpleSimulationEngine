
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

#include "PointCloudComparator.h"
#include "SphereTreeND.h"

//#include "testUtils.h"

// ============ Global Variables

MolecularWorld world;

bool   converged   = false;
//double fmaxConverg = 0.00001;


PointCloudComparator comp;
TypePointsComparator compT;

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

    double * getOptPointer_pos(){ for(int i=0; i<world.optimizer->n; i++){ printf("%i %f \n",i,world.optimizer->pos[i]);}; return world.optimizer->pos; }
    int      getOptSize       (){ return world.optimizer->n;   }

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

    int getNAtoms() { return world.getNAtoms(); };
    int getAtomPos  ( double * buff ){ return world.getAtomPos  ( (Vec3d*)buff ); };
    int getAtomTypes( int    * buff ){ return world.getAtomTypes(         buff ); };

    int atoms2map(){world.atoms2map();};

    double collisionForce( int imol, double * pose, double * fpose ){
        world.collisionForce(  imol,  *((Vec3d*)pose), *((Quat4d*)(pose+3)),  *((Vec3d*)fpose), *((Quat4d*)(fpose+3)) );
        double ftot = 0;
        for(int i=0; i<7; i++){  double f = fpose[i]; ftot+=f*f; }
        return ftot;
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

    void initComparator( int n, double * points ){
        comp.allocate( n );
        comp.setGrid ( 1.0, 0.5 ); // rmax is smaller than rmin that is not mistake: rmax maximim distance between corresponding atoms, rmin minimal distance between different atoms
        comp.setRefPoints( n, (Vec3d*)points );
    }

    double compDistance( double * points ){
        //return comp.dist_overlap_N2( (Vec3d*)points );
        //return comp.dist_maxmin_N2( (Vec3d*)points );
        return comp.dist_maxmin_N2( (Vec3d*)points );
    }

    void initComparatorT( int n, double * points, int * ptypes ){
        int ntyp = compT.setupTypes  ( n, ptypes );                  printf( "DEBUG 1 %i\n", ntyp );
        compT.setRefPoints( (Vec3d*)points, ptypes );     printf( "DEBUG 2\n" );
        compT.initRefGrid ( 1.0, 0.5 );                   printf( "DEBUG 3\n" );
    }

    double compDistanceT( int n, double * points, int * types ){
        return compT.dist( n, (Vec3d*)points, types );
    }

    void getPlaneWaveDescriptor( double * center_, int np, double * points, int nk,  double * ks, double * coefs ){
        Vec3d center; center.set( center_[0], center_[1], center_[2] );
        return getPlaneWaveDescriptor( center, np, (Vec3d*)points, nk,  (Vec3d*)ks, coefs );
    }

}
