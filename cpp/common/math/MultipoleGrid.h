
#ifndef  MultipoleGrid_h
#define  MultipoleGrid_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec3.h"

#include "grids3D.h"
#include "Multipoles.h"

class Atom{  public:
    Vec3d  pos;
    double Q;
};

class MultipoleGrid{ public:
    CubeGridRuler ruler;

    int perCell = 4;
    int nAtomMax;

    Atom   * atoms;
    double * coefs;

    int * cellNs;     //  index of
    int * cell2atom;  //  index of  first atom in cell
    int * atom2cell;

    void allocate( int nAtomMax_, Vec3i nxyz ){
        nAtomMax = nAtomMax_; ruler.setn(nxyz);
        cellNs    = new int [ruler.ntot];
        cell2atom = new int [ruler.ntot];
        coefs     = new double[ruler.ntot*perCell];
        atoms     = new Atom[nAtomMax_];
        atom2cell = new int [nAtomMax_];

    }

    void projectPointCharges( int n, Vec3d * ps, double * Qs ){
        Vec3d poff = (Vec3d){ruler.step*0.5,ruler.step*0.5,ruler.step*0.5};
        for(int i=0; i<n; i++){
            Vec3d dpos,pos = ps[i];
            Vec3i ipos;
            ruler.pos2box( pos, ipos, dpos );
            int icell = ruler.ixyz2i(ipos);
            getMultipole ( dpos - poff, Qs[i], 1, coefs+(icell*4) );
        }
    };

    //void getForce( const Vec3d& pos, double Q, Vec3d& force ){}

    void atomsToCells( int n, Vec3d * ps, double * Qs ){
        Vec3d poff = (Vec3d){ruler.step*0.5,ruler.step*0.5,ruler.step*0.5};
        printf("DEBUG 1 \n");
        for(int i=0; i<ruler.ntot; i++){ cellNs[i] = 0; }
        printf("DEBUG 2 \n");
        for(int i=0; i<n; i++){
            Vec3d dpos,pos = ps[i];
            Vec3i ipos;
            ruler.pos2box( pos, ipos, dpos );
            int icell = ruler.ixyz2i(ipos);
            //printf("%i (%g,%g,%g) (%i,%i,%i) %i \n", i, pos.x, pos.y, pos.z, ipos.x, ipos.y, ipos.z, icell );
            getMultipole ( dpos - poff, Qs[i], 1, coefs+(icell*4) );
            cellNs   [icell]++;
            atom2cell[i] = icell;
        }
        printf("DEBUG 3 \n");
        int nsum = 0;
        for(int i=0; i<ruler.ntot; i++){
            cell2atom[i]  = nsum;
            nsum         += cellNs[i];
            cellNs[i]     = 0;
        }
        cell2atom[ruler.ntot] = n;
        printf("DEBUG 4 \n");
        //for(int i=0; i<ncell; i++){ printf( "cell2atom[%i] = %i \n", i, cell2atom[i] ); }
        for(int i=0; i<n; i++){
            int icell     = atom2cell[i];
            int ni        = cellNs   [icell];
            int i_        = cell2atom[icell] + ni;
            cellNs[icell] = ni+1;
            atoms[i_]     = (Atom){ps[i],Qs[i]};
        }
        printf("DEBUG 5 \n");
        //for(int i=0; i<n; i++){ ps[i]=pos_[i]; }
    }

};

#endif

