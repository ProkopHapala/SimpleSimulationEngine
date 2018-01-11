#ifndef PBCsystem_h
#define PBCsystem_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "grids3D.h"

inline void triclinicBounds( const Mat3d& lvec, const Vec3d& pmin, const Vec3d& pmax, Vec3d& cmin, Vec3d& cmax ){
    //Vec3d p;, c000,c001,c010,c011,c000,c001,c010,c011;
    Vec3d cs[8];
    lvec.dot_to( {pmin.x,pmin.y,pmin.z}, cs[0] );
    lvec.dot_to( {pmin.x,pmin.y,pmax.z}, cs[1] );
    lvec.dot_to( {pmin.x,pmax.y,pmin.z}, cs[2] );
    lvec.dot_to( {pmin.x,pmax.y,pmax.z}, cs[3] );
    lvec.dot_to( {pmax.x,pmin.y,pmin.z}, cs[4] );
    lvec.dot_to( {pmax.x,pmin.y,pmax.z}, cs[5] );
    lvec.dot_to( {pmax.x,pmax.y,pmin.z}, cs[6] );
    lvec.dot_to( {pmax.x,pmax.y,pmax.z}, cs[7] );
    cmin={ 1e+300, 1e+300, 1e+300};
    cmax={-1e+300,-1e+300,-1e+300};
    for(int i=0; i<8; i++){
        cmin.setIfLower  (cs[i]);
        cmax.setIfGreater(cs[i]);
    }
}


class PBCsystem{ public:
    Mat3d  lvec;
    Mat3d  ilvec;
    Vec3i  npbc;

    CubeGridRuler ruler;

    int natoms;
    int   *type = NULL;
    Vec3d *pos  = NULL;

    int nAtomPBCmax=0;
    int nAtomPBC   =0;
    int * cellNs      = NULL;  //  index of
    int * cell2atom   = NULL;  //  index of  first atom in cell
    int   * atom2cell = NULL;
    int   * pbc_iatom = NULL;
    Vec3d * pbc_pos   = NULL;


    /*
    void allocateCells( ){
        cellNs    = new int [ruler.ntot  ];
        cell2atom = new int [ruler.ntot  ];
        //atoms     = new Atom[natoms*npbc ];

        nAtomPBCmax = natoms*npbc.a*npbc.b*npbc.c;
        atom2cell   = new int   [ nAtomPBCmax ];
        pbc_iatom   = new int   [ nAtomPBCmax ];
        pbc_pos     = new Vec3d [ nAtomPBCmax ];
    }
    */

    inline void setCell( const Mat3d& lvec_ ){ lvec=lvec_; lvec.invert_to( ilvec ); };

    void initRuler( double Rcut ){
        // find how large cubic grid we need to facilitate all interactions
        Vec3d pmin = {0.0,0.0,0.0};
        Vec3d pmax = lvec.a + lvec.b + lvec.c; // assuming all atoms are within the cell
        /*
        pmin={ 1e+300, 1e+300, 1e+300};
        pmax={-1e+300,-1e+300,-1e+300};
        for(int i=0; i<natoms; i++){ // is this worth it ? we can simply use corners of the pbc cell
            cmin.setIfLower  (pos[i]);
            cmax.setIfGreater(pos[i]);
        }
        */
        pmin.add( {Rcut,Rcut,Rcut} );
        pmax.sub( {Rcut,Rcut,Rcut} );
        ruler.setup( pmin, pmax, Rcut );
    }

    void atomsToCells(){
        Vec3d poff = (Vec3d){ruler.step*0.5,ruler.step*0.5,ruler.step*0.5};
        printf("DEBUG 1 \n");
        for(int i=0; i<ruler.ntot; i++){ cellNs[i] = 0; }
        printf("DEBUG 2 \n");

        Vec3d pmax = ruler.pos0 + ruler.span;
        Vec3d pcell;
        nAtomPBC=0;


        nAtomPBCmax = natoms*npbc.a*npbc.b*npbc.c;
        int*   atom2cell  = new int   [ nAtomPBCmax ];
        int*   pbc_iatom_ = new int   [ nAtomPBCmax ];
        Vec3d* pbc_pos_   = new Vec3d [ nAtomPBCmax ];

        for(int ia=0; ia<ruler.n.a; ia++){
            for(int ib=0; ib<ruler.n.b; ib++){
                for(int ic=0; ic<ruler.n.c; ic++){
                    pcell = lvec.a*ia + lvec.b*ib + lvec.c*ic;
                    for(int i=0; i<natoms; i++){
                        Vec3d p = pos[i] + pcell;
                        //if( (p.isGreater(ruler.pos0)&&(ruler.span.isGreater(p+ruler.pos0)) )  ){
                        if( p.isBetween( ruler.pos0, pmax ) ){ // or we can directly check distance from inner bbox ?
                            Vec3d dpos;
                            Vec3i ipos;
                            //ruler.pos2box( p, ipos, dpos );
                            //int icell = ruler.ixyz2i(ipos);
                            int icell = ruler.icell(p);

                            cellNs     [icell    ]++;
                            atom2cell  [nAtomPBC] = icell;
                            pbc_iatom_ [nAtomPBC] = i;
                            pbc_pos_   [nAtomPBC] = p;
                            nAtomPBC++;
                        }
                    }
                }
            }
        }

        pbc_iatom  = new int   [ nAtomPBC ];
        pbc_pos    = new Vec3d [ nAtomPBC ];

        printf("DEBUG 3 \n");
        int nsum = 0;
        for(int i=0; i<ruler.ntot; i++){
            cell2atom[i]  = nsum;
            nsum         += cellNs[i];
            cellNs[i]     = 0;
        }
        cell2atom[ruler.ntot] = nAtomPBC;
        printf("DEBUG 4 \n");
        //for(int i=0; i<ncell; i++){ printf( "cell2atom[%i] = %i \n", i, cell2atom[i] ); }
        for(int i=0; i<nAtomPBC; i++){
            int icell     = atom2cell[i];
            int ni        = cellNs   [icell];
            int i_        = cell2atom[icell] + ni;
            cellNs[icell] = ni+1;
            //atoms[i_]     = (Atom){ps[i],Qs[i]};
            pbc_iatom[i_] = pbc_iatom_[i];
            pbc_pos  [i_] = pbc_pos_  [i];
        }
        printf("DEBUG 5 \n");
        delete [] pbc_iatom_;
        delete [] pbc_pos_;
        //for(int i=0; i<n; i++){ ps[i]=pos_[i]; }
    }


}; // RigidSubstrate


#endif
