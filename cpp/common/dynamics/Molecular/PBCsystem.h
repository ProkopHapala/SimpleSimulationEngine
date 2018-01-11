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

class SymFuncFF_interface{
    virtual void acumFF( int nbond, Vec3d *hbond, double* rbond, int *btype )=0;
};

class PBCsystem{ public:
    Mat3d  lvec;
    Mat3d  ilvec;
    //Vec3d  pbc0;
    Vec3i  npbc,pbc0,pbc1;
    Vec3d cmin,cmax; // debug

    CubeGridRuler ruler;
    double Rcut,R2cut;
    double R2min = 0.01;

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

    void initRuler( double Rcut_ ){
        Rcut=Rcut_; R2cut=Rcut*Rcut;
        // find how large cubic grid we need to facilitate all interactions
        //Vec3d pmin = {0.0,0.0,0.0};
        //Vec3d pmax = lvec.a + lvec.b + lvec.c; // assuming all atoms are within the cell

        Vec3d pmin={ 1e+300, 1e+300, 1e+300};
        Vec3d pmax={-1e+300,-1e+300,-1e+300};
        for(int i=0; i<natoms; i++){ // is this worth it ? we can simply use corners of the pbc cell
            pmin.setIfLower  (pos[i]);
            pmax.setIfGreater(pos[i]);
        }
        //printf( "pmin %f %f %f\n",pmin.x,pmin.y,pmin.z);
        //printf( "pmax %f %f %f\n",pmax.x,pmax.y,pmax.z );

        pmin.sub( {Rcut,Rcut,Rcut} );
        pmax.add( {Rcut,Rcut,Rcut} );
        ruler.setup( pmin, pmax, Rcut );

        printf( "ruler.ntot %i \n", ruler.ntot );
        // alocate cells
        if(cellNs)    delete [] cellNs;    cellNs    = new int[ruler.ntot  ];
        if(cell2atom) delete [] cell2atom; cell2atom = new int[ruler.ntot  ];

        //printf( "pmin %f %f %f\n",pmin.x,pmin.y,pmin.z);
        //printf( "pmax %f %f %f\n",pmax.x,pmax.y,pmax.z );

        //Vec3d cmin,cmax;
        triclinicBounds( ilvec, pmin, pmax, cmin, cmax );
        pbc0.a = (int)(cmin.a-1); pbc1.a = (int)(cmax.a+1);
        pbc0.b = (int)(cmin.b-1); pbc1.b = (int)(cmax.b+1);
        pbc0.c = (int)(cmin.c-1); pbc1.c = (int)(cmax.c+1);
        npbc=pbc1-pbc0;
        //npbc.a = (int)fmax( -cmin.a, cmax.a )+1;
        //npbc.b = (int)fmax( -cmin.b, cmax.b )+1;
        //npbc.c = (int)fmax( -cmin.c, cmax.c )+1;

    }

    void atomsToCells(){
        // allocations
        nAtomPBCmax = natoms*npbc.a*npbc.b*npbc.c;
        printf( "nAtomPBCmax %i \n", nAtomPBCmax );
        if(atom2cell) delete [] atom2cell; atom2cell = new int[ nAtomPBCmax ];
        int*   pbc_iatom_ = new int   [ nAtomPBCmax ];
        Vec3d* pbc_pos_   = new Vec3d [ nAtomPBCmax ];
        for(int i=0; i<ruler.ntot; i++){ cellNs[i] = 0; }
        printf("DEBUG 2 \n");
        Vec3d pcell;
        nAtomPBC=0;
        for(int ia=pbc0.a; ia<pbc1.a; ia++){
            for(int ib=pbc0.b; ib<pbc1.b; ib++){
                for(int ic=pbc0.c; ic<pbc1.c; ic++){
                    pcell = lvec.a*ia + lvec.b*ib + lvec.c*ic;
                    for(int i=0; i<natoms; i++){
                        Vec3d p = pos[i] + pcell;
                        //if( (p.isGreater(ruler.pos0)&&(ruler.span.isGreater(p+ruler.pos0)) )  ){
                        if( p.isBetween( ruler.pos0, ruler.pmax ) ){ // or we can directly check distance from inner bbox ?
                            //Vec3d dpos;
                            //Vec3i ipos;
                            //ruler.pos2box( p, ipos, dpos );
                            //int icell = ruler.ixyz2i(ipos);
                            int icell = ruler.icell(p);
                            //printf( "%i (%i,%i,%i) %i \n", nAtomPBC, ia, ib, ic, i );
                            cellNs     [icell   ]++;
                            atom2cell  [nAtomPBC] = icell;
                            pbc_iatom_ [nAtomPBC] = i;
                            pbc_pos_   [nAtomPBC] = p;
                            nAtomPBC++;
                        }
                    }
                }
            }
        }
        //pbc_pos = pbc_pos_;
        //printf("finished \n");

        if(pbc_iatom) delete [] pbc_iatom; pbc_iatom  = new int   [ nAtomPBC ];
        if(pbc_pos  ) delete [] pbc_pos;   pbc_pos    = new Vec3d [ nAtomPBC ];

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

    void symfunc_R  ( int type0, int type1,            double r1                       ){};
    void symfunc_Ang( int type0, int type1, int type2, double r1, double r2, double ca ){}; // just placeholder

    void gatherSymFunctions( ){
        //int icells[27];
        Vec3d  *hbond     = new Vec3d [nAtomPBC];
        double *lbond     = new double[nAtomPBC];
        int    *neighType = new   int [nAtomPBC];

        for(int i=0; i<natoms; i++ ){
            Vec3d p  = pos[i];
            Vec3i ip = ruler.ipcell(p);
            int ii=0;
            int nbonds=0;
            int typei = type[i];
            for( int ix=-1; ix<2; ix++ ){ for( int iy=-1; iy<2; iy++ ){ for( int iz=-1; iz<2; iz++ ){
                int icell = ruler.ixyz2i({ip.x+ix,ip.y+iy,ip.z+iz});
                int j0  = cell2atom[icell];
                int j1  = j0+cellNs[icell];
                for( int j=j0; j<j1; j++ ){
                    Vec3d d = pbc_pos[j] - p;
                    double r2 = d.norm2();
                    if( (r2 < R2cut)&&(r2>R2min) ){
                        double r = sqrt(r2);
                        lbond    [nbonds] = r;
                        hbond    [nbonds].set_mul(d,1/r);
                        neighType[nbonds] = type[pbc_iatom[j]];
                        nbonds++;
                    }
                }
            }}}
            // TODO: this should be made idependent (class call?,function pointer?, lambda function? )
            for(int ib=0; ib<nbonds; ib++){
                Vec3d  hi    = hbond[ib];
                double ri    = lbond[ib];
                int    type1 = neighType[ib];
                symfunc_R( typei, type1, ri );
                for(int jb=0; jb<nbonds; jb++){
                    if(ib==jb) continue;
                    double ca = hi.dot( hbond[jb] );
                    symfunc_Ang( typei, type1, neighType[jb], ri, lbond[jb], ca );
                }
            }
        }
        delete [] hbond;
        delete [] lbond;
        delete [] neighType;
    }


}; // RigidSubstrate


#endif
