
#ifndef AtomicConfiguration_h
#define AtomicConfiguration_h

#include "Vec3.h"

#include "PBCsystem.h"

class FastAtomicMetric{

    int   natoms;
    int   * types=0;
    Vec3d * pos  =0;

    double Rcut;

    CubeGridRuler ruler;
    int  nprj=0;
    int  *  cellNs         = 0;  //  index of
    int  ** cell2atoms     = 0;
    int  *  cell2atomBuff  = 0;  //  index of  first atom in cell
    //int   * atom2cell = 0;

    //void makeCell2atom(  ){}

    void initRuler( Vec3d pmin, Vec3d pmax, double step ){
        Rcut=Rcut;
        ruler.setup( pmin, pmax, step );
        _realloc(cellNs    ,ruler.ntot);
        _realloc(cell2atoms,ruler.ntot);
    };

    void projectAtoms( int natoms_, int* types_, Vec3d* pos_, double Rcut ){
        //
        //  NOTE : all symmetry folding should be here ( periodicity, inversion.. )
        //
        natoms=natoms_; types=types_; pos=pos_;
        int* atom2cells = new int[natoms*8];
        int* atomNs     = new int[natoms  ];
        int* atom2cells_=atom2cells;
        nprj=0;
        for(int i=0; i<natoms; i++){
            int nc = ruler.overpSphere( pos[i], Rcut, atom2cells_ );
            atom2cells_+=nc;
            nprj       +=nc;
            atomNs[i]   =nprj;
        }
        _realloc( cell2atomBuff, nprj);
        for(int i=0; i<ruler.ntot; i++){ cellNs[i]=0;             };
        for(int i=0; i<nprj;       i++){ cellNs[atom2cells[i]]++; };
        int i0=0;
        for(int i=0; i<ruler.ntot; i++){ cell2atoms[i]=cell2atomBuff+i0; i0+=cellNs[i]; };
        for(int i=0; i<ruler.ntot; i++){ cellNs[i]=0; };
        int iatom=0;
        for(int i=0; i<nprj;       i++){
            int ic=atom2cells[i];
            int j=cellNs[ic];
            if(i>atomNs[iatom])iatom++;
            cell2atoms[ic][j] = iatom;
            cellNs[ic]++;
        };
        delete [] atom2cells;
        delete [] atomNs;
    };

    double dist( int * types_, Vec3d * pos_ ){
        double R2    = 0;
        double R2cut = sq(Rcut);
        for(int i=0; i<natoms; i++){
            int   ti = types_[i];
            Vec3d pi = pos_  [i];
            double r2min = 0;
            int ic = ruler.icell( pi );
            for(int jj=0; jj<cellNs[ic]; jj++){
                int j = cell2atoms[ic][jj];
                if(ti==types[j]){
                    double r2 = (pos[j]-pi).norm2();
                    if( (r2<R2cut)&&(r2<r2min) )r2min=r2;
                }
            }
            R2+=r2min;
        }
        return R2;
    };

};

class AtomicConfiguration{
    int n;
    int   * types=0;
    Vec3d * pos  =0;

    double dist_On2(const AtomicConfiguration& p){
        double R2 = 0;
        for(int i=0; i<p.n; i++){
            int   ti=p.types[i];
            Vec3d pi=p.pos  [i];
            double r2min = 0;
            for(int j=0; j<n; j++){
                if(ti==types[j]){
                    double r2 = (pos[j] - pi).norm2();
                    if(r2<r2min) r2min=r2;
                }
            }
            R2 += r2min;
        }
        return R2;
    }

};

#endif
