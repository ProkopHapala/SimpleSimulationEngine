
#ifndef AtomicConfiguration_h
#define AtomicConfiguration_h

#include "Vec3.h"
#include "PBCsystem.h"
#include "globals.h"


class AtomicConfiguration{ public:
    int   natoms;
    int   * types=0;
    Vec3d * pos  =0;

    double dist(const AtomicConfiguration& p){
        double R2 = 0;
        for(int i=0; i<p.natoms; i++){
            int   ti=p.types[i];
            Vec3d pi=p.pos  [i];
            double r2min = 1e+300;
            for(int j=0; j<natoms; j++){
                //printf( "%i %i(%f,%f,%f) | %i %i(%f,%f,%f)\n", i, p.types[i], p.pos[i].x,p.pos[i].y,p.pos[i].z,   j, types[j], pos[j].x,pos[j].y,pos[j].z  );
                if(ti==types[j]){
                    double r2 = (pos[j] - pi).norm2();   //printf( "%i %i %i %f \n", i,  j,  p.types[i], r2 );
                    if(r2<r2min) r2min=r2;
                }
            }
            R2 += r2min;
            //printf( "=== %i %i %f \n", i, p.types[i], r2min );
        }
        return R2;
    }

    void realloc(int n){ natoms=n; _realloc(types,natoms); _realloc(pos,natoms); }

    void bind( int n, int* types_, Vec3d* pos_ ){ natoms=n; types=types_; pos=pos_; }

    void copyOf(const AtomicConfiguration& p){
        if(natoms!=p.natoms)realloc(p.natoms);
        //for(int i=0; i<natoms; i++){ pos[i]=p.pos[i]; types[i] = p.types[i]; }
        memcpy( types, p.types, sizeof(int)  *natoms );
        memcpy( pos,   p.pos,   sizeof(Vec3d)*natoms );
    }

    AtomicConfiguration(){};
    AtomicConfiguration(int n){};
    AtomicConfiguration(const AtomicConfiguration& p){ copyOf(p);  };

};

class FastAtomicMetric : public AtomicConfiguration { public:

    //int   natoms;
    //int   * types=0;
    //Vec3d * pos  =0;

    double Rcut=0;

    CubeGridRuler ruler;
    int  nprj=0;
    int  *  cellNs         = 0;  //  index of
    int  ** cell2atoms     = 0;
    int  *  cell2atomBuff  = 0;  //  index of  first atom in cell
    //int   * atom2cell = 0;

    // this may be just local in function
    int* atom2cells = 0;
    int* atomNs     = 0;

    //void makeCell2atom(  ){}

    void initRuler( Vec3d pmin, Vec3d pmax, double step ){ 
        ruler.setup( pmin, pmax, step );
        _realloc(cellNs    ,ruler.ntot);
        _realloc(cell2atoms,ruler.ntot);
        Rcut=ruler.step*0.5-1e-3;
    };

    //void project( int natoms_, int* types_, Vec3d* pos_, double Rcut ){
    void toCells(){
        //
        //  NOTE : all symmetry folding should be here ( periodicity, inversion.. )
        //
        //natoms=natoms_; types=types_; pos=pos_;
        //Rcut=Rcut_;
        //int* atom2cells = new int[natoms*8];
        //int* atomNs     = new int[natoms  ];
        DEBUG
        DBG( "toCells natoms, atom2cells, atomNs %i %i %i \n", natoms, atom2cells, atomNs );
        //printf( "toCells natoms, atom2cells, atomNs %i %i %i \n", natoms, atom2cells, atomNs );
        _realloc(atom2cells,natoms*8);
        _realloc(atomNs    ,natoms  );
        int* atom2cells_=atom2cells;
        nprj=0;
        DEBUG
        //printf("--- toCells 1 \n");
        for(int i=0; i<ruler.ntot; i++){ cellNs[i]=0;             };
        for(int i=0; i<natoms; i++){
            int nc = ruler.overlap_Sphere( pos[i], Rcut, atom2cells_ );
            nprj       +=nc;
            atomNs[i]   =nc;
            for(int j=0; j<nc; j++){
                cellNs[atom2cells_[j]]++;
                //printf( "%i ", atom2cells_[j] );
            };
            //printf( "|atom %i N %i %i \n", i, atomNs[i], nprj );
            atom2cells_+=nc;
        }
        DEBUG
        //printf( "nprj %i \n", nprj );
        _realloc( cell2atomBuff, nprj);
          // printf("toCells 2 \n");
        //for(int i=0; i<nprj;       i++){ cellNs[atom2cells[i]]++;
        //    //printf(" %i %i %i \n", i, atom2cells[i], cellNs[atom2cells[i]] );
        //};   // printf("toCells 3 \n");
        int i0=0;
        for(int i=0; i<ruler.ntot; i++){ cell2atoms[i]=cell2atomBuff+i0; i0+=cellNs[i]; cellNs[i]=0; // this can be much simplified if we allocate uper bound for cell2atomBuff
            //printf( "%i %i %i \n", i, cell2atoms[i]-cell2atomBuff, cellNs[i] );
        };
        DEBUG
        //printf("toCells 4 \n");
        //for(int i=0; i<ruler.ntot; i++){ cellNs[i]=0; };                                   printf("toCells 5 \n");
        int ia=-1,nprj_=0;
        for(int i=0; i<nprj; i++){
            if(i>=nprj_){ ia++; nprj_+=atomNs[ia]; };
            int ic=atom2cells[i];
            int j =cellNs[ic];
            cell2atoms[ic][j] = ia;
            cellNs[ic]++;
            //printf( " %i [%i][%i] ia %i N %i \n", i, ic, j, ia, cellNs[ic] );
        };
        DEBUG
        for(int i=0; i<nprj; i++){ };
        //printf("toCells 6 \n");
        //delete [] atom2cells;
        //delete [] atomNs;
        //exit(0);
    };

    inline void toCells(double Rcut_){ Rcut=Rcut_; toCells(); }

    double dist( int n, int * types_, Vec3d * pos_ )const{
        double R2    = 0;
        double R2cut = sq(Rcut);
        for(int i=0; i<n; i++){
            int   ti = types_[i];
            Vec3d pi = pos_  [i];
            double r2min = 1e+300;
            int ic = ruler.icell( pi );
            //printf( "== i %i ic %i cellN %i \n", i, ic, cellNs[ic] );
            for(int jj=0; jj<cellNs[ic]; jj++){
                int j = cell2atoms[ic][jj];
                //printf( "%i %i %i %i  \n", i, j, ic, types[j] );
                if(ti==types[j]){
                    double r2 = (pos[j]-pi).norm2();
                    //printf( "%i %i %f  \n", i, j, r2 );
                    if( (r2<R2cut)&&(r2<r2min) )r2min=r2;
                }
            }
            R2+=r2min;
        }
        return R2;
    };

    int findNeighs( const Vec3d& p, double Rcut, int* out )const{
        double R2cut = Rcut*Rcut;
        int nfound=0;
        int ic = ruler.icell(p);
        //if((ic<0)||(ic>ruler.ntot)) return -1;
        //printf( "CellN[%i] = %i \n", ic, cellNs[ic] );
        for(int jj=0; jj<cellNs[ic]; jj++){
            int j = cell2atoms[ic][jj];
            double r2 = (pos[j]-p).norm2();
            //printf( "%i %f %f \n", j, r2, R2cut );
            if(r2<R2cut){ out[nfound]=j; nfound++;};
        }
        return nfound;
    };

    inline double dist( const AtomicConfiguration& p )const{ return dist( p.natoms, p.types, p.pos ); }


};

#endif
