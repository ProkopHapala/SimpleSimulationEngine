
#ifndef MMFF_h
#define MMFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"

#include "integerOps.h"
#include <unordered_map>
#include <vector>
#include <string>
// forcefield parameters

//inline BondType2id
//inline BondTypeDec

#define RSAFE   1.0e-4f
#define R2SAFE  1.0e-8f
#define F2MAX   10.0f

#define SIGN_MASK 2147483648

//inline int invIndex( int i ){ return i^SIGN_MASK; }

// compied from RigidMolecule.h      ============>   TODO : make common lib of inter-molecular forcefields formulas

// ================================
// ====   Forcefield Functions
// ================================

inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double q ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*q*-14.3996448915f )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceQ( const Vec3d& dp, Vec3d& f, double q ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double fr   = ( ir*q*-14.3996448915f )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceLJ( const Vec3d& dp, Vec3d& f, double r0, double eps ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceExp( const Vec3d& dp, Vec3d& f, double r0, double eps, double alpha ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double r    = dp.norm();
    double E    = eps*exp( alpha*(r-r0) );
    double fr   = alpha*E/(r+RSAFE);
    f.add_mul( dp, fr );
    //f.add_mul( dp, 1/(dp.norm2()+R2SAFE) ); // WARRNING DEBUG !!!!
}

inline Vec3d REQ2PLQ( Vec3d REQ, double alpha ){
    double eps   = sqrt(REQ.y);
    double expar = exp(-alpha*REQ.x);
    double CP =    eps*expar*expar;
    double CL = -2*eps*expar;
    //printf( "REQ2PLQ: %g %g %g  ->  %g %g\n", REQ.x, eps, alpha,   CP, CL );
    return (Vec3d){ CP, CL, REQ.z };
}

inline Vec3d getForceSpringPlane( const Vec3d& p, const Vec3d& normal, double c0, double k ){
    double cdot = normal.dot(p) - c0;
    return normal * (cdot * k);
}

inline Vec3d getForceHamakerPlane( const Vec3d& p, const Vec3d& normal, double c0, double e0, double r0 ){
    // https://en.wikipedia.org/wiki/Lennard-Jones_potential
    //printf(  " normal %g %g %g \n", normal.x, normal.y, normal.z );
    double cdot = normal.dot(p) - c0;
    double ir   = r0/cdot;
    double ir3  = ir*ir*ir;
    double f    = e0*(ir/r0)*ir3*(ir3-1);
    //printf( "%g %g %g %g %g %g %g \n", f, cdot, ir, ir3, e0, c0, r0  );
    return normal * f;
}

inline Vec3d getForceSpringRay( const Vec3d& p, const Vec3d& hray, const Vec3d& ray0, double k ){
    Vec3d dp; dp.set_sub( p, ray0 );
    double cdot = hray.dot(dp);
    dp.add_mul(hray,-cdot);
    return dp*k;
}

int pickParticle( int n, Vec3d * ps, Vec3d& ray0, Vec3d& hRay, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int i=0; i<n; i++){
        double ti = raySphere( ray0, hRay, R, ps[i] );
        if(ti<tmin){ imin=i; tmin=ti; }
    }
    return imin;
}

inline uint64_t getBondTypeId( uint16_t at1, uint16_t at2, uint8_t order ){
    if (at1>at2){ SWAP(at1,at2,uint16_t); }
    return pack64( at1, at2, order, 0 );
}

class BondType{ public:
    double length;
    double stiffness;
    uint16_t at1,at2;
    uint8_t  order;
};


// ==========  AtomType

class AtomType{ public:
    char      name[8];
    uint8_t   iZ;         // proton number
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint8_t   sym;        // sp, sp2, sp3 ... tetrahedral, triangle, linear, kink, square, octahedral
    uint32_t  color;
    double    RvdW;
    double    EvdW;

    char* toString( char * str ){
        sprintf( str, "printf: %s %i %i %i %i %lf %lf %x", name,  iZ,   neval,  valence,   sym,    RvdW, EvdW,   color );
        return str;
    }

    void fromString( char * str ){
        int iZ_, neval_, valence_, sym_;
        //char sclr[6];
        sscanf( str, " %s %i %i %i %i %lf %lf %x \n", name, &iZ_, &neval_, &valence_, &sym_,  &RvdW, &EvdW, &color );
        iZ=iZ_; neval=neval_; valence=valence_; sym=sym_;
        //printf( "printf: %s %i %i %i %i %lf %lf %x \n", name,  iZ,   neval_,  valence,   sym,    RvdW, EvdW,   color );
        //char ss[256]; printf("%s\n", toString(ss) );
    }

};

// ======================
// ====   MMFFparams
// ======================

class MMFFparams{ public:

    // http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html

    std::unordered_map<std::string,int> atypNames;
    std::vector<AtomType>               atypes;

    std::unordered_map<uint64_t,BondType> bonds;

    double default_bond_length      = 2.0;
    double default_bond_stiffness   = 1.0;

    int loadAtomTypes(char * fname){

        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;

        AtomType atyp;
        int i;
        for(int i; i<0xFFFF; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            atyp.fromString( line );
            atypes.push_back(atyp);
            atypNames[atyp.name] = atypes.size()-1;
        }
        return i;
    }

    int loadBondTypes(char * fname){

        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        BondType bt;
        //line = fgets( buff, 1024, pFile ); //printf("%s",line);
        //sscanf( line, "%i %i\n", &n );
        int i;
        for(int i; i<0xFFFF; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            //printf("%s",line);
            sscanf(  line, "%i %i %i %lf %lf\n", &bt.at1, &bt.at2, &bt.order, &bt.length, &bt.stiffness );
            printf(        "%i %i %i %lf %lf\n",  bt.at1,  bt.at2,  bt.order,  bt.length,  bt.stiffness );
            uint64_t id = getBondTypeId( bt.at1, bt.at2, bt.order );
            //printf( ":: (%i,%i,%i) -> %i \n", bt.at1, bt.at2, bt.order, id );
            //bt.at1--; bt.at2--;
            bonds[id]=bt;
        }
        return i;
    }

    bool getBondParams( int atyp1, int atyp2, int btyp, double& l0, double& k ){
        uint64_t id  = getBondTypeId( atypes[atyp1].iZ, atypes[atyp2].iZ, btyp );
        //printf( "(%i,%i,%i) -> %i \n", atypes[atyp1].iZ, atypes[atyp2].iZ, btyp, id );
        //uint64_t id  = getBondTypeId( atyp1, atyp2, btyp );
        auto it = bonds.find(id);
        if   ( it == bonds.end() ){ l0=default_bond_length; k=default_bond_stiffness; return false;}
        else                      { l0=it->second.length;   k=it->second.stiffness;   return true; }
    }

    void fillBondParams( int nbonds, Vec2i * bond2atom, int * bondOrder, int * atomType, double * bond_0, double * bond_k ){
        printf("fillBondParams: %i\n", nbonds);
        for(int i=0; i<nbonds; i++){
            Vec2i ib = bond2atom[i];
            getBondParams( atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
            //printf( "%i (%i %i) %i %g %g \n", i, atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
        }
    }

};

// ======================
// ====   GridFF
// ======================

class GridFF{ public:
    GridShape   grid;
    Vec3d  *FFPauli  = NULL;
    Vec3d  *FFLondon = NULL;
    Vec3d  *FFelec   = NULL;
    //Vec3d  *FFtot    = NULL; // total FF is not used since each atom-type has different linear combination

    int  natoms     = 0;
    int    * atypes = NULL;
    Vec3d  * apos   = NULL;   // atomic position
    Vec3d  * aREQs  = NULL;
    //Vec3d  * aPLQ = NULL;

    double alpha    = -1.6;


    void allocateFFs(){
        int ntot = grid.getNtot();
        //FFtot    = new Vec3d[ntot];
        FFPauli  = new Vec3d[ntot];
        FFLondon = new Vec3d[ntot];
        FFelec   = new Vec3d[ntot];
    }

    void allocateAtoms(int natoms_){
        natoms = natoms_;
        atypes = new int  [natoms];
        apos   = new Vec3d[natoms];
        aREQs  = new Vec3d[natoms];
    }

    int loadCell( char * fname ){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        fscanf( pFile, "%lf %lf %lf", &grid.cell.a.x, &grid.cell.a.y, &grid.cell.a.z );
        fscanf( pFile, "%lf %lf %lf", &grid.cell.b.x, &grid.cell.b.y, &grid.cell.b.z );
        fscanf( pFile, "%lf %lf %lf", &grid.cell.c.x, &grid.cell.c.y, &grid.cell.c.z );
        grid.updateCell();
    }

    int loadXYZ( char * fname, MMFFparams& params ){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        sscanf( line, "%i %i\n", &natoms );
        printf("%i \n", natoms );
        allocateAtoms(natoms);
        line = fgets( buff, 1024, pFile );
        for(int i=0; i<natoms; i++){
            char at_name[8];
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            double Q;
            int nret = sscanf( line, "%s %lf %lf %lf %lf\n", at_name, &apos[i].x, &apos[i].y, &apos[i].z, &Q );
            if( nret<5 ) Q = 0.0d;
            //printf(                  "%s %lf %lf %lf %lf\n", at_name,  apos[i].x,  apos[i].y,  apos[i].z,  Q );
            aREQs[i].z = Q;
            // atomType[i] = atomChar2int( ch );
            auto it = params.atypNames.find( at_name );
            if( it != params.atypNames.end() ){
                atypes[i] = it->second;
                aREQs[i].x = params.atypes[atypes[i]].RvdW;
                //aLJq[i].y = params.atypes[atypes[i]].EvdW;
                aREQs[i].y = sqrt( params.atypes[atypes[i]].EvdW );
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atypes[i] = -1;
                //aLJq[i].x  = 0.0d;
                //aLJq[i].x  = 0.0d;
            }
            printf(  "%i : >>%s<< %i (%g,%g,%g) %i (%g,%g,%g)\n", i, at_name, atypes[i],  apos[i].x,  apos[i].y,  apos[i].z, atypes[i], aREQs[i].x, aREQs[i].y, aREQs[i].z );
        }
        return natoms;
    }

    inline void addForce( Vec3d pos, Vec3d PLQ, Vec3d& f ){
        Vec3d gpos;
        grid.cartesian2grid(pos, gpos);
        f.add_mul( interpolate3DvecWrap( FFPauli,  grid.n, gpos ) , PLQ.x );
        f.add_mul( interpolate3DvecWrap( FFLondon, grid.n, gpos ) , PLQ.y );
        f.add_mul( interpolate3DvecWrap( FFelec,   grid.n, gpos ) , PLQ.z );
    }

    void init( Vec3i n, Mat3d cell, Vec3d pos0 ){
        grid.n     = n;
        grid.setCell(cell);
        grid.pos0  = pos0;
        allocateFFs();
    }

    void evalGridFFel(int natoms, Vec3d * apos, Vec3d * aREQs, Vec3d * FF ){
        //interateGrid3D( (Vec3d){0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p)->void{
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = (Vec3d){0.0,0.0,0.0};
            for(int ia=0; ia<natoms; ia++){ addAtomicForceQ( p-apos[ia], f, aREQs[ia].z ); }
            FF[ibuff]=f;
        });
    }

    void evalGridFFexp(int natoms, Vec3d * apos, Vec3d * aREQs, double alpha, double A, Vec3d * FF ){
        //interateGrid3D( (Vec3d){0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p){
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = (Vec3d){0.0,0.0,0.0};
            for(int ia=0; ia<natoms; ia++){
                //printf( " %i (%g,%g,%g) (%g,%g)\n", ia, apos[ia].x, apos[ia].y, apos[ia].z,  aLJq[ia].x, aLJq[ia].y  );
                addAtomicForceExp( p-apos[ia], f, aREQs[ia].x, aREQs[ia].y,    alpha );
                //addAtomicForceExp( p-apos[ia], f, aLJq[ia].x, aLJq[ia].y,    alpha*2 );
                //addAtomicForceExp( p-apos[ia], f, aLJq[ia].x, aLJq[ia].y*-2, alpha   );
            }
            //printf( " >> %i %i (%g,%g,%g) %g \n", ibuff, natoms, f.x, f.y, f.z, A  );
            FF[ibuff]=f*A;
            //printf( " %i (%g,%g,%g) \n", ibuff, p.x, p.y, p.z );
            //FF[ibuff]=p;
        });
    }

    void evalGridFFs(int natoms, Vec3d * apos, Vec3d * REQs ){
        //interateGrid3D( (Vec3d){0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p){
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d fp = (Vec3d){0.0d,0.0d,0.0d};
            Vec3d fl = (Vec3d){0.0d,0.0d,0.0d};
            Vec3d fe = (Vec3d){0.0d,0.0d,0.0d};
            for(int ia=0; ia<natoms; ia++){
                Vec3d dp; dp.set_sub( p, apos[ia] );
                Vec3d REQi = aREQs[ia];
                double r      = dp.norm();
                double ir     = 1/(r+RSAFE);
                double expar  = exp( alpha*(r-REQi.x) );
                double fexp   = alpha*expar*REQi.y*ir;
                fp.add_mul( dp, fexp*expar*2 );                    // repulsive part of Morse
                fl.add_mul( dp, fexp         );                    // attractive part of Morse
                fe.add_mul( dp, -14.3996448915d*REQi.z*ir*ir*ir ); // Coulomb
            }
            if(FFPauli)  FFPauli [ibuff]=fp;
            if(FFLondon) FFLondon[ibuff]=fl;
            if(FFelec)   FFelec  [ibuff]=fe;
        });
    }

    void evalGridFFs(int natoms, Vec3d * apos, Vec3d * REQs, Vec3i nPBC ){
        //interateGrid3D( (Vec3d){0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p){
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d fp = (Vec3d){0.0d,0.0d,0.0d};
            Vec3d fl = (Vec3d){0.0d,0.0d,0.0d};
            Vec3d fe = (Vec3d){0.0d,0.0d,0.0d};
            for(int ia=0; ia<natoms; ia++){
                Vec3d dp0; dp0.set_sub( p, apos[ia] );
                Vec3d REQi = aREQs[ia];
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    double r      = dp.norm();
                    double ir     = 1/(r+RSAFE);
                    double expar  = exp( alpha*(r-REQi.x) );
                    double fexp   = alpha*expar*REQi.y*ir;
                    fp.add_mul( dp, -fexp*expar*2 );                    // repulsive part of Morse
                    fl.add_mul( dp, -fexp         );                    // attractive part of Morse
                    fe.add_mul( dp, -14.3996448915d*REQi.z*ir*ir*ir ); // Coulomb
                }}}
            }
            if(FFPauli)  FFPauli [ibuff]=fp;
            if(FFLondon) FFLondon[ibuff]=fl;
            if(FFelec)   FFelec  [ibuff]=fe;
        });
    }

    void evalGridFFs( Vec3i nPBC){
        //evalGridFFexp( natoms, apos, aREQs, alpha*2,  1, FFPauli  );
        //evalGridFFexp( natoms, apos, aREQs, alpha  , 1, FFLondon );  // -2.0 coef is in  REQ2PLQ
        //evalGridFFel ( natoms, apos, aREQs,              FFelec   );
        //evalGridFFs( natoms, apos, aREQs );
        evalGridFFs( natoms, apos, aREQs, nPBC );
    }

    void evalCombindGridFF( Vec3d REQ, Vec3d * FF ){
        Vec3d PLQ = REQ2PLQ( REQ, alpha );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = (Vec3d){0.0d,0.0d,0.0d};
            if(FFPauli ) f.add_mul( FFPauli[ibuff],  PLQ.x );
            if(FFLondon) f.add_mul( FFLondon[ibuff], PLQ.y );
            if(FFelec  ) f.add_mul( FFelec[ibuff],   PLQ.z );
            FF[ibuff] =  f;
        });
    }

    void evalCombindGridFF_CheckInterp( Vec3d REQ, Vec3d * FF ){
        Vec3d PLQ = REQ2PLQ( REQ, alpha );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = (Vec3d){0.0d,0.0d,0.0d};
            //addForce( p, PLQ, f );
            addForce( p+(Vec3d){0.1,0.1,0.1}, PLQ, f );
            FF[ibuff] =  f;
        });
    }

    void evalFFline( int n, Vec3d p0, Vec3d p1, Vec3d PLQ, Vec3d * pos, Vec3d * fs ){
        Vec3d dp = p1-p0; dp.mul(1.0d/(n-1));
        Vec3d  p = p0;
        for(int i=0; i<n; i++){
            if(fs ){
                Vec3d fp = (Vec3d){0.0,0.0,0.0};
                Vec3d fl = (Vec3d){0.0,0.0,0.0};
                Vec3d fq = (Vec3d){0.0,0.0,0.0};
                for(int ia=0; ia<natoms; ia++){
                    Vec3d d = p-apos[ia];
                    addAtomicForceExp( d, fp, aREQs[ia].x, aREQs[ia].y, 2*alpha );
                    addAtomicForceExp( d, fl, aREQs[ia].x, aREQs[ia].y,   alpha );
                    //addAtomicForceQ  ( d, fq, aLJq[ia].z );
                    //printf( "%i %i %g  (%g,%g)   %g \n", i, ia, d.z, aLJq[ia].x, aLJq[ia].y, alpha  );
                }
                fs[i] = fp*PLQ.x + fl*PLQ.y + fq*PLQ.z;
            }
            if(pos)pos[i]=p;
            //printf("%i %20.10f %20.10f %20.10f    %20.10e %20.10e %20.10e\n" , i, pos[i].x, pos[i].y, pos[i].z, fs[i].x, fs[i].y, fs[i].z );
            p.add(dp);
        }
    }

    void evalFFlineToFile( int n, Vec3d p0, Vec3d p1, Vec3d REQ, const char * fname ){
        Vec3d * pos = new Vec3d[n];
        Vec3d * fs  = new Vec3d[n];
        evalFFline( n, p0, p1, REQ2PLQ(REQ, alpha), pos, fs );
        FILE* pfile;
        pfile = fopen(fname, "w" );
        for(int i=0; i<n; i++){
            fprintf( pfile, "%i %20.10f %20.10f %20.10f    %20.10e %20.10e %20.10e\n", i, pos[i].x, pos[i].y, pos[i].z, fs[i].x, fs[i].y, fs[i].z);
        }
        fclose(pfile);
        delete [] pos; delete [] fs;
    }

}; // RigidSubstrate


// ======================
// ====   MMFF
// ======================

class MMFF{ public:
    int  natoms=0, nbonds=0, nang=0, ntors=0;

    Vec2i  * bond2atom = NULL;
    double * bond_0    = NULL;  // [A]
    double * bond_k    = NULL;  // [eV/A] ?

    Vec2i  * ang2bond  = NULL;
    Vec3i  * ang2atom  = NULL;
    Vec2d  * ang_0     = NULL; // [1]
    double * ang_k     = NULL; // [eV/A^2]

    Vec3i  * tors2bond = NULL;
    Quat4i * tors2atom = NULL;
    Vec2d  * tors_0    = NULL; // [1]
    double * tors_k    = NULL; // [eV/A^2]

    // molecular gemeotry
    int    * atypes = NULL;
    Vec3d  * apos   = NULL;   // atomic position
    Vec3d  * aLJq   = NULL;
    Vec3d  * aPLQ   = NULL;   // this is used for grid-accelerated factorized potential
    double * lbond  = NULL;   // bond lengths
    Vec3d  * hbond  = NULL;   // normalized bond unitary vectors
    Vec3d  * aforce = NULL;

    //RigidSubstrate substrate;
    GridFF gridFF;


void allocate( int natoms_, int nbonds_, int nang_, int ntors_ ){
    natoms=natoms_; nbonds=nbonds_; nang=nang_; ntors=ntors_;
    if(atypes   ==NULL) atypes    = new int   [natoms];
    if(apos     ==NULL) apos      = new Vec3d [natoms];
    if(aforce   ==NULL) aforce    = new Vec3d [natoms];
    if(aLJq     ==NULL) aLJq      = new Vec3d [natoms];

    if(lbond    ==NULL) lbond     = new double[nbonds];
    if(hbond    ==NULL) hbond     = new Vec3d [nbonds];

    if(bond2atom==NULL) bond2atom = new Vec2i [nbonds];
    if(bond_0   ==NULL) bond_0    = new double[nbonds];
    if(bond_k   ==NULL) bond_k    = new double[nbonds];

    if(ang2bond ==NULL) ang2bond  = new Vec2i [nang];
    if(ang2atom ==NULL) ang2atom  = new Vec3i [nang];
    if(ang_0    ==NULL) ang_0     = new Vec2d [nang];
    if(ang_k    ==NULL) ang_k     = new double[nang];

    /*
    tors2bond = new Vec3i [ntors];
    tors2atom = new Quat4i[ntors];
    tors_0    = new double[ntors];
    tors_k    = new double[ntors];
    */
}

void deallocate(){
    delete [] apos;     delete [] aforce;   delete [] aLJq;
    delete [] lbond;    delete [] hbond;    delete [] bond2atom; delete [] bond_0; delete [] bond_k;
    delete [] ang2bond; delete [] ang2atom; delete [] ang_0;     delete [] ang_k;
    if(aPLQ) delete [] aPLQ;
}

int pickBond( Vec3d& ray0, Vec3d& hRay, double R ){
    double dist_min =  R;
    int    imin     = -1;
    for(int ib=0; ib<nbonds; ib++){
        Vec2i iat = bond2atom[ib];
        double t1,t2;
        double dist = rayLine( ray0, hRay, apos[iat.x], hbond[ib], t1, t2 );
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLine( ray0+(hRay*t1), apos[iat.x]+(hbond[ib]*t2) );
        //printf( "%i %g %g %g \n", ib, t1, t2, dist );
        if( (dist<dist_min) && (t2>0) && (t2<lbond[ib]) ){
            imin=ib; dist_min=dist;
        }
    }
    return imin;
}

void genPLQ(){
    if(aPLQ==NULL) aPLQ = new Vec3d[natoms];
    for(int i=0; i<natoms; i++){
        aPLQ[i] = REQ2PLQ( aLJq[i], gridFF.alpha );
        printf( "genPLQ %i (%g,%g,%g)->(%g,%g,%g) \n", i, aLJq[i].x, aLJq[i].y, aLJq[i].z,   aPLQ[i].x, aPLQ[i].y, aPLQ[i].z );
    }
    //exit(0);
}

void translate( Vec3d dpos){ for(int i=0; i<natoms; i++) apos[i].add(dpos); }

void ang_b2a(){
    for(int i=0; i<nang; i++){
        Vec2i ib = ang2bond[i];
        Vec2i b1,b2;
        b1 = bond2atom[ib.x];
        b2 = bond2atom[ib.y];
        if     ( b1.x == b2.x ){ ang2atom[i] = (Vec3i){ b1.y, b2.y, b1.x }; }
        else if( b1.x == b2.y ){ ang2atom[i] = (Vec3i){ b1.y, b2.x, b1.x }; ang2bond[i].y|=SIGN_MASK; }
        else if( b1.y == b2.x ){ ang2atom[i] = (Vec3i){ b1.x, b2.y, b1.y }; ang2bond[i].x|=SIGN_MASK; }
        else if( b1.y == b2.y ){ ang2atom[i] = (Vec3i){ b1.x, b2.x, b1.y }; ang2bond[i].x|=SIGN_MASK; ang2bond[i].y|=SIGN_MASK; }
        //printf( "%i (%i,%i)(%i,%i) (%i,%i,%i) (%i,%i) (%i,%i) \n", i, b1.x+1,b1.y+1, b2.x+1,b2.y+1, ang2atom[i].x+1,ang2atom[i].y+1,ang2atom[i].z+1,  ang2bond[i].x&0xFFFF, ang2bond[i].y&0xFFFF, ang2bond[i].x, ang2bond[i].y );
    }
}

void eval_bonds( const bool substract_LJq ){
    for(int ib=0; ib<nbonds; ib++){
        Vec2i iat = bond2atom[ib];
        Vec3d dp; dp.set_sub( apos[iat.y], apos[iat.x] );
        Vec3d f;
        //printf( "%i %i (%g,%g,%g) (%g,%g,%g) \n", iat.x, iat.y, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
        double l = dp.norm();
        f.set_mul(dp,1.0d/l);
        //printf( " %i (%i,%i) (%g,%g,%g) %g \n", ib, iat.x, iat.y, dp.x, dp.y, dp.z, l );

        lbond [ib] = l;
        hbond [ib] = f;

        f.mul( (l-bond_0[ib])*bond_k[ib] );

        if( substract_LJq ){
            double rij = aLJq[iat.x].x+aLJq[iat.y].x;
            double eij = aLJq[iat.x].y*aLJq[iat.y].y;
            double qq  = aLJq[iat.x].z*aLJq[iat.y].z;
            addAtomicForceLJQ( dp, f, rij, -eij, qq );
        }

        aforce[iat.x].add( f );
        aforce[iat.y].sub( f );
    }
    //exit(0);
}

void eval_angcos(){
    // simpler and faster angular forces cos(theta)
    // v      = -k*cos(a234) = k*dot(h23,h34)
    // dv/dr2 = ( h34 - dot(h23,h34)*h23 )/r23
    // dv/dr4 = ( dot(h23,h34)*h34 - h23 )/r34
    for(int ig=0; ig<nang; ig++){
        Vec2i ib = ang2bond[ig];
        Vec3i ia = ang2atom[ig];
        //Vec3d h1 = hbond[ib.x];
        //Vec3d h2 = hbond[ib.y];
        Vec3d h1,h2;

        if(ib.x&SIGN_MASK){ ib.x&=0xFFFF; h1 = hbond[ib.x]; h1.mul(-1.0d); }else{ h1 = hbond[ib.x]; };
        if(ib.y&SIGN_MASK){ ib.y&=0xFFFF; h2 = hbond[ib.y]; h2.mul(-1.0d); }else{ h2 = hbond[ib.y]; };

        /*
        Vec3d h1_,h2_;
        h1_ = apos[ia.x] - apos[ia.z]; h1_.normalize();
        h2_ = apos[ia.y] - apos[ia.z]; h2_.normalize();
        printf( " ig %i c1 %g c2 %g   (%i,%i) \n", ig, h1_.dot(h1), h2_.dot(h2),   ang2bond[ig].x,  ang2bond[ig].y );
        */

        double c = h1.dot(h2);
        //double s = sqrt(1-c*c);

        //printf( " %i (%i,%i) (%g,%g,%g) (%g,%g,%g) \n", ib, ib.x,ib.y, h1.x,h1.y,h1.z, h2.x,h2.y,h2.z );

        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = h2 - h1*c;
        hf2 = h1 - h2*c;
        //hf1 = h1*c - h2;
        //hf2 = h2*c - h1;

        //printf("ig %i ib (%i,%i) \n", ig, ib.x, ib.y );

        double fang = -ang_k[ig]/(1.02-c);
        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fang/lbond[ib.y] );

        //printf("ia (%i,%i,%i)\n", ia.x, ia.y, ia.z );

        /*
        glColor3f(0.0f,1.0f,0.0f);
        Draw3D::drawVecInPos( h1*0.25, apos[ia.z] );
        Draw3D::drawVecInPos( h2*0.25, apos[ia.z] );

        glColor3f(1.0f,0.0f,0.0f);
        Draw3D::drawVecInPos( hf1, apos[ia.x] );
        Draw3D::drawVecInPos( hf2, apos[ia.y] );
        */

        aforce[ia.x].add( hf1     );
        aforce[ia.y].add( hf2     );
        aforce[ia.z].sub( hf1+hf2 );

        //printf("ig %i\n", ig );
        //aforce[ia.x] ????
        // TODO : zero moment condition
    }

    /*
    evalGridFFs(int natoms, Vec3d * apos, Vec3d * aLJqs ){
        // TODO : this will be faster since we may reuse dpos, sqrt() function etc.
    }
    */

}

void eval_angles(){
    for(int ig=0; ig<nang; ig++){
        Vec2i ib = ang2bond[ig];
        Vec3d h1 = hbond[ib.x];
        Vec3d h2 = hbond[ib.y];
        Vec2d cs;
        cs.x = h1.dot(h2);
        cs.y = sqrt(1-cs.x*cs.x);

        Vec3d hf1,hf2; // unitary vectors of force
        double inv_sa = 1.0d/cs.y;
        hf1 = (h2 - h1*cs.x); // *inv_sa -> moved to fang for faster evaluation
        hf2 = (h1 - h2*cs.x);

        cs.mul_cmplx( ang_0[ig] );
        //E = 0.5*ang_k[ig]*cs.x*cs.x;
        double fang = (ang_k[ig]*cs.x*cs.y)*inv_sa;

        Vec3i ia = ang2atom[ig];

        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fang/lbond[ib.y] );

        aforce[ia.y].add( hf1     );
        aforce[ia.z].add( hf2     );
        aforce[ia.x].sub( hf1+hf2 );
        // TODO : check zero moment condition ... probably OK since pivot atom is in r=0 (COG, no torque by definition)
    }
}

// torsions should be special - for whole group - more atoms
void eval_torsion(){
    for(int it=0; it<ntors; it++){
        Vec3i ib = tors2bond[it];
        Vec3d h1 = hbond[ib.x];
        Vec3d h2 = hbond[ib.y];
        Vec3d ax = hbond[ib.z];

        // --- outproject axis
        h1.add_mul( ax, ax.dot(h1) );
        h2.add_mul( ax, ax.dot(h2) );

        Vec2d cs;
        cs.x = h1.dot(h2);
        cs.y = sqrt(1-cs.x*cs.x);

        double inv_sa = 1.0d/cs.y;
        Vec3d hf1,hf2; // unitary vectors of force
        hf1 = (h1 - h2*cs.x)*inv_sa;
        hf2 = (h2 - h1*cs.x)*inv_sa;

        cs.mul_cmplx( tors_0[it] );

        // TODO : torsion should be rather some fourier series or spline ? not just one rotation angle (?)
        double fang = ang_k[it] * cs.x * cs.y;

        Quat4i ia = tors2atom[it];

        hf1.mul( fang/lbond[ib.x] );
        hf2.mul( fang/lbond[ib.y] );

        aforce[ia.y].add( hf1 );
        aforce[ia.z].add( hf2 );
        aforce[ia.x].sub( hf1 ); // is this correct ?
        aforce[ia.w].sub( hf2 ); // is this correct ?
        //aforce[ia.x] ????
        // TODO : zero moment condition
    }
}

void eval_LJq_On2(){
    for(int i=0; i<natoms; i++){
        Vec3d ljq_i = aLJq[i];
        Vec3d pi    = apos[i];
        Vec3d f; f.set(0.0);
        for(int j=0; j<natoms; j++){
            if(i!=j){
                Vec3d& ljq_j = aLJq[j];
                double rij = ljq_i.x+ljq_j.x;
                double eij = ljq_i.y*ljq_j.y;
                double qq  = ljq_i.z*ljq_j.z;
                addAtomicForceLJQ( pi-apos[j], f, rij, -eij, qq );
            }
        }
        aforce[i].add(f);
    }
}

void eval_FFgrid(){
    for(int i=0; i<natoms; i++){
        gridFF.addForce( apos[i], aPLQ[i], aforce[i] );
    }
}

/*
void eval_spring( Vec3d hray, double k ){

};
*/

void printBondParams(){
    for( int i=0; i<nbonds; i++ ){
        printf( "%i (%i,%i) %g %g \n", i, bond2atom[i].x+1, bond2atom[i].y+1, bond_0[i], bond_k[i] );
    }
}

}; // MMFF


// ======================
// ====   MMFFBuilder
// ======================

#include "Molecule.h"

class MMFFAtom{ public:
    int type;
    Vec3d pos;
    Vec3d LJq;
};

class MMFFBond{ public:
    int    type;
    Vec2i  atoms;
    //double l0,k;
};

class MMFFAngle{ public:
    int type;
    Vec2i  bonds;
    //double a0,k;
};

class MMFFmol{ public:
    Molecule * mol;
    Vec3i    i0;
};

class MMFFBuilder{  public:
    std::vector<MMFFAtom>  atoms;
    std::vector<MMFFBond>  bonds;
    std::vector<MMFFAngle> angles;
    std::vector<MMFFmol>   mols;

    void insertMolecule( Molecule * mol, Vec3d pos, Mat3d rot ){
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        int nangle0 = angles.size();
        mols.push_back( (MMFFmol){mol, (Vec3i){natom0,nbond0,nangle0} } );
        for(int i=0; i<mol->natoms; i++){
            //Vec3d LJq = (Vec3d){0.0,0.0,0.0}; // TO DO : LJq can be set by type
            Vec3d LJq = (Vec3d){1.0,0.03,0.0}; // TO DO : LJq can be set by type
            Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
            atoms.push_back( (MMFFAtom){mol->atomType[i], p, LJq } );
        }
        for(int i=0; i<mol->nbonds; i++){
            bonds.push_back( (MMFFBond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}) } );
        }
        for(int i=0; i<mol->nang; i++){
            angles.push_back( (MMFFAngle){ 1, mol->ang2bond[i] + ((Vec2i){nbond0,nbond0}) } );
        }
    }


    void assignAtomTypes( MMFFparams * params ){
        for(int i=0; i<atoms.size(); i++){
            //mmff->aLJq [i]  = atoms[i].type;
            int ityp = atoms[i].type;
            atoms[i].LJq.x = params->atypes[ityp].RvdW;
            atoms[i].LJq.y = params->atypes[ityp].EvdW;
            atoms[i].LJq.z = 0;
            //atomTypes[i]  = atoms[i].type;
        }
    }

    void toMMFF( MMFF * mmff, MMFFparams * params ){
        mmff->deallocate();
        mmff->allocate( atoms.size(), bonds.size(), angles.size(), 0 );
        //int * atomTypes = new int[atoms.size()];
        //int * bondTypes = new int[bonds.size()];
        for(int i=0; i<atoms.size(); i++){
            mmff->atypes[i] = atoms[i].type;
            mmff->apos [i]  = atoms[i].pos;
            mmff->aLJq [i]  = atoms[i].LJq;
            //atomTypes[i]  = atoms[i].type;
        }
        for(int i=0; i<bonds.size(); i++){
            mmff->bond2atom[i] = bonds[i].atoms;
            Vec2i ib           = bonds[i].atoms;
            params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, mmff->bond_0[i], mmff->bond_k[i] );
            //bondTypes[i]       = bonds[i].type;
        }
        for(int i=0; i<angles.size(); i++){
            mmff->ang2bond[i] = angles[i].bonds;
            mmff->ang_0[i] = {1.0,0.0}; // TODO FIXME
            mmff->ang_k[i] = 0.5;       // TODO FIXME
        }
        //params.fillBondParams( mmff->nbonds, mmff->bond2atom, bondTypes, atomTypes, mmff->bond_0, mmff->bond_k );
        //delete [] atomTypes;
        //delete [] bondTypes;
    }
};

int write2xyz( FILE* pfile, MMFF * mmff, MMFFparams * params ){
    fprintf(pfile, "%i\n", mmff->natoms );
    fprintf(pfile, "#comment \n");
    for(int i=0; i<mmff->natoms; i++){
        int ityp   = mmff->atypes[i];
        Vec3d&  pi = mmff->apos[i];
        printf( "%i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
        fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n", params->atypes[ityp].name, pi.x,pi.y,pi.z );
    };
    return mmff->natoms;
}

int save2xyz( char * fname, MMFF * mmff, MMFFparams * params ){
    FILE* pfile = fopen(fname, "w");
    if( pfile == NULL ) return -1;
    int n = write2xyz(pfile, mmff, params );
    fclose(pfile);
    return n;
}

#endif
