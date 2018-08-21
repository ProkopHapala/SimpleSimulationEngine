#ifndef GridFF_h
#define GridFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Grid.h"
#include "Forces.h"
#include "MMFFparams.h"

class GridFF{ public:
    GridShape   grid;
    Vec3d  *FFPauli    = NULL;
    Vec3d  *FFLondon   = NULL;
    Vec3d  *FFelec     = NULL;
    
    Quat4f *FFPauli_f  = NULL;
    Quat4f *FFLondon_f = NULL;
    Quat4f *FFelec_f   = NULL;
    
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
    
    void setAtoms( int natoms_, Vec3d * apos_, Vec3d * REQs_ ){
        natoms = natoms_;
        //atypes = new int  [natoms];
        apos   = apos_;
        aREQs  = REQs_;
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
            for(int iat=0; iat<natoms; iat++){
                Vec3d dp0; dp0.set_sub( p, apos[iat] );
                Vec3d REQi = aREQs[iat];
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    double r      = dp.norm();
                    double ir     = 1/(r+RSAFE);
                    double expar  = exp( alpha*(r-REQi.x) );
                    double fexp   = alpha*expar*REQi.y*ir;
                    fp.add_mul( dp, -fexp*expar*2 );                    // repulsive part of Morse
                    fl.add_mul( dp, -fexp         );                    // attractive part of Morse
                    fe.add_mul( dp, 14.3996448915d*REQi.z*ir*ir*ir ); // Coulomb
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


#endif
