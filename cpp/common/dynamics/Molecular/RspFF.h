
/*

Rigid Atom sp-hybridization forcefield

*/

#ifndef RARFFarr_h
#define RARFFarr_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"


#define N_BOND_MAX 4
#define R1SAFE    1e-8

static const double sp3_hs[] = {
-0.57735026919, -0.57735026919, -0.57735026919,
+0.57735026919, +0.57735026919, -0.57735026919,
-0.57735026919, +0.57735026919, +0.57735026919,
+0.57735026919, -0.57735026919, +0.57735026919
};

static const double sp2_hs[] = {
+1.000000000000, -0.00000000000,  0.00000000000,
-0.500000000000, +0.86602540378,  0.00000000000,
-0.500000000000, -0.86602540378,  0.00000000000,
 0.00000000000,   0.00000000000,  1.00000000000     // electron - can be repulsive
};

static const double sp1_hs[] = {
+1.000000000000,  0.00000000000,  0.00000000000,
-1.000000000000,  0.00000000000,  0.00000000000,
 0.000000000000,  1.00000000000,  0.00000000000,    // electron - can be repulsive
 0.000000000000,  0.00000000000,  1.00000000000     // electron - can be repulsive
};

template<typename T>
void rotateVectors(int n, const Quat4TYPE<T>& qrot, Vec3TYPE<T>* h0s, Vec3TYPE<T>* hs ){
    Mat3TYPE<T> mrot;
    qrot.toMatrix(mrot);
    for( int j=0; j<n; j++ ){
        Vec3TYPE<T> h;
        //mrot.dot_to_T    ( h0s[j], h );
        mrot.dot_to      ( h0s[j], h );
        hs[j] = h;
        //ps[j].set_add_mul( pos, p_, r0 );
    }
}

/*
struct RigidAtomType{
    int    nbond = 4;  // number bonds

    double rbond0 =  0.5;
    double aMorse =  4.0;
    double bMorse = -0.7;

    double Epz    =  0.5;
    double c6     = -15.0;
    double R2vdW  =  8.0;
    Vec3d* bh0s = (Vec3d*)sp3_hs;

    inline void combine(const RigidAtomType& a, const RigidAtomType& b ){
        nbond   = a.nbond;
        aMorse  = a.aMorse * b.aMorse;  // TODO
        bMorse  = a.bMorse + b.bMorse;
        rbond0  = a.rbond0 + b.rbond0;

        Epz   = a.Epz * b.Epz;
        c6    = -(a.c6  * b.c6);
        R2vdW = a.R2vdW + a.R2vdW;
    }

    void print(){
        printf( "nbond  %i rbond0 %g\n", nbond, rbond0 );
        printf( "aMorse %g bMorse %g\n", aMorse, bMorse );
        printf( "c6     %g r2vdW  %g\n", c6, R2vdW );
        //exit(0);
    }
};
*/

class RARFF2arr{ public:

    double invRotMass = 2.0;

    int natom = 0;
    int nbond = 0;

    // atom properties
    //RigidAtomType** types = 0;

    int   * types   = 0;
    Vec3d * aREQs   = 0;

    Vec2i  * b2a    = 0;
    Vec2d  * bLKs   = 0;

    Vec3d  * poss   = 0;
    Vec3d  * forces = 0;
    Vec3d  * vels   = 0;

    Quat4d * qrots  = 0;
    Vec3d  * torqs  = 0;
    Vec3d  * omegas = 0;

    // aux
    //double *        ebonds=0;
    Vec3d  *        hbonds=0;
    Vec3d  *        fbonds=0;

    //double F2pos=0;
    //double F2rot=0;

    void realloc(int natom_){
        natom=natom_;
        //_realloc(atoms,natom   );
        _realloc(types  ,natom);
        _realloc(aREQs  ,natom);

        _realloc(bLKs   ,nbond);
        _realloc(b2a    ,nbond);

        _realloc(poss   ,natom);
        _realloc(forces ,natom);
        _realloc(omegas ,natom);

        _realloc(qrots  ,natom);
        _realloc(torqs  ,natom);
        _realloc(vels   ,natom);

        //_realloc(ebonds ,nbond);
        _realloc(hbonds , natom*N_BOND_MAX);
        _realloc(fbonds , natom*N_BOND_MAX);
    }

    void dealloc(){

        _dealloc(types  );
        _dealloc(aREQs  );

        _dealloc(bLKs   );
        _dealloc(b2a    );

        _dealloc(poss   );
        _dealloc(forces );
        _dealloc(omegas );

        _dealloc(qrots  );
        _dealloc(torqs  );
        _dealloc(vels   );

        //_dealloc(ebonds );
        _dealloc(hbonds );
        _dealloc(fbonds );
    }

    inline double bondEF( int ibond ){

        int ib = b2a[ibond].a;
        int jb = b2a[ibond].b;
        int ja = ib>>2;
        int ja = jb>>2;

        Vec3d  dij = poss[ja] - poss[ia];
        double r2  = dij.norm2() + R2SAFE;
        double rij = sqrt( r2 );
        Vec3d  hij = dij*(1/rij);

        Vec2d& lk = bLKs[ibond];

        double dl    = (rij-lk.x);
        double eb    = lk.y*dl;
        E           += eb;
        double fr    = lk.y*2*dl;

        Vec3d force=hij*fr;

        if( substract_LJq ){
            addAtomicForceLJQ( dp, f, aREQ[iat.x].x+aREQ[iat.y].x, -aREQ[iat.x].y*aREQ[iat.y].y, aREQ[iat.x].z*aREQ[iat.y].z );
            //addAtomicForceMorseQ( dp, f, aREQ[iat.x].x+aREQ[iat.y].x, -aREQ[iat.x].y*aREQ[iat.y].y, aREQ[iat.x].z*aREQ[iat.y].z, gridFF.alpha );
        }

        const Vec3d& hi = hbonds[ib];
        Vec3d&       fi = fbonds[ib];

        const Vec3d& hj = hbonds[jb];
        Vec3d&       fj = fbonds[jb];

        double ci       = hij.dot( hi );
        double cj       = hij.dot( hj );  // cj  = <hj|hij>
        //double cij      = hi .dot( hj );  // cij = <hj|hi>

        //if( (cj>0)||(cij>0) ) continue;
        // Erot = <hi,hij> - <hj,hij>
        // E = cos(x)                   : derivative vanish at x=pi
        // E     (1-sqrt( (1+cos(x))/2 ))*4 : correct this problem
        // dE =  dcos(x)/(sqrt( (1+cos(x)/2))

        double fci = Kss/sqrt((1+ci)*0.5);
        double fcj = Kss/sqrt((1+cj)*0.5);

        //double cc  = ci*cj*cij;
        //double cc2 = cc*cc;
        //double eb   = (ci - cj) * Kss;
        //double de  = 4*cc2*cc;

        //double eEb = e * Eb;
        //eis[ib]+=eEb*0.5;
        //ejs[jb]+=eEb*0.5;
        //E += eb;

        double ir3 = 1/(rij*r2);
        double xx  = dij.x*dij.x;
        double yy  = dij.y*dij.y;
        double zz  = dij.z*dij.z;

        double dxy=-dij.x*dij.y*ir3;
        double dxz=-dij.x*dij.z*ir3;
        double dyz=-dij.y*dij.z*ir3;
        double dxx=(yy+zz)*ir3;
        double dyy=(xx+zz)*ir3;
        double dzz=(xx+yy)*ir3;

        force.x += fci*( hi.x*dxx + hi.y*dxy + hi.z*dxz ) - fcj*( hj.x*dxx + hj.y*dxy + hj.z*dxz );
        force.y += fci*( hi.x*dxy + hi.y*dyy + hi.z*dyz ) - fcj*( hj.x*dxy + hj.y*dyy + hj.z*dyz );
        force.z += fci*( hi.x*dxz + hi.y*dyz + hi.z*dzz ) - fcj*( hj.x*dxz + hj.y*dyz + hj.z*dzz );

        fi.add_mul(hij,fci);
        fj.sub_mul(hij,fcj);

        //printf( "%i,%i nbi %i nbj %i \n", ia, ja, nbi, nbj );
        if( (nbi==3)&&(nbj==3) ){ // align pz-pz in sp2-sp2 pair
            const Vec3d& hi = his[3];
            const Vec3d& hj = hjs[3];
            Vec3d&       fi = fis[3];
            Vec3d&       fj = fjs[3];
            double cdot = hi.dot(hj);
            //double E    = Epz * ( cdot*cdot );
            double de = -2*cdot*type.Epz*Eb;
            //printf( "cdot %g Epz %g de %g \n", cdot, type.Epz, de );
            fi.add_mul(hj,de);
            fj.add_mul(hi,de);
            //force      // TODO: proper derivatives of energy
        }

        forces[ja].sub(force);
        forces[ia].add(force);
        return E;
    }

    double interEF(){
        double E = 0;
        for(int i=0; i<nbond; ++){ bondEF( i ); }
        /*
        for(int i=0; i<natom; i++){
            RigidAtomType& typei = *types[i];
            Vec3d           pi   = poss[i];
            for(int j=i+1; j<natom; j++){
            //for(int j=0; j<natom; j++){
                RigidAtomType pairType;
                pairType.combine( typei, *types[j] );
                E += pairEF( i, j, typei.nbond, types[j]->nbond, pairType );
            }
        }
        */
        return E;
    }

    void eval_LJq_On2(){
        for(int i=0; i<natoms; i++){
            const Vec3d& ljq_i = aREQ[i];
            const Vec3d& pi    = apos[i];
            Vec3d f; f.set(0.0);
            for(int j=0; j<natoms; j++){
                if(i!=j){ //  ToDo : can be up to twice faster if we do not do all pairs
                    const Vec3d& ljq_j = aREQ[j];
                    double rij = ljq_i.x+ljq_j.x;
                    double eij = ljq_i.y*ljq_j.y;
                    double qq  = ljq_i.z*ljq_j.z;
                    addAtomicForceLJQ( pi-apos[j], f, rij, -eij, qq );
                }
            }
            aforce[i].add(f);
        }
    }

    void cleanAux(){
        for(int i=0; i<natom; i++){
            vels  [i]=Vec3dZero;
            omegas[i]=Vec3dZero;
            torqs [i]=Vec3dZero;
            forces[i]=Vec3dZero;
        }
        int nb   = natom*N_BOND_MAX;
        for(int i=0; i<nb;   i++){
            ebonds  [i]=0;
            bondCaps[i]=-1;
            fbonds  [i].set(0.0);
        }
    }

    void cleanAtomForce(){
        for(int i=0; i<natom; i++){
            //atoms[i].cleanForceTorq();
            torqs [i]=Vec3dZero;
            forces[i]=Vec3dZero;
        }
        int nb   = natom*N_BOND_MAX;    for(int i=0; i<nb;   i++){ ebonds[i]=0; }
        //printf("\n"); for(int i=0; i<nb;   i++){ printf("%i ", bondCaps[i] ); }; printf("\n");
        int nval = natom*N_BOND_MAX*3;  for(int i=0; i<nval; i++){ ((double*)fbonds)[i]=0;}
    }

    double evalF2rot(){ double F2=0; for(int i=0; i<natom; i++){ F2+=torqs [i].norm2(); }; return F2; }
    double evalF2pos(){ double F2=0; for(int i=0; i<natom; i++){ F2+=forces[i].norm2(); }; return F2; }

    double projectBonds(){
        for(int i=0; i<natom; i++){
            rotateVectors( N_BOND_MAX, qrots[i], types[i]->bh0s, hbonds + i*N_BOND_MAX );
        }
    }

    void evalTorques(){
        for(int ia=0; ia<natom; ia++){
            //RigidAtom& atomi = atoms[ia];
            //Vec3d torq = Vec3dZero;
            //int nbi =  types[ia]->nbond;
            //for(int ib=0; ib<nbi; ib++){
            for(int ib=0; ib<N_BOND_MAX; ib++){
                int io = 4*ia+ib;
                fbonds[io].makeOrthoU(hbonds[io]);
                //printf( "ia %i ib %i f(%g,%g,%g)\n", ia, ib,  fbonds[io].x,fbonds[io].y,fbonds[io].z );
                torqs[ia].add_cross( hbonds[io], fbonds[io] );
            }
        }
    }

    void move(double dt){
        for(int i=0; i<natom; i++){
            //atoms[i].moveRotGD(dt*invRotMass);
            //atoms[i].movePosGD(dt);
            qrots[i].dRot_exact( dt, torqs[i] );  qrots[i].normalize();          // costly => debug
            //dRot_taylor2( dt, torqs[i] );   qrots[i].normalize_taylor3();  // fast => op
            poss[i].add_mul(forces[i],dt);
        }
    }

    void moveMDdamp(double dt, double damp){
        for(int i=0; i<natom; i++){
            //atoms[i].moveMDdamp( dt, invRotMass, damp);
            vels  [i].set_lincomb( damp, vels  [i], dt,            forces[i] );
            omegas[i].set_lincomb( damp, omegas[i], dt*invRotMass, torqs [i] );
            poss  [i].add_mul    ( vels[i], dt    );
            qrots [i].dRot_exact ( dt,  omegas[i] );
            qrots [i].normalize();
        }
    }

};

#endif
