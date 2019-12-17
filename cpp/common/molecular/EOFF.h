
#ifndef EOFF_h
#define EOFF_h


/*

Electron Orbital Force-Field
============================

General Iddea:

 - each atom has 1 core and 4 valence orbitals (like sp3 tetraherdon) which can be occupied by 0,1,2 electrons

 There are several kinds of forces:
  - intra-atomic
    - core-orb - keep  certain distance from core (radial spring, r2-ir2 )
    - orb-orb  - keeps orbitals spread far from each other
        - can depend on orbital occupation (?)
 - inter-atomic
    - core-core   repulsion ( perhaps A*exp(b*R) )
    - orb-orb     attraction - two orbitals from different atoms attract each other (tries to maximize overlap, share electros)
    - orb|orb-orb repulsion  - two orbitals which are close to each other repel all other orbitals in range (exclusion principle)


  orb|orb-orb  derivs:

    E     = exp( bij * lij ) * exp(  bijk * (lik+ljk)/2 ) = exp( bij*lij  +  bijk_*(lik+ljk) )
    dE/dlij = E * bij ;
    dE/dlik = E * bijk;
    dE/dljk = E * bijk;

    dlij/dxi = (xi-xj)/lij


*/

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
//#include "GridFF.h"

#define N_ATOM_ORB 4

struct EAFFelement{
    // repulsion between cores  E_core = core_A * exp( core_b * R)
    double E0;      // free atom energy
    double core_b;
    double core_A;

    double oc_r0 = 0.7;  // distance of orbitals from center of atom
    double oc_ir2;       // interaction between core and orbitals
    double oc_r2;
    double oo_ir2;       // interaction between core and orbitals
    //double C6;
};

inline double FE_ir2(const Vec3d d, Vec3d& f, double cir2){
    double r2  = d.norm2();
    double ir2 = 1/r2;
    double e = cir2*ir2;
    //f.add_mul( dic, (eic)/r2 );
    f.set_mul( d, e/r2 );
    return e;
}

inline double FE_r2ir2(const Vec3d d, Vec3d& f, double cr2, double cir2){
    double r2  = d.norm2();
    double ir2 = 1/r2;
    double e = cir2*ir2 + cr2*r2;
    //f.add_mul( dic, (eic)/r2 );
    f.set_mul( d, e/r2 );
    return e;
}

inline double setDistance(Vec3d& p, const Vec3d p0, double r0){
    Vec3d  d; d.set_sub(p,p0);
    double r = d.norm();
    p.set_add_mul(p0,d,r0/r);
    return r;
}

struct EOFFAtom{
    EAFFelement* type = 0;
    //uint8_t  occup;
    Vec3d p;
    Vec3d f;
    double E;
    //Vec3d v;
    uint8_t  occups[N_ATOM_ORB];
    Vec3d    orb_p [N_ATOM_ORB];
    Vec3d    orb_f [N_ATOM_ORB];
    //Vec3d  orb_v[N_ATOM_ORB];
    //evalEnergyForce( );

    double intraEF(){
        double Ein= 0;
        //Vec3d orb_h[N_ATOM_ORB];
        for(int i=0; i<N_ATOM_ORB; i++){
            // NOTE: orbital position fixed on radius
            // TODO: or we can use forcefield for that ?
            //Vec3d fic;
            //E += FE_r2ir2( orb_p[i] - p, fic, type->oc_r2, type->oc_ir2);
            //f       .sub(fic);
            //orb_f[i].add(fic);

            // TODO: we may rather use radial forcefield instead of radial ?
            for(int j=0; j<N_ATOM_ORB; j++){
                Vec3d fij;
                Ein += FE_ir2( orb_p[j]-orb_p[i], fij, type->oo_ir2);
                orb_f[i].sub(fij);
                orb_f[j].add(fij);
            }
        }
        //E+=Ein;
        return Ein;
    }

    void moveGD(double dt){
        p.add_mul(f,dt); // move core
        for(int i=0; i<N_ATOM_ORB; i++){    // move orbitals
            orb_p[i].add_mul(orb_f[i],dt);
            setDistance(orb_p[i],p,type->oc_r0); // Fix orbitals on radius
        }
    }

};

class EOFF{ public:

    //EAFFelement* FF;
    std::vector<EOFFAtom> atoms;
    const int nNeighMax = 64;
    const int nOOMax    = 64;

    double R2cut = 1.0; // total max rcut
    double R2oo  = 1.0;

    void interEF_On2(){
        int natoms = atoms.size();
        int neighs[nNeighMax];
        int nneigh = 0;
        for(int i=0; i<natoms; i++){
            EOFFAtom& atomi = atoms[i];
            for(int j=i+1; j<natoms; j++){
                EOFFAtom& atomj = atoms[j];
                Vec3d d   = atomj.p - atomi.p;
                double r2 = d.norm2();
                if( r2>R2cut ) continue;
                neighs[nneigh]=j;
                nneigh++;
                // core-core energy
                //inter orbital energy
            }
            ooEF( atomi, nneigh, neighs );
        }
    }

    double ooEF( EOFFAtom& atomk, int nneigh, int* neighs ){

        int    oas[nOOMax];
        int    oos[nOOMax];
        double  rs[nOOMax];
        //double  hs[nOOMax];
        //Vec3d  ds [nOOMax];
        //Vec3d  ds [nOOMax];
        //double r2s[nOOMax];

        const double aik = 1.5;
        const double bik = 1.0;

        const double aijk = 1.5;
        const double bijk = 1.0;
        const double bij  = 1.0;

        double Eoo = 0.0;

        for(int ko=0; ko<N_ATOM_ORB; ko++){ // all orbitals of pivot atom
            Vec3d pk = atomk.orb_p[ko];
            int noo  = 0;
            for(int ii=0; ii<nneigh; ii++ ){
                int i = neighs[ii];
                EOFFAtom& atomi = atoms[i];
                for(int io=0; io<N_ATOM_ORB; io++){
                    Vec3d pi  = atomi.orb_p[io];
                    Vec3d dik = atomi.orb_p[io] - pk;
                    double r2 = dik.norm2();
                    if(r2>R2oo) continue;

                    double rik = sqrt(r2);
                    double eik = aik*exp( -bik*rik );

                    //dik.mul(1/rik);
                    Vec3d fik  = dik*(eik/rik);

                    Eoo += eik;
                    atomk.orb_f[ko].sub(fik);
                    atomi.orb_f[io].add(fik);

                    oos[noo] = io;
                    oas[noo] = i;
                    rs [noo] = rik;
                    //hs [noo] = dik*(1/rik);

                    //rs [i]   = r;
                    //hs [i].set_mul(dik,r);
                    noo++;

                    /*
                    for(int jj=0; jj<nneigh; jj++ ){
                        int j = neighs[jj];
                        EOFFAtom& atomj = atoms[j];
                        for(int jo=0; jo<N_ATOM_ORB; jo++){
                            Vec3d pj  = atomj.orb_p[io] - pi;
                            Vec3d dik = pk - pi;
                        }
                    }
                    */

                }
            }

            // exclusion of bonds
            for(int ioo=0; ioo<noo; ioo++ ){
                EOFFAtom& atomi = atoms[oas[ioo]];
                int io      = oos[ioo];
                Vec3d  pi   = atomi.orb_p[io];
                Vec3d  dik  = pi - pk;
                double rik  = rs[ioo];
                for(int joo=ioo; joo<noo; joo++ ){
                    EOFFAtom& atomj = atoms[oas[joo]];
                    int jo      = oos[joo];
                    Vec3d pj    = atomj.orb_p[jo];
                    Vec3d dij   = pj - pi;
                    double r2ij = dij.norm2();
                    if(r2ij>R2oo) continue;

                    //Vec3d  pij   = pj + pi;
                    //Vec3d  dijk  = pij - pk;
                    //double r2ijk = dijk.norm();
                    //double sij   = exp(bij*rij); // overlap ij
                    //double eijk  = sij * exp( *rijk );

                    Vec3d  djk  = pj - pk;
                    double rjk  = rs[joo];
                    double rij  = sqrt(r2ij);
                    double eijk = aijk * exp( bij*rij + bijk*(rik+rjk) );

                    Eoo        += eijk;

                    Vec3d fik = dik*(eijk/rik);
                    Vec3d fjk = djk*(eijk/rjk);
                    Vec3d fij = dij*(eijk/rij);

                    atomk.orb_f[ko].sub(fik);
                    atomk.orb_f[ko].sub(fjk);
                    atomi.orb_f[io].add(fik);
                    atomi.orb_f[io].sub(fij);
                    atomj.orb_f[jo].add(fik);
                    atomj.orb_f[jo].add(fij);
                }
                //atomk.orb_f[ko].sub(fik);
                //atomi.orb_f[io].add(fik);
            }
        } // ko
        return Eoo;
    }


};



#endif
