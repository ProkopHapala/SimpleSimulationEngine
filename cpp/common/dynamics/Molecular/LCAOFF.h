
#ifndef LCAOFF_h
#define LCAOFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
//#include "GridFF.h"

struct TBparams{
    // decay
    double bss;
    double bsp;
    double bpp;
    double bp2;
    // amplitude
    double ass;
    double asp;
    double app;
    double ap2;
}

struct LCAOatom{
    int    type;
    Vec3d  pos;
    Vec3d  fpos;
    Quat4d  orbs[4]; // orbital coefficient
    Quat4d forbs[4]; // orbital coefficient foces

    void orthoNormForces(double knorm, double korth){
        for(int i=0; i<4; i++){
            Quat4d& forbi      = forbs[i];
            const Quat4d& orbi =  orbs[i];
            double cii = orbi.norm2();
            forbi.add_mul( orbi, (1-cii)*knorm );
            for(int j=i+1; j<4; j++){
                Quat4d fij;
                fij.set_mul( orbi,  korth*orbi.dot( orbs[j] ));
                forbi   .sub( fij );
                forbs[j].add( fij );
            }
        }
    }

    void moveOrbsGD(double dt){ for(int i=0; i<4; i++){ orbs[i].add_mul(forbs[i], dt); } }
    void move      (double dt){ pos.add_mul(fpos,dt) }

}

// https://chem.libretexts.org/Textbook_Maps/Physical_and_Theoretical_Chemistry_Textbook_Maps/Map%3A_Physical_Chemistry_(McQuarrie_and_Simon)/09%3A_The_Chemical_Bond%3A_Diatomic_Molecules/09.3%3A_The_Overlap_Integral
inline double radial(double r){ return exp(r)*(1+r+0.33333*r*r); }

inline void dradial(double r, double amp, double& f, double& df ){ 
    double expar = amp*exp(r);
    f  = expar*(1 +   r + 0.33333333*r*r ); 
    df = expar*(6 + 5*r +            r*r ); 
}

void bondEnergy( LCAOatom& a, LCAOatom& b, const TBparams& par ){

    //Quat4d S[4];  // overlap matrix
    //double S[16];

    // ==== Build overlap matrix S

    Vec3d  h = b.pos - a.pos;
    double r = h.normalize();

    double Rss = radial(r*par.bss)*ass;
    double Rsp = radial(r*par.bsp)*asp;
    double Rpp = radial(r*par.bpp)*app;
    double Rp2 = radial(r*par.bp2)*ap2;

    /*
    S[0 ] = Rp2 * (3*h.x*h.x - 1);
    S[1 ] = Rpp * (3*h.x*h.y);
    S[2 ] = Rpp * (3*h.x*h.z);
    S[3 ] = Rsp * h.x;

    S[4 ] = Rp2 * (3*h.y*h.x - 1);
    S[5 ] = Rpp * (3*h.y*h.y);
    S[6 ] = Rpp * (3*h.y*h.z);
    S[7 ] = Rsp * h.y;

    S[8 ] = Rp2 * (3*h.z*h.x - 1);
    S[9 ] = Rpp * (3*h.z*h.y);
    S[10] = Rpp * (3*h.z*h.z);
    S[11] = Rsp * h.z;

    S[12] = Rsp * (3*h.x*h.x - 1);
    S[13] = Rsp * (3*h.x*h.y);
    S[14] = Rsp * (3*h.x*h.z);
    S[15] = Rss;
    */

    double Sxx   = Rp2 * (3*h.x*h.x - 1);
    double Sxy   = Rpp * (3*h.x*h.y    );
    double Sxz   = Rpp * (3*h.x*h.z    );
    //double Sxs = Rsp * h.x;

    //double Syx = Rp2 * (3*h.y*h.x  );
    double Syy   = Rpp * (3*h.y*h.y - 1);
    double Syz   = Rpp * (3*h.y*h.z    );
    //double Sxs = Rsp * h.x;

    //double Szx = Rp2 * (3*h.z*h.x);
    //double Szy = Rpp * (3*h.z*h.y);
    double Szz   = Rpp * (3*h.z*h.z - 1);
    //double Szs = Rsp * h.x;

    double Ssx = Rsp * h.x;
    double Ssy = Rsp * h.y;
    double Ssz = Rsp * h.z;
    double Sss = Rss;

    // ==== Eval energy

    Quat4d asum = Quat4dZero;
    asum.add(a.orbs[0]); 
    asum.add(a.orbs[1]); 
    asum.add(a.orbs[2]); 
    asum.add(a.orbs[3]);

    Quat4d bsum = Quat4dZero;
    bsum.add(b.orbs[0]); 
    bsum.add(b.orbs[1]); 
    bsum.add(b.orbs[2]); 
    bsum.add(b.orbs[3]);

    double E = 
    + Sss*( asum.s*bsum.s )
    + Sxx*( asum.x*bsum.x )
    + Syy*( asum.y*bsum.y )
    + Szz*( asum.z*bsum.z )

    + Ssx*( asum.s*bsum.x + asum.x*bsum.s ) 
    + Ssy*( asum.s*bsum.y + asum.y*bsum.s ) 
    + Ssz*( asum.s*bsum.z + asum.z*bsum.s ) 

    + Sxy*( asum.x*bsum.y + asum.y*bsum.x ) 
    + Sxz*( asum.x*bsum.z + asum.z*bsum.x ) 
    + Syz*( asum.y*bsum.z + asum.z*bsum.y ) ;

}

//double bondEF( LCAOatom& a, LCAOatom& b,  const TBparams& par, Vec3d& force, double* fa, double* fa ){
double bondEF( LCAOatom& a, LCAOatom& b, const TBparams& par){

    //Quat4d S[4];  // overlap matrix
    //double S[16];

    // ==== Build overlap matrix S

    Vec3d  h = b.pos - a.pos;
    double r = h.normalize();

    double Rss,Rsp,Rpp,Rp2;
    double dRss,dRsp,dRpp,dRp2;

    dradial(r*par.bss,par.ass,Rss,dRss);
    dradial(r*par.bsp,par.asp,Rsp,dRsp);
    dradial(r*par.bpp,par.app,Rpp,dRpp);
    dradial(r*par.bp2,par.ap2,Rp2,dRp2);

    double hxx = sq(h.x);
    double hyy = sq(h.y);
    double hzz = sq(h.z);
    double hxy = h.x*h.y;
    double hxz = h.x*h.z;
    double hyz = h.y*h.z;
    double hxyz = hxy * h.z;

    double Ssx = Rsp * h.x;
    double Ssy = Rsp * h.y;
    double Ssz = Rsp * h.z;
    double Sss = Rss;

    double Sxx = Rp2 * (3*hxx - 1);
    double Sxy = Rpp * (3*hxy    );
    double Sxz = Rpp * (3*hxz    );

    double Syy = Rpp * (3*hyy - 1);
    double Syz = Rpp * (3*hyz    );

    double Szz = Rpp * (3*hzz - 1);

    /*
    double Ssx_x = Rsp * (hyy+hzz)*ir +   dRsp * hxx ;
    double Ssx_y = Rsp * (hxy    )*ir +   dRsp * hyx ;
    double Ssx_z = Rsp * (hxz    )*ir +   dRsp * hzx ;

    double Sxx_x = Rp2 * (hyy+hzz)*h.x*ir*6   +   dRsp * h.x * (3*hxx - 1);
    double Sxx_y = Rp2 * (hxx    )*h.y*ir*6   +   dRsp * h.y * (3*hxx - 1);
    double Sxx_z = Rp2 * (hxx    )*h.z*ir*6   +   dRsp * h.z * (3*hxx - 1);

    double Sxy_x = Rp2 * (1-2*hxx)*hy*ir*3   +   dRsp * h.x * 3*hxy;
    double Sxy_y = Rp2 * (1-2*hyy)*hx*ir*3   +   dRsp * h.y * 3*hxy;
    double Sxy_z = Rp2 *  hxyz       *ir*6   +   dRsp * h.z * 3*hxy;
    */

    // ==== Eval energy
    // TODO : relaxation of LCAO coefs should do several sub-steps within one nucleai step

    Quat4d asum = Quat4dZero; for(int i=0; i<4; i++){ asum.add(a.orbs[i]); }
    Quat4d bsum = Quat4dZero; for(int i=0; i<4; i++){ bsum.add(b.orbs[i]); }

    double E = 
    + Sss*( asum.s*bsum.w )
    + Sxx*( asum.x*bsum.x )
    + Syy*( asum.y*bsum.y )
    + Szz*( asum.z*bsum.z )

    + Ssx*( asum.s*bsum.x + asum.x*bsum.w ) 
    + Ssy*( asum.s*bsum.y + asum.y*bsum.w ) 
    + Ssz*( asum.s*bsum.z + asum.z*bsum.w ) 

    + Sxy*( asum.x*bsum.y + asum.y*bsum.x ) 
    + Sxz*( asum.x*bsum.z + asum.z*bsum.x ) 
    + Syz*( asum.y*bsum.z + asum.z*bsum.y ) ;

    // forces on atomic position

    // forces on LCAO expansion coefs
    for(int i=0; i<4; i++){
        a.forbs[i].x += a.forbs[i].x * ( Ssx*bsum.w + Sxx*bsum.x + Sxy*bsum.y  + Sxz*bsum.z );
        a.forbs[i].y += a.forbs[i].y * ( Ssy*bsum.w + Sxy*bsum.x + Syy*bsum.y  + Syz*bsum.z );
        a.forbs[i].z += a.forbs[i].z * ( Ssz*bsum.w + Sxz*bsum.x + Syz*bsum.y  + Szz*bsum.z );
        a.forbs[i].w += a.forbs[i].w * ( Sss*bsum.w + Ssx*bsum.x + Ssy*bsum.y  + Ssz*bsum.z );

        b.forbs[i].x += b.forbs[i].x * ( Ssx*asum.w + Sxx*asum.x + Sxy*asum.y  + Sxz*asum.z );
        b.forbs[i].y += b.forbs[i].y * ( Ssy*asum.w + Sxy*asum.x + Syy*asum.y  + Syz*asum.z );
        b.forbs[i].z += b.forbs[i].z * ( Ssz*asum.w + Sxz*asum.x + Syz*asum.y  + Szz*asum.z );
        b.forbs[i].w += b.forbs[i].w * ( Sss*asum.w + Ssx*asum.x + Ssy*asum.y  + Ssz*asum.z );
    }

    return E;

}











class LCAOFF{ public:
    LCAOatom* atoms = 0;

};



#endif
