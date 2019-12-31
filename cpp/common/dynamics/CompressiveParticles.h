
#ifndef EFF_h
#define EFF_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"

#include "geom3D.h"

#include "InteractionsGauss.h"

const double const_Rgas    = 8.31446261815324; // J/(K*mol)
const double const_Vsphere = M_PI*4./3.;






double eosForce_IdealGas( double r, double U ){
    // pV=nRT=U
    // dU/dV = p = nRT/V
    // dU/dr = ( dU/dV )*( dV/dr ) = nRT/((4/3)*pi*r^3)   *  ( 4*pi*r^2 ) = nRT/(3*r)
    return U/(3*r);
}

/*
double compressAdiabat( double dR, double kappa, double& R, double& U ){
    // adiabatic process   https://en.wikipedia.org/wiki/Adiabatic_process
    // W = a*n*R*T1*( (V2/V1)^(1-k) - 1 ) = a*p1*V1*( (V2/V1)^(1-k) - 1 ) = U*( (V2/V1)^(1-k) - 1 )
    double V1 = R*R*R;
    R += dR;
    double V2 = R*R*R;
    U *= pow( V2/V1, 1-kappa );
}
*/

// state variable are Radius(R) and pressure(p)
double compressAdiabat( double dR, double kappa, double& R, double& p ){
    // p*V^kappa = const
    // p1*V1^k = p2*V2^k
    // p2 = p1 * (V1/V2)^k
    // p2 = p1 * (r1^3/r2^3)^k
    // p2 = p1 * (r1/r2)^(3*k)
    //double V1 = R*R*R;
    double oR=R;
    R += dR;
    //double V2 = R*R*R;
    double factor = pow( oR/R, 3*kappa );
    p *= factor;
    printf( "factor %g [1] p %g [GPa]\n",  factor,  p*1e-9  );
}


double overlap( const Vec3d& dR, double si, double sj, double pi, double pj, Vec3d& fp, double& fsi, double& fsj ){
    // E = pV = (p1+p2)*(V1+V2)  = p1*V1 + p2*V2   +    p1*V1 + p2*V2;

    double KK = 4; // 3 is normal

    double l = dR.norm();
    double E = 0;
    double dl = (si + sj) - l;
    if( dl>0 ){
        //double ci = dl/si;
        //double cj = dl/sj;
        //double Vi = si*si*si;
        //double Vj = sj*sj*sj;
        //double cVi = si*si*dl; // Vi*ci =  (si*si*si)*(dl/si)  =  si*si*dl  ;
        //double cVj = si*si*dl;
        //E = const_Vsphere*( pi*cVj + pj*cVi );
        double ksi = const_Vsphere*pj*(si*si);
        double ksj = const_Vsphere*pi*(sj*sj);
        //double K = const_Vsphere*( pi*sj2 + pj*si2 );
        double K=ksi+ksj;
        E = dl*K;
        fp.set_mul( dR, K/l );
        fsi =  KK*ksi; //  pj*const_Vsphere *   d(  si*si*( (si + sj) - l ) ) /dsi
        fsj =  KK*ksj;

        printf( "K %g pi,j %g %g si,j %g %g -> fsi,j %g %g \n", dl, pi, pj, si, sj, fsi, fsj );
        //printf( "dl %g ksij %g %g fp(%g,%g,%g) \n", dl, ksi, ksj, fp.x, fp.y, fp.z  );

        // = pushing one pressurized sphere into the other
    }else{
        fp  = Vec3dZero;
        fsi = 0; //  pj*const_Vsphere *   d(  si*si*( (si + sj) - l ) ) /dsi
        fsj = 0;
    }
    return E;
}


class CompressiveParticles{ public:

    double kappa = 1.6666; // https://en.wikipedia.org/wiki/Heat_capacity_ratio    monoatomic 5/3=1.666, diatomic 1.4
    double Kwall = 1.0e+11;

    int n;
    double* mass = 0;

    Vec3d*  pos  = 0;
    double* Rs   = 0;
    //double* Us  = 0; // internal energy
    double* press = 0;  // pressures

    Vec3d*  fpos = 0;
    double* fRs  = 0;

    Vec3d*  vpos = 0;
    double* vRs  = 0;

    int nwall;
    Plane3D* walls = 0;

    double Eone=0,Epair=0,Ewall=0,Ek=0;


void realloc(int n_, int nwall_){
    n=n_;
    _realloc( mass, n );

    _realloc( press, n );

    _realloc( pos, n  );
    _realloc( Rs,  n  );

    _realloc( fpos, n );
    _realloc( fRs,  n );

    _realloc( vpos, n );
    _realloc( vRs,  n );

    nwall=nwall_;
    _realloc( walls, nwall );
}

void dealloc(){
    delete [] mass;
    delete [] press;
    delete [] pos;
    delete [] Rs;
    delete [] fpos;
    delete [] fRs;
    delete [] vpos;
    delete [] vRs;
    //delete [] vradius;
}

/*
double evalOneBody(){
    Eone=0;
    for(int i=0; i<n; i++){
        fradius += eosForce_IdealGas( Rs[i], Us[i] );
        Eone = Us;
    }
    return Eone;
}
*/

double evalPairs(){
    for(int i=0; i<n; i++){
        Vec3d  posi   = pos[i];
        double Ri     = Rs [i];
        double pressi = press[i];
        for(int j=0; j<i; j++){
            //Vec3d  fp=Vec3dZero;
            Vec3d  fp;
            double fri,frj;
            Epair += overlap( pos[j]-posi, Ri, Rs[j], pressi, press[j], fp, fri, frj );
            fpos[i].sub(fp);
            fpos[j].add(fp);
            fRs [i]-=fri;
            fRs [j]-=frj;
        }
    }
    return Epair;
}

double evalWalls(){
    for(int i=0; i<n; i++){
        Vec3d  posi   = pos[i];
        double Ri     = Rs [i];
        double pressi = press[i];
        Vec3d  fp = Vec3dZero;
        double fr = 0;
        double Ki = Kwall * (Ri*Ri);
        for(int j=0; j<nwall; j++){
            const Plane3D& w = walls[j];
            //double dl  = -w.normal.dot(posi) - Ri - w.C;
            double dl  = -w.normal.dot(posi) + w.C + Ri;
            if(dl>0){
                double f = dl * Ki;
                fp.add_mul( w.normal, f );
                fr      -= f;
                //glColor3f( 0.0,1.0,0.0); Draw3D::drawVecInPos(  w.normal* f, posi );
                //glColor3f( 0.0,1.0,0.0); Draw3D::drawVecInPos(  w.normal, posi );
            }
        }
        fpos[i].add(fp);
        fRs [i]+=fr;
    }
    return Ewall;
}


double eval(){
    return
    //+ evalOneBody()
    + evalPairs()
    + evalWalls()
    ;
}

void clearVelocity(){
    for(int i=0; i<n; i++){ vpos[i]=Vec3dZero; vRs[i]=0; }
}

void clearForce(){
    for(int i=0; i<n; i++){ fpos[i]=Vec3dZero; fRs[i]=0; }
}

double moveMD(double dt){
    double cR    = 4.0;
    double rDamp = fmin( 0.5, dt/1e-6 );
    Ek = 0;
    for(int i=0; i<n; i++){
        double mi = mass[i];
        double dp = dt/mi;

        vpos[i].add_mul( fpos[i], dp );
        pos [i].add_mul( vpos[i], dt );
        Ek   += 0.5*mi*vpos[i].norm2();

        double fpress = press[i]*sq(Rs[i])*(4*M_PI); // pressure force   (dU/dr)

        vRs[i] *= (1-rDamp);
        vRs[i] += ( fRs[i] + fpress )*(dp*cR);
        printf( " fRs[%i] %g[GPs] fpress %g[GPa] vR %g[km/s] \n", i, fRs[i]*1e-9, fpress*1e-9, vRs[i]*1e-3 );
        //vRs[i] +=  fRs[i]*(dp*cR);
        //compressAdiabat( vRs[i]*dt, kappa, Rs[i], Us[i] );
        compressAdiabat( vRs[i]*dt, kappa, Rs[i], press[i] );
        //Rs  [i] += vRs[i]*dt;
        //Rs  [i] += fRs[i]*dt;
        //Us  [i] += f;

    }
    return Ek;
}

void setMoles(int i, double nMols, double molarMass, double temperature, double radius){
    // Ideal gas : pV = nRT   =>    p = (n/V)RT = (rho/M)RT ~ rho^1
    // Fermi gas : p  ~ rho^(5/3)
    double m = nMols * molarMass;
    double U = nMols * const_Rgas * temperature;
    double V = const_Vsphere * radius*radius*radius;
    mass [i] = m;
    Rs   [i] = radius;
    press[i] = U/V;
}

void setDensity(int i, double density, double molarMass, double temperature, double radius ){
    // Ideal gas : pV = nRT   =>    p = (n/V)RT = (rho/M)RT ~ rho^1
    // Fermi gas : p  ~ rho^(5/3)
    //double m = nMols * molarMass;
    double V = const_Vsphere * radius*radius*radius;
    mass [i] = density * V;
    Rs   [i] = radius;
    press[i] = const_Rgas*temperature * (density/molarMass);

    printf( "[%i] rho %g [kg/l] M %g [g/mol] T %g [K] r %g[m]\n", i, density*1e-3, molarMass*1e+3, temperature, radius );
    printf( " ->  mass %g [kg] V %g [l] U %g [MJ] p %g [GPa] \n",    mass[i], V*1e+3, press[i]*V*1e-6, press[i]*1e-9 );
}

};

#endif
