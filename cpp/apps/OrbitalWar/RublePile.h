
#ifndef  RublePile_h
#define  RublePile_h

#include <vector>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "geom3D.h"
#include "Body.h"

#include "Noise.h"
#include "SphereSampling.h"
#include "DrawSphereMap.h"

//#include "SpaceCraft.h"
//#include "EditSpaceCraft.h"


double powerlaw( double yrnd, double beta, double vmin, double vmax ){
    // http://mathworld.wolfram.com/RandomNumber.html
    // https://stackoverflow.com/questions/17882907/python-scipy-stats-powerlaw-negative-exponent/46065079#46065079
    double b    = beta+1;
    double pmin = pow(vmin,b);
    double pmax = pow(vmax,b);
    return pow( (pmax-pmin)*yrnd + pmin, 1./b );
    // return ( (k_max**(-gamma+1) - k_min**(-gamma+1))*y  + k_min**(-gamma+1.0))**1.0/(-gamma + 1.0)
}

double logerp( double x, double ymax, double ymin ){
    return ymin*pow( ymax/ymin, x );
}


class Boulder: public RigidBody{ public:

    int nsamp=0;
    int npix =0;
    float   hscale;
    float*  heights = 0;

    double radius;
    Vec3d  span;

    Boulder()=default;
    Boulder(int nsamp_, int nCrater_){ make( nsamp_, nCrater_ ); }
    ~Boulder(){ _dealloc(heights); }

    inline void setRadius(double radius_){
        radius = radius_;
        span.set(radius);
    }

    void realloc(int nsamp_){
        nsamp = nsamp_;
        npix  = 10*nsamp*nsamp;
        _realloc( heights, npix );
    }

    void make( int nsamp_, int nCrater_ ){
        realloc(nsamp_);
        makeSurf(nCrater_);
    }

    void makeSurf( int nCrater ){
        Vec3d*  craterPos = new Vec3d[nCrater];
        double* craterSz  = new double[nCrater];
        for(int i=0; i<nCrater; i++){
            craterPos[i].fromRandomSphereSample();
            craterSz[i] = randf()+0.4;
        }
        Vec2d dab = (Vec2d){ 1.0/nsamp, 1.0/nsamp };

        uint32_t quadrupole_seed = rand();

        for(int iface=0; iface<10; iface++ ){
            for(int ia=0; ia<nsamp; ia++ ){
                for(int ib=0; ib<nsamp; ib++ ){
                    Vec3d p;
                    SphereSampling::icosa2cartes( (Vec2i){nsamp,nsamp}, iface, ia*dab.a, ib*dab.b, p );
                    int i = iface*nsamp*nsamp + ia*nsamp + ib;
                    p.normalize();
                    double h=0;
                    //h  += Noise::getQuadruRand( p, quadrupole_seed )*2;
                    h += Noise::getCraterHeight( p, nCrater, 1.0, craterPos, craterSz )*0.3;
                    heights[i] = h; //+ sin( h*4 );
                }
            }
        }
        delete [] craterPos;
        delete [] craterSz;
    }

};

class RublePile{ public:
    std::vector<Boulder> boulders;
    SpaceBody* body = 0;
    double radius = 0;
    RublePile() = default;
    //RublePile(){}

    char* infoToBuff(char* sout){
        sout += sprintf( sout, "Name   %s \n", body->name.c_str()  );
        sout += sprintf( sout, "Radius %g  %g [m]   \n", radius, body->radius );
        sout += sprintf( sout, "Volume %g [m^3] \n", radius*radius*radius*4./3.*M_PI );
        sout += sprintf( sout, "Mass   %g [kg]  \n", body->mass   );
        sout += sprintf( sout, "Stype  `%s`  \n", body->Stype  );
        sout += sprintf( sout, " ## Orbit : \n" );
        sout += sprintf( sout, "semi_major  [AU]  %g \n", body->orbit->semi_major   );
        sout += sprintf( sout, "eccentricity  %g \n", body->orbit->eccentricity );
        double inc = acos( body->orbit->rot.c.z );
        sout += sprintf( sout, "inclination   %g rad (%g deg) \n", inc, inc*180.0/M_PI );
        return sout;
    };

    void makeRublePile( int nboulders, int nsamp, int nCrater, double mass, double mmin, double mmax, double beta ){
        double density = 2000.0; // [kg/m^3]
        std::vector<double> masses(nboulders);
        // distribution of masses is

        radius = pow( (mass/density)/(M_PI*4./3.), 1./3.);

        double msum = 0;
        for(int i=0; i<nboulders; i++){
            //double m  = powerlaw( randf(), beta, mmin, mmax );
            double m  = logerp( randf(), mmin, mmax );
            msum     += m;
            masses[i] = m;
        }
        double mscale = mass/msum;
        boulders.resize(nboulders);
        for(int i=0; i<nboulders; i++){
            double m = masses[i]*mscale;
            double V = m/density;
            double R = pow( V/(M_PI*4./3.), 1./3.);
            //printf( "Boulder[%i] mi %g[kg]  %g[kg] %g[m^3] %g[m] \n", masses[i], i, m, V, R );
            boulders[i].setRadius( R );
            boulders[i].hscale   = 0.25;
            boulders[i].make( nsamp, nCrater );
            boulders[i].pos.fromRandomCube( radius );
            boulders[i].rotMat=Mat3dIdentity;
        }
        //exit(0);
    }

    void relaxBoulders( double k, int niter ){
        for(Boulder& b : boulders){ b.vel=Vec3dZero; }
        for(int iter=0; iter<niter; iter++){
            for(Boulder& b : boulders){ b.force=Vec3dZero; }
            // get force
            for(int i=0; i<boulders.size(); i++){
                Vec3d  pi=boulders[i].pos;
                Vec3d& fi=boulders[i].force;
                double ri=boulders[i].radius;

                fi.add_mul( pi, -0.5 );

                for(int j=0; j<i; j++){
                    double rij = ri + boulders[j].radius;
                    Vec3d  d  = boulders[j].pos - pi;
                    double r  = d.norm();
                    double dr = r - rij;
                    //printf( "%i,%i dr %g |  rij %g  r %g \n", i,j, dr,  rij,  r );
                    if(dr<0){
                        //printf( "%i,%i dr %g \n", i,j, dr );
                        Vec3d f; f.set_mul( d, dr/r );
                        boulders[j].force.sub(f);
                        fi.add(f);
                    }
                }

            }
            for(Boulder& b : boulders){
                b.vel.mul    ( 1-k );
                b.vel.add_mul( b.force, k );
                b.pos.add    ( b.vel );
            }
        }
    }

};

#endif
