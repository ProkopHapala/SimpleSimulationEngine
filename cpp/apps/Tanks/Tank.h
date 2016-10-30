
#ifndef Tank_h
#define Tank_h

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Object3D.h"
#include "Terrain25D.h"
#include "Warrior3D.h"
#include "Projectile3D.h"

class Wheel{
    public:
    double k,m,damp;
    Vec3d  lpos;
    double h,vh,fh;

    inline void   clean_temp(){ fh=0.0; }
    inline double getSpringForce(){return k*(lpos.y-h); }
    inline void move( double dt){
        vh     += ( vh*damp + fh + getSpringForce() )*dt;
        lpos.z += vh*dt;
    }
    inline void init( Vec3d lpos_, double k_, double m_, double damp_ ){
        k=k_; m=m_; damp=damp_; lpos.set( lpos); h=lpos.y; vh=0.0; fh=0.0;
    }
};

class Tank : public Warrior3D {
	public:
    int nwheel = 0;
    Wheel * wheels;

	void makeWheels( int n, double xmin, double xmax, double width, double height, double k, double m, double damp ){
        wheels = new Wheel[n*2];
        double dx = (xmax-xmin)/(n-1);
        for(int i=0;i<n;i++){
            int i2 = i<<1;
            wheels[i2  ].init( {dx*i2, height, width*+0.5}, k, m, damp );
            wheels[i2+1].init( {dx*i2, height, width*-0.5}, k, m, damp );
        }
	};

    void interactTerrain( Terrain25D * terrain ){
        Vec3d gpos,normal;
        Vec2d dv;
        normal.set(0.0,1.0,0.0);
        for(int i=0; i<nwheel; i++){
            wheels[i].clean_temp();
            rotMat.dot_to( {wheels[i].lpos.x,wheels[i].h,wheels[i].lpos.z}, gpos );
            double h  = terrain->eval( {gpos.x,gpos.z}, dv );
            double dh = gpos.y - h;
                if( dh  < 0.0 ){
                    //w->force.add( {0.0,gravity,0.0} );
                    //w->force.add( {dv.x, dh*(-1-0.5*w->vel.y), dv.y} );
                    //wheels[i].fh += rotMat.b.dot( normal ) * dh;
                    wheels[i].fh += dh;
                    //force.add( {vel.x*landDrag,0.0,vel.z*landDrag} );
                }
            apply_force( rotMat.b*wheels[i].getSpringForce(), gpos-pos );
        }
    }

    void update( double dt ){
        for(int i=0; i<nwheel; i++){ wheels[i].move( dt); }
        move(dt);
    }

};

#endif
