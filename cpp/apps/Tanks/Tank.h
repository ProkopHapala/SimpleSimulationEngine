
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

class Wheel{ public:
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

class ArmorPlate{ public:
    double thickness;
    double mass;
    Vec3d  normal;
    int material;
};

class VehicleBlock : public KinematicBody, public Mesh { public:
    std::vector<ArmorPlate> armor;

    double getMaxArmor(){
        double maxThickness = 0.0;
        for(ArmorPlate& pl : armor){
            maxThickness = fmax(pl.thickness, maxThickness);
        }
        return maxThickness;
    }

    double getArmorMass( double density ){
        double m = 0.0;
        for(int i=0; i<polygons.size(); i++){
            double D  = armor[i].thickness*1e-3;
            double S  = polygonArea( i, &armor[i].normal );
            double mi = D*S*density;
            armor[i].mass = mi;
            m+=mi;
            printf( "%i D %f S %f dm %f n (%g,%g,%g)%g \n", i, D, S, D*S*density, armor[i].normal.x, armor[i].normal.y, armor[i].normal.z, armor[i].normal.norm() );
        }
        return m;
    }

    int loadArmor( char * fname ){
        printf("VehicleBlock.loadArmor %s\n", fname );
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){ printf("cannot find %s\n", fname ); return -1; }
        char buff [1024];
        char mat_name[64];
        char * line;
        int nl;
        for( int i=0; i<polygons.size(); i++ ){
            double thickness;
            line = fgets( buff, 1024, pFile );
            sscanf( line, "%s %lf", mat_name, &thickness );
            printf( "%i %s %lf", mat_name, &thickness );
            armor.push_back((ArmorPlate){thickness});
        }
    }

    void fromString( char * str ){
        printf("VehicleBlock.fromString %s\n", str );
        //char* fname1,fname2;
        char fname1[64];
        char fname2[64];
        sscanf( str, "%s %s\n", fname1, fname2 );
        fromFileOBJ(fname1);
        loadArmor(fname2);
    }
};

class Tank : public Warrior3D { public:
    int nwheel = 0;
    Wheel * wheels;

    VehicleBlock hull;
    VehicleBlock turret;

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

    int fromFile( char * fname ){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){ printf("cannot find %s\n", fname ); return -1; }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); hull  .fromString( line );
        line = fgets( buff, 1024, pFile ); turret.fromString( line );
    }

};

#endif
