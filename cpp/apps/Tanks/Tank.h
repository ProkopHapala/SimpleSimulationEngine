
#ifndef Tank_h
#define Tank_h

#include "appliedPhysics.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Object3D.h"
#include "Terrain25D.h"
#include "Warrior3D.h"
#include "Projectile3D.h"
#include "Turret.h"

class Wheel{ public:
    double k,mass,damp;
    Vec3d  lpos;
    double h,vh,fh;

    inline void   clean_temp(){ fh=0.0; }
    inline double getSpringForce(){return k*(lpos.y-h); }
    inline void move( double dt){
        vh     += ( vh*damp + fh + getSpringForce() )*dt;
        lpos.z += vh*dt;
    }
    inline void init( Vec3d lpos_, double k_, double mass_, double damp_ ){
        k=k_; mass=mass_; damp=damp_; lpos=lpos_; h=lpos.y; vh=0.0; fh=0.0;
    }
};

class ArmorPlate{ public:
    double thickness;
    double mass;
    //Vec3d  normal;
    int material;
};

class VehicleBlock : public KinematicBody, public Mesh { public:
    std::vector<ArmorPlate> armor;
    int glo_captions;
    int glo_armor;

    double ray( const Vec3d &ray0_, const Vec3d &hRay_, Vec3d& normal, int& imin ){
        //ray0 -= lpos;
        Vec3d hRay = lrot.dot(hRay_);
        Vec3d ray0 = lrot.dot(ray0_);
        return Mesh::ray( ray0, hRay, normal, imin );
        //normal     = lrot.dotT(normal);
    }

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
            //double S  = polygonArea( i, &armor[i].normal );
            double S  = polygonArea( i, NULL );
            double mi = D*S*density;
            armor[i].mass = mi;
            m+=mi;
            //printf( "%i D %f S %f dm %f n (%g,%g,%g)%g \n", i, D, S, D*S*density, armor[i].normal.x, armor[i].normal.y, armor[i].normal.z, armor[i].normal.norm() );
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

class Tank : public RigidBody, public Warrior3D { public:
    int nwheel = 0;
    Wheel        * wheels;
    VehicleBlock   hull;
    VehicleBlock   turret;

    double elevation = 0.0;

    double maxPower= 1.0;

    double power_gear[2] = {0.0,0.0};

	void makeWheels( int n, double xmin, double xmax, double width, double height, double k, double mass, double damp ){
        nwheel = n*2;
        wheels = new Wheel[nwheel];
        double dx = (xmax-xmin)/(nwheel-1);
        //printf( "makeWheels n %i xmin %g xmax %g width %g height %g k %g m %g damp %g \n", n, xmin, xmax, width, height, k, m, damp );
        for(int i=0;i<n;i++){
            int i2 = i<<1;
            wheels[i2  ].init( {xmin+dx*i2, height, width*+0.5}, k, mass, damp );
            wheels[i2+1].init( {xmin+dx*i2, height, width*-0.5}, k, mass, damp );
            //printf( "wheel %i (%g,%g,%g)\n", i2  , wheels[i2  ].lpos.x, wheels[i2  ].lpos.y, wheels[i2  ].lpos.z   );
            //printf( "wheel %i (%g,%g,%g)\n", i2+1, wheels[i2+1].lpos.x, wheels[i2+1].lpos.y, wheels[i2+1].lpos.z   );
        }
        //exit(0);
	};

	void rotateTurret( double angle ){
		double ca   = cos(angle);
		double sa   = sin(angle);
        turret.lrot.c.rotate_csa( ca, sa, turret.lrot.b );
        turret.lrot.a.rotate_csa( ca, sa, turret.lrot.b );
	}

    void rotateTurretToward( const Vec3d& dir ){
        Mat3d grot;
        turret.globalRotT(rotMat, grot);

        double sa,ca,sa_;

        double sa0=gun_rot.dot( grot.b  );
        //sa_  =  dir.dot( grot.b ) - sa0;
        //if( sa_>0 ){ sa=sa0+0.01; }else{sa=sa0-0.01; };
        sa = sa0 + clamp( dir.dot( grot.b )-sa0, -0.01, 0.01 );
        ca = sqrt( 1 - sa*sa );
        gun_rot = grot.b*sa + grot.a*ca;

        sa = clamp( -dir.dot( grot.c ), -0.01, 0.01 );
        ca = sqrt(1 - sa*sa);
        turret.lrot.c.rotate_csa( ca, sa, turret.lrot.b );
        turret.lrot.a.rotate_csa( ca, sa, turret.lrot.b );

	}

	double ray( const Vec3d &ray0_, const Vec3d &hRay_, int& ipl, VehicleBlock*& block, double& effthick, Vec3d& normal ){
        //ray0.sub(pos);
        Vec3d hRay = rotMat.dot( hRay_ );
        Vec3d ray0 = rotMat.dot( ray0_-pos );
        //Vec3d hRay = rotMat.dotT( hRay_ );
        //Vec3d ray0 = rotMat.dotT( ray0_-pos );

        Vec3d normal2;
	    int itr,itr2;
	    ipl=-1;
        block = &hull;

        double t  = hull  .ray( ray0, hRay, normal,  itr  );
        double t2 = turret.ray( ray0, hRay, normal2, itr2 );
        if( (itr2>=0)&&(t2<t) ){
            itr=itr2;
            block = &turret;
            t=t2;
            normal=normal2;
            //normal     = lrot.dotT(normal);
        }
        normal.normalize();
        normal = rotMat.dot(normal);
        if( itr>=0 ){
            ipl = block->tri2poly[itr];
            double thick    = block->armor[ipl].thickness;
            //double cdot     = block->armor[ipl].normal.dot( hRay );
            normal = block->lrot.dotT(normal);
            double cdot     = normal.dot( hRay_ );
            effthick = thick/-cdot;
        }
        return t;
	}

	inline void getWheelPos( int i, Vec3d& gpos ){
        rotMat.dot_to_T( {wheels[i].lpos.x,wheels[i].h,wheels[i].lpos.z}, gpos );
        gpos.add(pos);
	}

    void interactTerrain( Terrain25D * terrain ){
        //Mat3d mrot; qrot.toMatrix_T(mrot);
        //Mat3d mrot; qrot.toMatrix_T(mrot);
        Mat3d mrot = rotMat;

        double wheelLock[2]={0.0,0.0};
        wheelLock[0] = 1-clamp(fabs(power_gear[0]), 0.0, 0.7);
        wheelLock[1] = 1-clamp(fabs(power_gear[1]), 0.0, 0.7);
        //printf( "wheelLock %g %g\n", wheelLock[0], wheelLock[1], power_gear[0], power_gear[1] );
        for(int i=0; i<nwheel; i++){
            Vec3d gdp,gv,gp;
            velOfPoint(wheels[i].lpos, gv, gdp);
            gp.set_add(gdp,pos);
            Vec2d dv;
            double h  = terrain->eval( {gp.x,gp.z}, dv);
            double dh = h-gp.y;
            if(dh>0){
                Vec3d f;
                double fnormal = dh*wheels[i].k; if( gv.y>0 ) fnormal*=wheels[i].damp;
                f.set( -dv.x*fnormal, fnormal, -dv.y*fnormal );
                f.add_mul( mrot.c, -5.0*gv.dot(mrot.c) );
                f.add_mul( mrot.a,  5.0*(power_gear[i&1]*maxPower - gv.dot(mrot.a)*wheelLock[i&1] ) );
                apply_force(f,gdp);
            }
        }
    }

    //virtual void interact( Terrain25D * terrain ){ interactTerrain( terrain ); }

    //virtual void asRigidBody(){}

    virtual RigidBody* asRigidBody( ){ return static_cast<RigidBody*>(this); };

    virtual void move_warrior( double dt, Vec3d& wind_speed, Vec3d& gravity, Terrain25D * terrain ){
        clean_temp();
        force.add( gravity*mass );
        if(terrain){ interactTerrain( terrain ); }
        //w->landed = collideWithWorld ( w->pos, w->vel, w->surf );
        move( dt );
    };

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
