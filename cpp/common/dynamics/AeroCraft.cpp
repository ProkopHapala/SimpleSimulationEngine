
#include <math.h>
#include <cstdlib>
#include <stdio.h>

//#include <SDL2/SDL_opengl.h>
#include "globals.h"

#include "AeroCraft.h" // THE HEADER

void AeroCraft::applyAeroForces( const Vec3d& vwind ){
    Vec3d vair = vwind - vel;
    //panels[0].DEBUGsurf = true;
    for( int i=0; i<nPanels; i++ ){
        //printf( " %i %i \n", i, nPanels );
        //panels[i].applyForceSimple( vair );
        panels[i].applyForce( vair );
        //printf(  "panel %i force (%g,%g,%g) torq (%g,%g,%g) \n", i,  force.x,force.y,force.z,   torq.x,torq.y,torq.z );
    }

    //Mat3 rmat; rmat.setT(rotMat);
    totalThrust.set(0.0d);
    for( int i=0; i<nPropelers; i++ ){
        if( propelers[i].power > 0 ){
            Vec3d gdpos,gdir;
            rotMat.dot_to_T( propelers[i].dir,  gdir  );
            //printf( "prop.fw %g \n", rotMat.c.dot(gdir)  );
            rotMat.dot_to_T( propelers[i].lpos, gdpos );
            double vdir   = -gdir.dot( vair );  // printf("vdir %g \n", vdir);
            double thrust = propelers[i].getThrust(vdir);
            //glColor3d(1.0f,1.0f,0.0f); Draw3D::drawVecInPos( gdir, gdpos+pos );  Draw3D::drawPointCross( gdpos+pos, 0.5 );
            gdir.mul(thrust);
            //glColor3d(1.0f,1.0f,0.0f); Draw3D::drawVecInPos( gdir, gdpos+pos );  Draw3D::drawPointCross( gdpos+pos, 0.2 );
            totalThrust.add(gdir);
            apply_force( gdir, gdpos );
            //printf( " thrust.fw %g |thrust| %g %g \n", rotMat.c.dot(gdir), gdir.norm(), thrust  );
        }
    }

}

double AeroCraft::getTotalPower() const{
    double power=0;
    for(int i=0; i<nPropelers; i++){ power+=propelers[i].power; }
    return power;
};

int AeroCraft::fromFile( const char * fname ){
    printf(" AeroCraft::fromFile: >>%s<<\n", fname );
    const int nbuf = 1024;
    char buf  [nbuf];
    FILE * pFile=0;
    pFile = fopen (fname,"r");
    if(pFile==0){ printf("ERROR AeroCraft::fromFile(%s): cannot open file\n", fname ); exit(0); }
    Vec3d Ispan;
    fscanf (pFile, " %lf %lf %lf %lf\n", &mass, &Ispan.x, &Ispan.y, &Ispan.z );
    printf(        " %lf %lf %lf %lf\n",  mass,  Ispan.x,  Ispan.y,  Ispan.z );
    fscanf ( pFile, "%i\n", &nPanels);
    panels = new AeroSurface[nPanels];
    for(int i=0; i<nPanels; i++){
        fgets( buf, nbuf, pFile);
        //AeroSurface * a = new AeroSurface();
        //panels[i].fromString( buf );
        panels[i].fromStringPolarModel( buf );
        panels[i].craft = this;
    }
    int ileftAirelon,irightAirelon,ielevator,irudder;
    fscanf ( pFile, "%i %i %i %i\n", &ileftAirelon,&irightAirelon,&ielevator,&irudder);
    printf( "%i %i %i %i\n", ileftAirelon,irightAirelon,ielevator,irudder );
    leftAirelon =&panels[ileftAirelon -1];
    rightAirelon=&panels[irightAirelon-1];
    elevator    =&panels[ielevator    -1];
    rudder      =&panels[irudder      -1];

    fscanf ( pFile, "%i\n", &nPropelers);
    propelers = new Propeler[nPropelers];
    for(int i=0; i<nPropelers; i++){
        fgets( buf, nbuf, pFile);
        propelers[i].fromString( buf );
        //propelers[i].craft = this;
    }

    //Ibody.invert_to( invIbody );
    //qrot.setOne();
    //qrot.toMatrix(rotMat);
    L.set(0,0,0);
    //setMass( mass );
    setInertia_box(mass, Ispan);
    vel.set(0,0,0);
    pos.set(0,0,0);
    clean_temp();
    //update_aux(); // MUST BE CALLED BEFORE SIMULATION STARTS !!!

    printf("AeroCraft loaded\n");
    fclose(pFile);
    return nPanels;
}


