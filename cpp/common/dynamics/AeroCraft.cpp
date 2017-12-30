
#include <math.h>
#include <cstdlib>
#include <stdio.h>

//#include <SDL2/SDL_opengl.h>

#include "AeroCraft.h" // THE HEADER

void AeroCraft::applyAeroForces( const Vec3d& vwind ){
    Vec3d vair = vwind - vel;
    //panels[0].DEBUGsurf = true;
    for( int i=0; i<nPanels; i++ ){
        //printf( " %i %i \n", i, nPanels );
        //panels[i].applyForceSimple( vair );
        panels[i].applyForce( vair );
    }
    //Mat3 rmat; rmat.setT(rotMat);
    totalThrust.set(0.0d);
    for( int i=0; i<nPropelers; i++ ){
        if( propelers[i].power > 0 ){
            Vec3d gdpos,gdir;
            rotMat.dot_to( propelers[i].dir,  gdir  );
            rotMat.dot_to( propelers[i].lpos, gdpos );
            double vdir   = -gdir.dot( vair );  // printf("vdir %g \n", vdir);
            double thrust = propelers[i].getThrust(vdir);
            //glColor3d(1.0f,1.0f,0.0f); Draw3D::drawVecInPos( gdir, gdpos+pos );  Draw3D::drawPointCross( gdpos+pos, 0.5 );
            gdir.mul(thrust);
            //glColor3d(1.0f,1.0f,0.0f); Draw3D::drawVecInPos( gdir, gdpos+pos );  Draw3D::drawPointCross( gdpos+pos, 0.2 );
            totalThrust.add(gdir);
            apply_force( gdir, gdpos );
        }
    }
}

int AeroCraft::fromFile( const char * fname ){
    const int nbuf = 1024;
    char buf  [nbuf];
    FILE * pFile;
    pFile = fopen (fname,"r");
    printf(" loading molTypes from: >>%s<<\n", fname );

    fscanf (pFile, " %lf %lf %lf %lf\n", &mass, &Ibody.xx, &Ibody.yy, &Ibody.zz );
    printf( " %lf %lf %lf %lf\n", mass, Ibody.xx, Ibody.yy, Ibody.zz );

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

    Ibody.invert_to( invIbody );
    qrot.setOne();
    qrot.toMatrix(rotMat);
    L.set(0,0,0);
    setMass( mass );
    vel.set(0,0,0);
    pos.set(0,0,0);
    clean_temp();

    printf("AeroCraft loaded\n");
    fclose(pFile);
    return nPanels;
}


