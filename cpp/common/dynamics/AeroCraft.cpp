
#include <math.h>
#include <cstdlib>
#include <stdio.h>

//#include <SDL2/SDL_opengl.h>

#include "AeroCraft.h" // THE HEADER

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
    pos.set(0,500,0);
    clean_temp();

    printf("AeroCraft loaded\n");
    fclose(pFile);
    return nPanels;
}


