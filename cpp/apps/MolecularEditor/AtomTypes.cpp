
#include <cstdio>
#include <cstdlib>

#include "AtomTypes.h"  // THE HEADER

bool AtomTypes::loadFromFile( char const* filename ){
    printf(" loading AtomTypes from: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"r");
    fscanf (pFile, "%i", &ntypes);
    //printf("ntypes %i \n", ntypes );
    Zs       = new    int[ ntypes ];
    vdwRs    = new double[   ntypes ];
    vdwEs    = new double[   ntypes ];
    names    = new char* [   ntypes ];
    colors   = new uint32_t[   ntypes ];

    char hexstring[8];
    for (int i=0; i<ntypes; i++){
        names[i] = new char[6];
        fscanf (pFile, " %lf %lf %i %s %s", &vdwRs[i], &vdwEs[i], &Zs[i], names[i], hexstring );
        colors[i] = (uint32_t)strtol(hexstring, NULL, 16);
        printf( "%i %f %f %i %s  %s %i\n", i, vdwRs[i], vdwEs[i], Zs[i], names[i], hexstring, colors[i] );
    }
    fclose (pFile);
    return 0;
}

AtomTypes::AtomTypes( char const* filename ){  loadFromFile( filename );  };



