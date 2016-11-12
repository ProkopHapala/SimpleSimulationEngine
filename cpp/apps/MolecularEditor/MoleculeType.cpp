
#include <cstdio>
#include <cstdlib>

#include "MoleculeType.h" // THE HEADER

void MoleculeType::allocateAtoms( int n ){
    atypes = new    int[ natoms ];
    xyzs   = new  Vec3d[ natoms ];
    Qs     = new double[ natoms ];
}

MoleculeType::MoleculeType( int natoms_ ){
    natoms  = natoms_;
    allocateAtoms( natoms );
}

// --- function implementation

bool MoleculeType::loadFromFile_bas( char const* filename ){
    printf(" filename: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"r");
    fscanf (pFile, " %i\n", &natoms);
    printf("natoms %i \n", natoms );
    allocateAtoms( natoms );
    char buf[256];
    for (int i=0; i<natoms; i++){
        double q=0;
        //int ai; double x,y,z;
        //int nw = fscanf (pFile, " %i %lf %lf %lf %lf", &atypes[i], &xyzs[i].x, &xyzs[i].y, &xyzs[i].z, &q );
        //int nw = fscanf (pFile, " %i %lf %lf %lf %lf\n", &ai, &x, &y, &z, &q );   atypes[i]=ai; xyzs[i].x=x; xyzs[i].y=y;  xyzs[i].z=z;
        fgets( buf, 256, pFile); //printf( ">%s<\n", buf );
        int nw = sscanf (buf, " %i %lf %lf %lf %lf", &atypes[i], &xyzs[i].x, &xyzs[i].y, &xyzs[i].z, &q );
        //int nw = sscanf (buf, " %i %lf %lf %lf %lf", &ai, &x, &y, &z, &q );   atypes[i]=ai; xyzs[i].x=x; xyzs[i].y=y;  xyzs[i].z=z;

        if( nw > 4 ){  Qs[i] = q; }else{ Qs[i] = 0; }
        atypes[i]--;
        //printf( " %i %f %f %f %f    %f %i \n", ai, x, y, z, Qs[i],    q, nw );
        printf( " %i %f %f %f %f    %f %i \n", atypes[i], xyzs[i].x, xyzs[i].y, xyzs[i].z, Qs[i],    q, nw );
    }
    fclose (pFile);
    return 0;
}

MoleculeType::MoleculeType( char const* filename ){ loadFromFile_bas( filename ); };

int MoleculeType::findBonds( double sc ){
    nbonds = 0;
    int * bonds_ = new int[ 2*natoms*natoms ];
    for (int i=0; i<natoms; i++ ){
        //printf( "%i %i %i \n", i, atypes[i], natoms );
        for (int j=i+1; j<natoms; j++ ){
            double rijmax = sc * ( typeList->vdwRs[atypes[i]] + typeList->vdwRs[atypes[j]] );
            Vec3d d; d.set_sub( xyzs[i], xyzs[j] );
            double rij2 = d.norm2();
            //printf( " %i %i %i %i %f %f    %f %f   %f %f %f \n", i, j, atypes[i], atypes[j], rij2, rijmax, typeList->vdwRs[atypes[i]],typeList->vdwRs[atypes[j]], d.x, d.y, d.z  );
            if( rij2 < ( rijmax*rijmax ) ){
                bonds_[nbonds]=i; bonds_[nbonds+1]=j;
                nbonds+=2;
            }
        }
    }
    if ( bonds != NULL ){ delete bonds; printf( " delete bonds \n"); }
    bonds = new  int[ nbonds ];
    for (int i=0; i<nbonds; i++ ){  bonds[i] = bonds_[i]; };
    delete bonds_;
    return nbonds;
}

/*
int MoleculeType::drawAtom( int i, int nsphere, float atomscale, Uint32 color ){
    setColorInt32( color );
    int nvert = drawSphere_oct( nsphere, atomscale*typeList->vdwRs[atypes[i]], xyzs[i] );
    return nvert;
}

int MoleculeType::drawBond( int i, int j, int nstick, float bondwidth  ){
    Vec3f ai,aj;
    convert( xyzs[i], ai );
    convert( xyzs[j], aj );
    int nvert = drawCylinderStrip( nstick, bondwidth, bondwidth, ai, aj );
    return nvert;
}

int MoleculeType::makeViewCPK ( int nsphere, int nstick, float atomscale, float bondwidth ){
    if( viewlist > 0 ) {	glDeleteLists( viewlist, 1 );	}
    int nvert = 0;
    viewlist = glGenLists(1);
    glNewList( viewlist , GL_COMPILE );
        glShadeModel ( GL_SMOOTH );
        for( int i=0; i<natoms; i++ ){	nvert+= drawAtom( i, nsphere, atomscale, typeList->colors[ atypes[i] ] );	}
        glColor3f( 0.2f, 0.2f, 0.2f );
        for( int ib=0; ib<nbonds; ib+=2 ){	nvert+= drawBond( bonds[ib], bonds[ib+1], nstick, bondwidth );	}
    glEndList();
    printf( " nvert %i \n", nvert );
    return viewlist;
}
*/

void MoleculeType::toCOG_minmax(){
    Vec3d pmin,pmax; pmin.set( 100000,100000,100000 );  pmax.set( -100000,-100000,-100000 );
    for (int i=0; i<natoms; i++){
        pmin.x = fmin( pmin.x, xyzs[i].x  ); pmin.y = fmin( pmin.y, xyzs[i].y ); pmin.z = fmin( pmin.z, xyzs[i].z );
        pmax.x = fmax( pmax.x, xyzs[i].x  ); pmax.y = fmax( pmax.y, xyzs[i].y ); pmax.z = fmax( pmax.z, xyzs[i].z );
    }
    pmin.add(pmax); pmin.mul(0.5);
    for (int i=0; i<natoms; i++){
        xyzs[i].sub( pmin );
    }
}

void MoleculeType::toCOG_average(){
    Vec3d cog; cog.set(0.0);
    for (int i=0; i<natoms; i++){  cog.add( xyzs[i] ); }
    cog.mul( 1.0/natoms );
    for (int i=0; i<natoms; i++){  xyzs[i].sub( cog ); }
}
