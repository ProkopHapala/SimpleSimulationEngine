
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "forceField.h"

#include "MolecularWorld.h" // THE HEADER


void MolecularWorld::initParams( ){
    // optimization parameters
    int nparams = (3+4)*nmols;
    //optimizer   = new OptimizerFIRE( nparams, new double[nparams], new double[nparams], new double[nparams], NULL );
    //optimizer   = new OptimizerDerivs( nparams, new double[nparams], new double[nparams], new double[nparams], NULL );
    //optimizer->dt = 0.05;
    //optimizer->damping = 0.05;
    optimizer = new DynamicOpt();
    optimizer->allocate( nparams );
    //optimizer->bindArrays( int n_, double * pos_, double * vel_, double * force_ );
    optimizer->initOpt( 0.01, 0.1 );

    int irot = 3*nmols;
    pos   = (Vec3d* )&optimizer->pos[0];
    vpos  = (Vec3d* )&optimizer->vel[0];
    fpos  = (Vec3d* )&optimizer->force[0];
    rot   = (Quat4d*)&optimizer->pos[irot];
    vrot  = (Quat4d*)&optimizer->vel[irot];
    frot  = (Quat4d*)&optimizer->force[irot];
}

void MolecularWorld::initTPoints(){
    // point temporary
    nptmp = 0;
    for (int i=0; i<nMolTypes; i++){ if( molTypes[i].natoms > nptmp ) nptmp = molTypes[i].natoms; } // size according to larges molecule type
    Tps_i = new Vec3d[nptmp];
    Tps_j = new Vec3d[nptmp];
    fs_i  = new Vec3d[nptmp];
    fs_j  = new Vec3d[nptmp];
    //atypeList = moleculeTypes[0]->typeList;
}

int MolecularWorld::loadMolTypes( char const* dirName, char const* fileName ){
    char str  [256];
    char fname[256];
    FILE * pFile;
    strcpy(fname,dirName ); strcat(fname,fileName);
    pFile = fopen (fname,"r");
    printf(" loading molTypes from: >>%s<<\n", fname );
    fscanf ( pFile, "%i", &nMolTypes);
    molTypes = new MoleculeType[nMolTypes];
    for(int i=0; i<nMolTypes; i++){
        fscanf ( pFile, "%s", &str );
        strcpy(fname,dirName ); strcat(fname,str);
        //MoleculeType * molType = new MoleculeType();
        //molType->loadFromFile_bas( fname );
        //molType->typeList = &atomTypes;
        //moleculeTypes[i]= molType;
        molTypes[i].loadFromFile_bas( fname );
        molTypes[i].typeList = &atomTypes;
    }
    return nMolTypes;
}

int MolecularWorld::loadInstances( char const* filename ){
    printf(" loading molecular instances from: >>%s<<\n", filename );

    FILE * pFile;
    pFile = fopen (filename,"r");
    fscanf ( pFile, " %i", &nmols);
    printf ( "nmols %i \n", nmols );
    //printf ( " PhysicalSystem 1 \n" );
    instances = new MoleculeType*[nmols];
    initParams( );
    //printf( " PhysicalSystem 2 \n" );

    for (int i=0; i<nmols; i++){
        int itype;
        Mat3d M;
        fscanf (pFile, " %i %lf %lf %lf     %lf %lf %lf     %lf %lf %lf", &itype, &pos[i].x, &pos[i].y, &pos[i].z,    &M.ax, &M.ay, &M.az,    &M.bx, &M.by, &M.bz );
        printf(  " %i %f %f %f     %f %f %f     %f %f %f \n", itype, pos[i].x, pos[i].y, pos[i].z,    M.ax, M.ay, M.az,    M.bx, M.by, M.bz );
        itype--;
        instances[i] = &molTypes[itype];

        M.a.normalize();
        M.b.add_mul( M.a, -M.a.dot( M.b ) );
        M.b.normalize();
        M.c.set_cross( M.a, M.b );
        rot[i].fromMatrix( M );
/*
        printf(" matrix:  \n");
        printVec(M.a);
        printVec(M.b);
        printVec(M.c);
        printf(" ortonormality:  \n");
        printf( " a.a b.b c.c %f %f %f \n", M.a.norm2(),  M.b.norm2(),  M.c.norm2() );
        printf( " a.b a.c b.c %f %f %f \n", M.a.dot(M.b), M.a.dot(M.c), M.b.dot(M.c) );
        printf( " det(M) %f \n", M.determinant() );
        printf(" matrix -> quat \n");
        rot[i].fromMatrix( M );
        printf(" q: %f %f %f %f  qnorm %f  \n", rot[i].x, rot[i].y, rot[i].z, rot[i].w,   rot[i].norm2() );
        printf(" quat -> matrix  \n");
        rot[i].toMatrix( M );
        printVec(M.a);
        printVec(M.b);
        printVec(M.c);
*/
       // molecules[ i ] = &molTypeList[ itype-1 ];
        //printf("=====\n");
    }
    return nmols;
}

bool MolecularWorld::fromDir( char const* dirName, char const* atom_fname, char const* mol_fname, char const* instance_fname ){

    //printf("dirName: >>%s<<\n", dirName);
    char fname[256];
    strcpy(fname,dirName ); strcat(fname,atom_fname);
    atomTypes.loadFromFile( fname );
    loadMolTypes(dirName,mol_fname);
    strcpy(fname,dirName ); strcat(fname,instance_fname);
    loadInstances( fname );

    //printf("dirName: >>%s<< fileName: >>%s<<\n", dirName, fname );

    //fclose (pFile);
    //printf( "PhysicalSystem 3 \n" );
    initTPoints();
    //printf( "PhysicalSystem 4 \n" );
    return true;
}

/*
MolecularWorld::MolecularWorld( int nmols_, MoleculeType ** molecules_ ){
    nmols = nmols_;
    molecules = molecules_;
    initParams ( );
    initTPoints( );
}
*/

void MolecularWorld::makeFF( ){
    makeLJparams( atomTypes.ntypes, atomTypes.vdwRs, atomTypes.vdwEs, C6s, C12s );
}


// =========== Rotation optimization

void MolecularWorld::transformPoints( const Vec3d& pos, const Quat4d& rot, int npoints, Vec3d * points, Vec3d * Tpoints ){
    Mat3d T;
    //q.toMatrix_unitary2( T );
    rot.toMatrix( T);
    //printVec( T.a );
    //printVec( T.b );
    //printVec( T.c );
    //printf( " %f %f %f \n", T.a.dot(T.b), T.a.dot(T.c), T.b.dot(T.c) );
    //printf( " %f %f %f \n", T.a.dot(T.a), T.b.dot(T.b), T.c.dot(T.c) );
    for( int i=0; i<npoints; i++ ){
        Vec3d Tp;
        T.dot_to_T(   points[i],  Tp );
        //T.dot_to(   points[i],  Tp );
        Tpoints[i].set_add( pos, Tp  );
    }
}

void MolecularWorld::forceFromPoints( int npoints, Vec3d * points, Vec3d * forces,  const Quat4d& q,  Vec3d& fp, Quat4d& fq ){
    //printf( "forceFromPoints ----------------\n" );
    for( int i=0; i<npoints; i++ ){
        //printf( " %i   %f %f %f   %f %f %f \n", i,   points[i].x, points[i].y, points[i].z,   forces[i].x, forces[i].y, forces[i].z  );
        q .addForceFromPoint( points[i], forces[i], fq );
        fp.add( forces[i] );
        //printf( " %i   %f %f %f   %f %f %f %f \n", i,   forces[i].x, forces[i].y, forces[i].z,    fp.x, fp.y, fp.z,  fq.x, fq.y, fq.z, fq.w  );
    }
}

void MolecularWorld::cleanPointForce( int npoints, Vec3d * forces ){	for( int i=0; i<npoints; i++ ){  forces[i].set(0.0);  }	}

void MolecularWorld::assembleForces( ){
    nInteractions = 0;
    for (int i=0; i<nmols; i++){
        //printf("DEBUG 2.1\n");
        MoleculeType * moli = instances[i];
        int npi = moli->natoms;
        //printf("DEBUG 2.2\n");
        transformPoints( pos[i], rot[i], npi, moli->xyzs, Tps_i );
        cleanPointForce( moli->natoms, fs_i );
        //printf("DEBUG 2.3\n");
        for (int j=0; j<i; j++){
            MoleculeType * molj = instances[j];
            int npj = molj->natoms;
            //printf("DEBUG 2.4\n");
            transformPoints( pos[j], rot[j], npj, molj->xyzs, Tps_j );
            //cleanPointForce( moli->natoms, fs_i );
            //printf("DEBUG 2.5\n");
            cleanPointForce( molj->natoms, fs_j );
            //printf("DEBUG 2.6\n");

            nInteractions +=
            interMolForceLJE(
                npi, moli->atypes, moli->Qs, Tps_i, fs_i,
                npj, molj->atypes, molj->Qs, Tps_j, fs_j,
                atomTypes.ntypes, C6s, C12s
            );
            //printf("DEBUG 2.7\n");
            forceFromPoints( npj, molj->xyzs, fs_j,  rot[j], fpos[j], frot[j] );
            //forceFromPoints( moli->natoms, moli->xyzs, fs_i,  rot[i],  fpos[i], frot[i] );
            //exit(0);
        }
        //printf("DEBUG 2.8\n");
        //forceMolSurf( surf_z0, surf_zMin, surf_Emin, surf_hat, npi,  Tps_i, fs_i );
        //printf("DEBUG 2.9\n");
        forceFromPoints( npi, moli->xyzs, fs_i,  rot[i],  fpos[i], frot[i] );
        //printf("DEBUG 2.10\n");
    }
}

void MolecularWorld::rigidOptStep( ){
    //printf("DEBUG 1\n");
    for (int i=0; i<nmols; i++){
        fpos[i].set(0.0); frot[i].set(0.0,0.0,0.0,0.0);   // set all forces to zero
        rot[i].normalize();                               // keep quaternion normalized, otherwise unstable !!!
    }
    //printf("DEBUG 2\n");
    assembleForces( );
    //printf("DEBUG 3\n");
    for (int i=0; i<nmols; i++){
        double qfq = rot[i].dot( frot[i] );
        frot[i].add_mul( rot[i], -qfq ); // out-project component which would harm unitarity of quaternion
    }
    //printf("DEBUG 4\n");
    //optimizer->move();
    //optimizer->optStep();
    //for(int i=0;i<optimizer->n; i++ ){ printf( " %i %g %g %g \n", i, optimizer->pos[i], optimizer->vel[i], optimizer->force[i] ); }
    optimizer->move_FIRE();
    //printf("DEBUG 5\n");
}

// =========== view utils

/*
void draw(){
    for (int i=0; i<nmols; i++){
        if( molecules[i]->viewlist > 0 ){
            Mat3d rotmat;
            //rot[i].toMatrix_unitary2( rotmat );
            //rot[i].toMatrix_unitary( rotmat );
            rot[i].toMatrix( rotmat );
            float glMat[4*4];
            toGLMat( pos[i], rotmat, glMat );
            glPushMatrix();
            glMultMatrixf( glMat );
            glCallList   ( molecules[i]->viewlist );
            glPopMatrix();
        }
    };
}
*/
