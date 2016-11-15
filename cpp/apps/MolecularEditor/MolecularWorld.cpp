
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "forceField.h"

#include  "Draw3D.h"

#include "MolecularWorld.h" // THE HEADER


// =========================================
// =========== Rotation optimization
// =========================================

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

int MolecularWorld::applyLinkerForce( ){
    for( int il=0; il<nLinkers; il++ ){
        Mat3d T;
        Vec3d dgpi,dgpj, dp;

        MolecularLink& li = linkers[il];
        int i = li.i;
        rot[i].toMatrix( T);
        T.dot_to( li.posi, dgpi );

        int j = li.j;
        rot[j].toMatrix(T);
        T.dot_to( li.posj, dgpj );

        dp       = dgpi + pos[i] - dgpj - pos[j];
        //Vec3d f = linkers[il].getForce(dp);
        Vec3d f = radialSpringForce( dp, li.k, li.l0 );



        //Draw3D::drawLine( dgpi + pos[i],   dgpj + pos[j] );

        //printf( "%i   %g %g %g   %g %g %g  \n", i, f.x,f.y,f.z, dp.x,dp.y,dp.z   );

        // TODO - we can optimize this if we use dGlobPos instead of LocPos
        fpos[i].add( f );   rot[i].addForceFromPoint( li.posi, f, frot[i] );
        f.mul(-1);
        fpos[j].add( f );   rot[j].addForceFromPoint( li.posj, f, frot[j] );
    }
    return nLinkers;
}

// BUG TO DO : when anchoring atom is not in COG simulation gets unstable !!!   -- tested on Si68+Butyl
int MolecularWorld::applyBondForce( ){
    for( int ib=0; ib<nBonds; ib++ ){
        Mat3d T;
        Vec3d lpi,lpj,dgpi,dgpj, dp;

        MolecularBond& bi =  bonds[ib];

        int i = bi.imol;
        rot[i].toMatrix( T);
        lpi = instances[i]->xyzs[bi.iatom];
        T.dot_to( lpi, dgpi );

        int j = bi.jmol;
        rot[j].toMatrix( T);
        lpj = instances[j]->xyzs[bi.jatom];
        T.dot_to( lpj, dgpj );

        dp       = dgpi + pos[i] - dgpj - pos[j];
        int ityp = instances[i]->atypes[bi.iatom];
        int jtyp = instances[j]->atypes[bi.jatom];

        int iff =  ityp*atomTypes.ntypes + jtyp;
        //printf( "applyBondForce: %i %i  %i %i  %i  %g %g \n", bi.iatom, bi.jatom, ityp, jtyp, iff, atomTypes.ntypes,  C6s[iff], C12s[iff] );
        Vec3d  f;
        if(nonCovalent){
            forceLJE( dp, -C6s[iff], -C12s[iff], 0.0, f );   // undo non-covalent force
            //f.add( radialSpringForce( dp, bi.k, bi.l0 ) );
            f.add( radialBondForce( dp, bi.k, bi.l0, bi.dlmax) );
        }else{
            f = radialBondForce( dp, bi.k, bi.l0, bi.dlmax);
        };

        //if(fmax<1e-4) printf(" ib %i r %f  dp (%g,%g,%g)\n", ib, dp.norm(), dp.x, dp.y, dp.z );

        //Draw3D::drawLine( dgpi + pos[i],   dgpj + pos[j] );

        //printf( "%i   %g %g %g   %g %g %g  \n", i, f.x,f.y,f.z, dp.x,dp.y,dp.z   );

        // TODO - we can optimize this if we use dGlobPos instead of LocPos
        fpos[i].add( f );   rot[i].addForceFromPoint( lpi, f, frot[i] );
        f.mul(-1);
        fpos[j].add( f );   rot[j].addForceFromPoint( lpj, f, frot[j] );
    }
    return nLinkers;
}



int MolecularWorld::checkBonds( double flmin, double flmax ){
    printf("checking bond lengths\n");
    for( int ib=0; ib<nBonds; ib++ ){
        Mat3d T;
        Vec3d lpi,lpj,dgpi,dgpj, dp;

        MolecularBond& bi =  bonds[ib];

        int i = bi.imol;
        rot[i].toMatrix( T);
        lpi = instances[i]->xyzs[bi.iatom];
        T.dot_to( lpi, dgpi );

        int j = bi.jmol;
        rot[j].toMatrix( T);
        lpj = instances[j]->xyzs[bi.jatom];
        T.dot_to( lpj, dgpj );

        dp       = dgpi + pos[i] - dgpj - pos[j];
        double l = dp.norm();
        if( l > flmax*bi.l0 ) printf( "bond too long  %i %i  %i %i | %g instead %g %g\n", ib, i, j, bi.iatom, bi.jatom, l, bi.l0, l/bi.l0 );
        if( l < flmin*bi.l0 ) printf( "bond too short %i %i  %i %i | %g instead %g %g\n", ib, i, j, bi.iatom, bi.jatom, l, bi.l0, l/bi.l0 );
        //printf( "%i %i  %i %i %g %g\n", ib, i, j, bi.iatom, bi.jatom, dp.norm(), bi.l0 );
    }
    return nLinkers;
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
    // points
    if( nonCovalent ){
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
    } }
    // linkers
    if(linkers) nInteractions += applyLinkerForce( );
    if(bonds)   nInteractions += applyBondForce( );
    //exit(0);
    /*
    for(int i=0; i<nmols; i++){
        printf( "%i %g %g %g    %g %g %g %g\n", i, fpos[i].x,fpos[i].y,fpos[i].z,   frot[i].x,frot[i].y,frot[i].z,frot[i].w );
    }
    exit(0);
    */
}

void MolecularWorld::rigidOptStep( ){
    //printf("DEBUG 1\n");
    for (int i=0; i<nmols; i++){
        fpos[i].set(0.0d); frot[i].set(0.0d,0.0d,0.0d,0.0d);   // set all forces to zero
        rot[i].normalize();                               // keep quaternion normalized, otherwise unstable !!!
    }
    //printf("DEBUG 2\n");
    assembleForces( );
    //printf("DEBUG 3\n");
    for (int i=0; i<nmols; i++){
        // out-project component which would harm unitarity of quaternion
        double cdot;
        cdot = rot[i].dot( frot[i] );  frot[i].add_mul( rot[i], -cdot );
        cdot = rot[i].dot( vrot[i] );  vrot[i].add_mul( rot[i], -cdot );
        // apply constrains
        int i2 = i<<1;
        if(constrains[i2  ]){ fpos[i].set(0.0d); vpos[i].set(0.0d); }
        if(constrains[i2+1]){ frot[i].set(0.0d); vrot[i].set(0.0d); }
    }
    //printf("DEBUG 4\n");
    //optimizer->move();
    //optimizer->optStep();
    //for(int i=0;i<optimizer->n; i++ ){ printf( " %i %g %g %g \n", i, optimizer->pos[i], optimizer->vel[i], optimizer->force[i] ); }
    optimizer->move_FIRE();
    fmax = optimizer->getFmaxAbs( );
    //printf("DEBUG 5\n");
}

// =========================================
// =========== initialization and I/O
// =========================================

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
    invMpos  = (Vec3d* )&optimizer->invMasses[0];

    rot   = (Quat4d*)&optimizer->pos[irot];
    vrot  = (Quat4d*)&optimizer->vel[irot];
    frot  = (Quat4d*)&optimizer->force[irot];
    invMrot  = (Quat4d*)&optimizer->invMasses[irot];

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
    fclose(pFile);
    return nMolTypes;
}

int MolecularWorld::loadInstances( char const* fileName ){
    printf(" loading molecular instances from: >>%s<<\n", fileName );

    FILE * pFile;
    pFile = fopen (fileName,"r");
    fscanf ( pFile, " %i", &nmols);
    printf ( "nmols %i \n", nmols );
    //printf ( " PhysicalSystem 1 \n" );

    instances  = new MoleculeType*[nmols];
    constrains = new bool       [2*nmols];
    initParams( );
    //printf( " PhysicalSystem 2 \n" );

    for (int i=0; i<nmols; i++){
        int i2 = i<<1;
        int itype;
        Mat3d M;
        //int c1,c2;
        fscanf (pFile, " %i   %lf %lf %lf     %lf %lf %lf     %lf %lf %lf  %i %i\n",
            &itype,
            &pos[i].x, &pos[i].y, &pos[i].z,
            &M.ax, &M.ay, &M.az,
            &M.bx, &M.by, &M.bz,
            &constrains[i2], &constrains[i2+1]
        );
        printf(  " %i   %f %f %f     %f %f %f     %f %f %f %i %i\n", itype, pos[i].x, pos[i].y, pos[i].z,    M.ax, M.ay, M.az,    M.bx, M.by, M.bz, constrains[i2], constrains[i2+1] );
        itype--;
        instances[i] = &molTypes[itype];

        //printf("M.a %f %f %f\n", M.ax,M.ay,M.az);
        //printf("M.b %f %f %f\n", M.bx,M.by,M.bz);

        M.a.normalize();
        M.b.add_mul( M.a, -M.a.dot( M.b ) );
        M.b.normalize();
        M.c.set_cross( M.a, M.b );
        //printf("M.a %f %f %f\n", M.ax,M.ay,M.az);
        //printf("M.b %f %f %f\n", M.bx,M.by,M.bz);
        //printf("M.c %f %f %f\n", M.cx,M.cy,M.cz);
        rot[i].fromMatrix( M );

        invMpos[i].set(1.0d);
        invMrot[i].set(1.0d/instances[i]->Rmax);

    }
    fclose(pFile);
    return nmols;
}

int MolecularWorld::loadLinkers( char const* fileName ){
    printf(" loading linkers from: >>%s<<\n", fileName );
    FILE * pFile;
    pFile = fopen (fileName,"r");
    fscanf ( pFile, " %i\n", &nLinkers);
    linkers = new MolecularLink[nLinkers];
    for (int i=0; i<nLinkers; i++){
        MolecularLink& li = linkers[i];
        fscanf (pFile, " %i %i   %lf %lf %lf     %lf %lf %lf     %lf %lf\n",
            &li.i, &li.j,
            &li.posi.x, &li.posi.y, &li.posi.z,
            &li.posj.x, &li.posj.y, &li.posj.z,
            &li.k, &li.l0
        );
        //li.i--; li.j--;  // uncoment this if instances numbered from 1 rather than from 0
        printf ( " %i %i   %lf %lf %lf     %lf %lf %lf     %lf %lf\n", li.i, li.j,  li.posi.x, li.posi.y, li.posi.z, li.posj.x, li.posj.y, li.posj.z, li.k, li.l0 );
    }
    fclose(pFile);
    return nLinkers;
}


int MolecularWorld::loadBonds( char const* fileName ){
    printf(" loading bonds from: >>%s<<\n", fileName );
    FILE * pFile;
    pFile = fopen (fileName,"r");
    fscanf ( pFile, " %i\n", &nBonds);
    bonds = new MolecularBond[nBonds];
    for (int i=0; i<nBonds; i++){
        MolecularBond& bi = bonds[i];
        fscanf (pFile, " %i %i   %i %i  %lf %lf %lf\n", &bi.imol, &bi.jmol,  &bi.iatom, &bi.jatom, &bi.k, &bi.l0, &bi.dlmax );
        //li.i--; li.j--;  // uncoment this if instances numbered from 1 rather than from 0
        printf ( " %i %i   %i %i   %lf %lf %lf %lf\n", bi.imol, bi.jmol, bi.iatom, bi.jatom, bi.k, bi.l0, bi.dlmax );
    }
    fclose(pFile);
    return nLinkers;
}

bool MolecularWorld::fromDir( char const* dirName, char const* atom_fname, char const* mol_fname, char const* instance_fname ){

    printf("dirName: >>%s<< atom_fname: >>%s<< mol_fname: >>%s<< instance_fname: >>%s<<\n", dirName, atom_fname, mol_fname, instance_fname );

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

int  MolecularWorld::saveInstances( char const* fileName ){
    //printf(" loading molecular instances from: >>%s<<\n", fileName );
    FILE * pFile;
    pFile = fopen (fileName,"w");
    fprintf ( pFile, "%i\n", nmols );
    for (int i=0; i<nmols; i++){
        int i2 = i<<1;
        Mat3d M;
        rot[i].toMatrix(M);
        int itype = instances[i] - molTypes + 1 ; // TO DO : this is a bit strange hack
        fprintf (pFile, " %i   %lf %lf %lf     %lf %lf %lf     %lf %lf %lf  %i %i\n",
            itype,
            pos[i].x, pos[i].y, pos[i].z,
            M.ax, M.ay, M.az,
            M.bx, M.by, M.bz,
            constrains[i2], constrains[i2+1]
        );
    }
    fclose(pFile);
    return nmols;
};

int MolecularWorld::exportAtomsXYZ(  FILE * pFile, const char * comment ){
    //printf( "exportAtomsXYZ \n");
    int natoms = 0;
    for (int i=0; i<nmols; i++){ natoms += instances[i]->natoms; }
    fprintf(pFile, "%i\n", natoms  );
    fprintf(pFile, "%s\n", comment );
    for (int i=0; i<nmols; i++){
        MoleculeType * moli = instances[i];
        int npi = moli->natoms;
        transformPoints( pos[i], rot[i], npi, moli->xyzs, Tps_i );
        for (int j=0; j<npi; j++){
            fprintf( pFile, " %s %3.6f %3.6f %3.6f\n", atomTypes.names[moli->atypes[j]], Tps_i[j].x, Tps_i[j].y, Tps_i[j].z );
            //printf( "DEBUG %s %3.6f %3.6f %3.6f\n", atomTypes.names[moli->atypes[j]], Tps_i[j].x, Tps_i[j].y, Tps_i[j].z );
        }
    }
    return natoms;
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
