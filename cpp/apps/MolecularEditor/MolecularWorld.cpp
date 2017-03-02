
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "forceField.h"

//#include  "Draw3D.h"

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
            forceLJE( dp, -ForceField::C6s[iff], -ForceField::C12s[iff], 0.0, f );   // undo non-covalent force
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

int MolecularWorld::nonBondingFroces_N2( ){
    int nInteractions_ = 0;
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

            nInteractions_ +=
            interMolForceLJE(
                npi, moli->atypes, moli->Qs, Tps_i, fs_i,
                npj, molj->atypes, molj->Qs, Tps_j, fs_j
                //atomTypes.ntypes, C6s, C12s
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
    return nInteractions_;
}



int MolecularWorld::nonBondingFroces_bbox( ){
    int nInteractions_ = 0;
    for (int i=0; i<nmols; i++){
        //printf("DEBUG 2.1\n");
        MoleculeType * moli = instances[i];
        int npi = moli->natoms;
        //printf("DEBUG 2.2\n");
        transformPoints( pos[i], rot[i], npi, moli->xyzs, Tps_i );
        cleanPointForce( moli->natoms, fs_i );
        //printf("DEBUG 2.3\n");

        Vec3d pmin,pmax;
        pmin.set(+1e+300,+1e+300,+1e+300);
        pmax.set(-1e+300,-1e+300,-1e+300);
        for(int iatom=0; iatom<npi; iatom++ ){
            Vec3d& p = Tps_i[iatom];
            if(p.x<pmin.x){ pmin.x=p.x;} if(p.x>pmax.x){ pmax.x=p.x;}
            if(p.y<pmin.y){ pmin.y=p.y;} if(p.y>pmax.y){ pmax.y=p.y;}
            if(p.z<pmin.z){ pmin.z=p.z;} if(p.z>pmax.z){ pmax.z=p.z;}
            //printf( "%i p: (%3.3f,%3.3f,%3.3f) \n", ia, p.x, p.y, p.z );
        }
        pmin.sub({Rcut*1.01,Rcut*1.01,Rcut*1.01});
        pmax.add({Rcut*1.01,Rcut*1.01,Rcut*1.01});

        for (int j=0; j<i; j++){
            MoleculeType * molj = instances[j];
            int npj = molj->natoms;
            //printf("DEBUG 2.4\n");
            transformPoints( pos[j], rot[j], npj, molj->xyzs, Tps_j );
            //cleanPointForce( moli->natoms, fs_i );
            //printf("DEBUG 2.5\n");
            cleanPointForce( molj->natoms, fs_j );
            //printf("DEBUG 2.6\n");

            for(int jatom=0; jatom<npj; jatom++ ){
                Vec3d& pj = Tps_j[jatom];
                if( (pj.x<pmin.x)||(pj.y<pmin.y)||(pj.z<pmin.z)||(pj.x>pmax.x)||(pj.y>pmax.y)||(pj.z>pmax.z) ) continue;

                int atyp  = molj->atypes[jatom];
                int ityp0 = ForceField::ntypes*atyp;
                double qa = molj->Qs[jatom];

                for(int iatom=0; iatom<npi; iatom++){

                    Vec3d dR;
                    dR.set_sub( Tps_i[iatom], pj );
                    //Draw3D::drawLine(Tps_i[iatom], pj);

                    nInteractions_++;
                    double r2 = dR.norm2();
                    if( r2>ForceField::Rcut2 ) continue;

                    //Draw3D::drawLine(Tps_i[iatom], pj);

                    int btyp   = moli->atypes[iatom];
                    int ityp   = ityp0 + btyp;
                    double C6  = ForceField::C6s [ ityp ];
                    double C12 = ForceField::C12s[ ityp ];
                    //printf( "interMolForceLJE %i %i  %i %i  %i %i  %g %g \n", ia, ib, atyp, btyp, ityp, ntypes,   C6, C12 );
                    Vec3d f;

                    //printf( " %i %i %i %i \n", ia, ib, atyp, btyp );

                    //forceLJE( dR, C6, C12, qa*Qbs[ib], f );

                    double qq   = qa*moli->Qs[iatom];
                    double fcut = ForceField::FcutLJ[ityp] + ForceField::FcutCoulomb[ityp]*qq;
                    //ForceField::forceLJE( dR, C6, C12, qq, -fcut, f );
                    ForceField::forceLJE( dR, C6, C12, qq, 0.0, f );

                    fs_i[iatom].add( f ); fs_j[jatom].sub( f );

                }
            }

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
    return nInteractions_;
}

int MolecularWorld::nonBondingFroces_buf( ){
    int nInteractions_ = 0;

    //Draw3D::drawBBox( {0.0,1.0,2.0}, {3.0,4.0,5.0} );

    for (int imol=0; imol<nmols; imol++){
        MoleculeType * moli = instances[imol];
        int npi = moli->natoms;
        transformPoints( pos[imol], rot[imol], npi, moli->xyzs, Tps_i );

        boxbuf.clear( );

        Vec3d pmin,pmax;
        pmin.set(+1e+300,+1e+300,+1e+300);
        pmax.set(-1e+300,-1e+300,-1e+300);
        for(int iatom=0; iatom<npi; iatom++ ){
            Vec3d& p = Tps_i[iatom];
            if(p.x<pmin.x){ pmin.x=p.x;} if(p.x>pmax.x){ pmax.x=p.x;}
            if(p.y<pmin.y){ pmin.y=p.y;} if(p.y>pmax.y){ pmax.y=p.y;}
            if(p.z<pmin.z){ pmin.z=p.z;} if(p.z>pmax.z){ pmax.z=p.z;}
            //printf( "%i p: (%3.3f,%3.3f,%3.3f) \n", ia, p.x, p.y, p.z );
        }
        boxbuf.pos0 = pmin; boxbuf.pos0.sub({Rcut*1.1,Rcut*1.1,Rcut*1.1});

        //printf( "pmin (%3.3f,%3.3f,%3.3f) pmax (%3.3f,%3.3f,%3.3f)\n", pmin.x, pmin.y, pmin.z, pmax.x, pmax.y, pmax.z );
        //printf( "N (%3.3f,%3.3f,%3.3f)\n", (pmax.x-pmin.x)*boxbuf.invStep, (pmax.x-pmin.x)*boxbuf.invStep, (pmax.x-pmin.x)*boxbuf.invStep );

        //for(int iatom=0; iatom<npi; iatom++ ){
        for(int iatom=0; iatom<1; iatom++ ){
            Vec3d dpos; Vec3i ipos;
            boxbuf.pos2box( Tps_i[iatom], ipos,dpos );
            //printf( "%i  (%3.3f,%3.3f,%3.3f) (%i,%i,%i)\n", ia, Tps_i[ia].x, Tps_i[ia].y, Tps_i[ia].z, ipos.x,ipos.y,ipos.z);
            boxbuf.insert( iatom, ipos, dpos, Rcut );
            //Draw3D::drawSphere_oct(3,Rcut,Tps_i[iatom]);
            //boxbuf.insert(ia,Tps_i[ia],Rcut);
        }
        //exit(0);

        /*
        for(int ibox=0; ibox<boxbuf.NXYZ; ibox++){
            int nib = boxbuf.counts[ibox];
            if(nib==0) continue;
            Vec3i ipos;
            boxbuf.i2xyz(ibox,ipos.x,ipos.y,ipos.z);
            Vec3d p0,p1;
            p0.set(ipos.x*boxbuf.step+boxbuf.pos0.x,
                   ipos.y*boxbuf.step+boxbuf.pos0.y,
                   ipos.z*boxbuf.step+boxbuf.pos0.z);
            p1.set_add( p0, {boxbuf.step,boxbuf.step,boxbuf.step} );
            Draw3D::drawBBox( p0, p1 );
            //printf("%i %i \n", ibox, nib );
            for( int im=0; im<nib; im++ ){
                int iatom = boxbuf.get(ibox,im);
                Draw3D::drawPointCross( Tps_i[iatom], 0.5 );
            }
        }
        */
        //exit(0);
        //return;
        cleanPointForce( moli->natoms, fs_i );

        for (int jmol=0; jmol<imol; jmol++){
            MoleculeType * molj = instances[jmol];
            int npj = molj->natoms;
            //printf("DEBUG 2.4\n");
            transformPoints( pos[jmol], rot[jmol], npj, molj->xyzs, Tps_j );
            //cleanPointForce( moli->natoms, fs_i );
            //printf("DEBUG 2.5\n");
            cleanPointForce( molj->natoms, fs_j );
            //printf("DEBUG 2.6\n");

            for(int jatom=0; jatom<npj; jatom++ ){
                Vec3d& pj = Tps_j[jatom];
                Vec3i ipos; Vec3d dpos;
                boxbuf.pos2box(pj,ipos,dpos);
                if( !boxbuf.validIndex( ipos ) ) continue;

                int atyp  = molj->atypes[jatom];
                int ityp0 = ForceField::ntypes*atyp;
                double qa = molj->Qs[jatom];

                int ixyz = boxbuf.xyz2i(ipos.x,ipos.y,ipos.z);
                int nib  = boxbuf.counts[ixyz];
                for(int im=0; im<nib; im++){
                    int iatom = boxbuf.get(ixyz,im);

                    Vec3d dR;
                    dR.set_sub( Tps_i[iatom], pj );
                    //Draw3D::drawLine(Tps_i[iatom], pj);

                    nInteractions_++;
                    double r2 = dR.norm2();
                    if( r2>ForceField::Rcut2 ) continue;

                    //Draw3D::drawLine(Tps_i[iatom], pj);

                    int btyp   = moli->atypes[iatom];
                    int ityp   = ityp0 + btyp;
                    double C6  = ForceField::C6s [ ityp ];
                    double C12 = ForceField::C12s[ ityp ];
                    //printf( "interMolForceLJE %i %i  %i %i  %i %i  %g %g \n", ia, ib, atyp, btyp, ityp, ntypes,   C6, C12 );
                    Vec3d f;

                    //printf( " %i %i %i %i \n", ia, ib, atyp, btyp );

                    //forceLJE( dR, C6, C12, qa*Qbs[ib], f );

                    double qq   = qa*moli->Qs[iatom];
                    double fcut = ForceField::FcutLJ[ityp] + ForceField::FcutCoulomb[ityp]*qq;
                    //ForceField::forceLJE( dR, C6, C12, qq, -fcut, f );
                    ForceField::forceLJE( dR, C6, C12, qq, 0.0, f );

                    fs_i[iatom].add( f ); fs_j[jatom].sub( f );

                }
            }

            forceFromPoints( npj, molj->xyzs, fs_j,  rot[jmol], fpos[jmol], frot[jmol] );
        }
        //printf("DEBUG 2.8\n");
        //forceMolSurf( surf_z0, surf_zMin, surf_Emin, surf_hat, npi,  Tps_i, fs_i );
        //printf("DEBUG 2.9\n");
        forceFromPoints( npi, moli->xyzs, fs_i,  rot[imol],  fpos[imol], frot[imol] );
    }
    return nInteractions_;
}

void MolecularWorld::assembleForces( ){
    nInteractions = 0;
    // points
    //if( nonCovalent ){ nInteractions += nonBondingFroces_N2(); }
    if( nonCovalent ){ nInteractions += nonBondingFroces_bbox(); }
    //if( nonCovalent ){ nInteractions += nonBondingFroces_buf(); }
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

    nAtomTot = 0;
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

        nAtomTot += instances[i]->natoms;

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
    if( !pFile ){ printf("cannot access >>%s<<\n", fileName); return 0; }
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

int MolecularWorld::loadSplines( char const* fileName ){
    char str  [256];
    printf(" loading splines from: >>%s<<\n", fileName );
    FILE * pFile;
    pFile = fopen (fileName,"r");
    fscanf ( pFile, " %i\n", &nSplines);
    splines = new GeneralSpline[nSplines];
    for (int i=0; i<nSplines; i++){
        fscanf ( pFile, "%s", &str );
        splines[i].loadFromFile(str);
    }
    fclose(pFile);
    return nSplines;
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




int MolecularWorld::getNAtoms(){
    int natoms = 0;
    for (int i=0; i<nmols; i++){ natoms += instances[i]->natoms; }
    return natoms;
}

int MolecularWorld::getAtomPos( Vec3d * buff ){
    //int natoms = getNAtoms();
    int iatom  = 0;
    for (int i=0; i<nmols; i++){
        MoleculeType * moli = instances[i];
        int npi = moli->natoms;
        transformPoints( pos[i], rot[i], npi, moli->xyzs, Tps_i );
        for (int j=0; j<npi; j++){
            buff[iatom].set( Tps_i[j].x, Tps_i[j].y, Tps_i[j].z );
            iatom++;
        }
    }
    return iatom;
}

int MolecularWorld::getAtomTypes( int * buff ){
    //int natoms = getNAtoms();
    int iatom  = 0;
    for (int i=0; i<nmols; i++){
        MoleculeType * moli = instances[i];
        int npi = moli->natoms;
        for (int j=0; j<npi; j++){
            buff[iatom] = moli->atypes[j];
            iatom++;
        }
    }
    return iatom;
}

int MolecularWorld::exportAtomsXYZ(  FILE * pFile, const char * comment ){
    //printf( "exportAtomsXYZ \n");
    int natoms = getNAtoms();
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




/*
MolecularWorld::MolecularWorld( int nmols_, MoleculeType ** molecules_ ){
    nmols = nmols_;
    molecules = molecules_;
    initParams ( );
    initTPoints( );
}
*/

void MolecularWorld::makeFF( ){
    //ForceField::init( 6.0, atomTypes.ntypes, atomTypes.vdwRs, atomTypes.vdwEs );
    ForceField::init( Rcut, atomTypes.ntypes, atomTypes.vdwRs, atomTypes.vdwEs );
/*
    makeLJparams( atomTypes.ntypes, atomTypes.vdwRs, atomTypes.vdwEs, C6s, C12s );
    ForceField::ntypes = atomTypes.ntypes;
    ForceField::C12s   = C12s;
    ForceField::C6s    = C6s;
*/

}

void MolecularWorld::setCutoff  ( double Rcut_){
    Rcut = Rcut_;
    boxbuf.step    = Rcut*2;
    boxbuf.invStep = 1/boxbuf.step;
    boxbuf.clear( );
    //ForceField::Rcut = Rcut; // TO DO : this is somwhat inconsistent
};


//void MolecularWorld::initBufBox( double Rmax ){}




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
