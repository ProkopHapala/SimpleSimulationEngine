#ifndef MolecularDraw_h
#define MolecularDraw_h

#include "AtomicConfiguration.h"

void colorRB( float f ){ glColor3f( 0.5+f, 0.5, 0.5-f ); }

void printPoses( int n, double * poses ){
    for( int i=0; i<n; i++ ){
        int i8 = i*8;
        //printf( "force[%04i] %g,%g,%g,%g|%g,%g,%g,%g\n",i, opt.force[i8+0], opt.force[i8+1], opt.force[i8+2], opt.force[i8+3],    opt.force[i8+4], opt.force[i8+5], opt.force[i8+6], opt.force[i8+7]  );
        printf( "[%04i] %g,%g,%g,%g | %g,%g,%g,%g \n",i, poses[i8+0], poses[i8+1], poses[i8+2], poses[i8+3],    poses[i8+4], poses[i8+5], poses[i8+6], poses[i8+7]  );
    }
}

void drawMapedPoints( const FastAtomicMetric& D, int itest ){
    //atomdist.natoms=1;
    //atomdist.pos[0]=cursor3D;
    //atomdist.toCells(atomdist.ruler.step*0.5-0.01);
    Draw3D::drawBBox( D.ruler.pos0, D.ruler.pmax );
    int j=0;
    for(int i=0; i<D.natoms; i++){
        //Draw3D::drawPointCross( atomdist.pos[i], atomdist.Rcut );
        //Draw3D::drawPointCross( atomdist.pos[i], 0.1 );
        bool b = ( i == (itest%D.natoms));
        if(b){ Draw3D::drawSphereOctLines( 16, D.Rcut, D.pos[i] ); }
        else { Draw3D::drawPointCross( D.pos[i], 0.1 ); }
        //printf("%i %i \n", i, D.atomNs[i] );
        for(int jj=0; jj<D.atomNs[i];jj++){
            if(b){
                int ic = D.atom2cells[j];
                Vec3i ip;  D.ruler.i2ixyz ( ic, ip );
                Vec3d  p = D.ruler.box2pos( ip, {0.0,0.0,0.0} );
                double d = D.ruler.step;
                Draw3D::drawBBox( p, p+(Vec3d){d,d,d} );
            }
            j++;
        }
    }
}

void drawNeighs( const FastAtomicMetric& D, Vec3d pos ){
    Draw3D::drawBBox( D.ruler.pos0, D.ruler.pmax );
    Draw3D::drawSphereOctLines(16,D.Rcut,pos);
    {
        //ip = atomdist.ruler.i cursor3D
        //Vec3i ip;  atomdist.ruler.i2ixyz ( icell, ip );
        Vec3i ip = D.ruler.ipcell( pos );
        Vec3d  p = D.ruler.box2pos( ip, {0.0,0.0,0.0} );
        double d = D.ruler.step;
        Draw3D::drawBBox( p, p+(Vec3d){d,d,d} );
    }
    //printf( "DEBUG 2 \n" );
    if( Box::pointIn( pos, D.ruler.pos0, D.ruler.pmax) ){
        int tmpIs[D.natoms];
        int nfound = D.findNeighs( pos, D.Rcut, tmpIs );
        //printf( "DEBUG 3 \n" );
        //printf( "nfound %i \n", nfound );
        for(int i=0; i<nfound; i++){
            Draw3D::drawLine( pos, D.pos[tmpIs[i]] );
        }
    }
    for(int i=0; i<D.natoms; i++){
        //Draw3D::drawPointCross( atomdist.pos[i], atomdist.Rcut );
        Draw3D::drawPointCross( D.pos[i], 0.1 );
    }
}


void drawPPRelaxTrj( int n, double dt, double damp, GridFF& gff, Vec3d pos, Vec3d PRQ ){
    Vec3d vel = (Vec3d){0.0,0.0,0.0};
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++){
        Vec3d f = (Vec3d){0.0,0.0,0.0};
        gff.addForce( pos, PRQ, f);
        vel.mul(damp);
        vel.add_mul( f  , dt);
        pos.add_mul( vel, dt );
        glVertex3f(pos.x,pos.y,pos.z);
        //printf( " %i (%g,%g,%g) (%g,%g,%g) \n", i, pos.x,pos.y,pos.z,  f.x,f.y,f.z );
    }
    glEnd();
    //exit(0);
}

void drawGridForceAlongLine( int n, GridFF& gff, Vec3d pos0, Vec3d dpos, Vec3d PRQ, double fsc ){
    Vec3d pos = pos0;
	for( int i=0; i<n; i++ ){
        Vec3d f = (Vec3d){0.0,0.0,0.0};
        gff.addForce( pos, PRQ, f);
        //printf( " %i (%g,%g,%g) (%g,%g,%g) \n", i, pos.x,pos.y,pos.z,  f.x,f.y,f.z );
        Draw3D::drawVecInPos( f*fsc, pos );
        Draw3D::drawPointCross( pos, 0.1 );
        pos.add(dpos);
	}
}

void plotSurfPlane( Vec3d normal, double c0, Vec2d d, Vec2i n ){
    Vec3d da,db;
    normal.getSomeOrtho( da,db );
    da.mul( d.a/da.norm() );
    db.mul( d.b/db.norm() );
    //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(normal, {0.0,0.0,0.0} );
    //glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos(da*10, {0.0,0.0,0.0} );
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(db*10, {0.0,0.0,0.0} );
    Draw3D::drawRectGridLines( n*2, (da*-n.a)+(db*-n.b) + normal*c0, da, db );
}

/*
void renderSubstrate( int n, Vec3d * points, GLenum mode ){
    //printf( "iso_points.size() %i \n", iso_points.size() );
    if( mode == GL_POINTS ){
        glBegin(GL_POINTS);
        for(int i=0; i<iso_points.size(); i++){ glVertex3f( iso_points[i].x, iso_points[i].y, iso_points[i].z      ); }
        glEnd();
    }
}
*/

void renderSubstrate_( const GridShape& grid, Vec3d * FF, double isoval, bool sign ){
    //printf( "iso_points.size() %i \n", iso_points.size() );
    int nxy = grid.n.x * grid.n.y;
    printf("nxy %i \n", nxy );
    Vec3d * pos     = new Vec3d[nxy];
    Vec3d * normals = new Vec3d[nxy];
    //printf( " -- DEBUG 1 \n" );
    //DEBUG
    getIsoSurfZ( grid, isoval, sign, FF, pos, normals );
    //printf( " -- DEBUG 2 \n" );
    //glEnable(GL_LIGHTING);
    //DEBUG
    for ( int ib=1; ib<grid.n.y; ib++ ){
        glBegin(GL_TRIANGLE_STRIP);
        for ( int ia=0; ia<grid.n.x; ia++ ){
            int ip1 = (ib-1)*grid.n.x + ia;
            int ip2 = (ib  )*grid.n.x + ia;
            //printf( "iba (%i,%i) pos (%g,%g,%g)\n", ib,ia, pos[ip1].x,pos[ip1].y,pos[ip1].z );
            //glColor3f(pos[ip1].z*5-2,1.0f,1.0f); glNormal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); glVertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //glColor3f(pos[ip2].z*5-2,1.0f,1.0f); glNormal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); glVertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
            glColor3f(0.7f,0.7f,0.7f); glNormal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); glVertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            glColor3f(0.8f,0.7f,0.7f); glNormal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); glVertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
        }
        glEnd();
    }
    DEBUG
    //printf( " -- DEBUG 3 \n" );
    delete [] pos;
    delete [] normals;
    //exit(0);
}

void renderSubstrate_( const GridShape& grid, Vec3d * FF, Vec3d * FFel, double isoval, bool sign ){
    //printf( "iso_points.size() %i \n", iso_points.size() );
    Vec3d * pos     = new Vec3d[grid.n.x * grid.n.y];
    Vec3d * normals = new Vec3d[grid.n.x * grid.n.y];
    //printf( " -- DEBUG 1 \n" );
    getIsoSurfZ( grid, isoval, sign, FF, pos, normals );

    //printf( " -- DEBUG 2 \n" );

    //glEnable(GL_LIGHTING);
    for ( int ib=1; ib<grid.n.y; ib++ ){
        glBegin(GL_TRIANGLE_STRIP);
        for ( int ia=0; ia<grid.n.x; ia++ ){
            int ip1 = (ib-1)*grid.n.x + ia;
            int ip2 = (ib  )*grid.n.x + ia;
            //printf( "iba (%i,%i) pos (%g,%g,%g)\n", ib,ia, pos[ip1].x,pos[ip1].y,pos[ip1].z );
            //glColor3f(pos[ip1].z*5-2,1.0f,1.0f); glNormal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); glVertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //glColor3f(pos[ip2].z*5-2,1.0f,1.0f); glNormal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); glVertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);

            Vec3d gpos,fel1,fel2;
            grid.cartesian2grid( pos[ip1], gpos); fel1 = interpolate3DvecWrap( FFel, grid.n, gpos );
            grid.cartesian2grid( pos[ip2], gpos); fel2 = interpolate3DvecWrap( FFel, grid.n, gpos );

            //glColor3f(0.7f,0.7f,0.7f); glNormal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); glVertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //glColor3f(0.8f,0.7f,0.7f); glNormal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); glVertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);

            //glColor3f( fel1.x, fel1.y, fel1.z ); glNormal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); glVertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //glColor3f( fel2.x, fel2.y, fel2.z ); glNormal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); glVertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);

            colorRB( fel1.z ); glNormal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); glVertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            colorRB( fel2.z ); glNormal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); glVertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
        }
        glEnd();
    }
    //printf( " -- DEBUG 3 \n" );
    delete [] pos;
    delete [] normals;
    //exit(0);
}

void viewSubstrate( int nx, int ny, int isoOgl, Vec3d a, Vec3d b ){
    glPushMatrix();
    for( int ix = -nx; ix<=nx; ix++ ){
        for( int iy = -ny; iy<=ny; iy++ ){
            Vec3d pos = a*ix + b*iy;
            glTranslatef(pos.x, pos.y, pos.z);
            glCallList(isoOgl);
            glTranslatef(-pos.x, -pos.y, -pos.z);
        }
    }
    glPopMatrix();
}

#endif
