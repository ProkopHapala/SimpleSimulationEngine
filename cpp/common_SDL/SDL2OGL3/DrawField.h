
#ifndef  DrawField_h
#define  DrawField_h

typedef Vec3d (*VecFieldFunc)(Vec3d R);

void plotVortexFilaments( int n, Vec3d* CPs, Vec3d dDir ){
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0,1.0,1.0);
    for(int i=0; i<n;i++){
        glVertex3f( CPs[i].x, CPs[i].y, CPs[i].z );
    }
    glEnd();
    glColor3f(0.75,0.75,0.75);
    glBegin(GL_LINES);
    for(int i=0; i<n;i++){
        Vec3d p = CPs[i];
        glVertex3f( p.x, p.y, p.z );
        p.add( dDir );
        glVertex3f( p.x, p.y, p.z );
    }
    glEnd();

};

void plotVecPlane( Vec2i n, Vec3d p0, Vec3d a, Vec3d b, double sz, double dt, VecFieldFunc func ){
    glBegin(GL_LINES);
    for(int ia=0; ia<n.a; ia++ ){
        for(int ib=0; ib<n.b; ib++ ){
            Vec3d p = p0 + a*ia + b*ib;
            if( sz>0 ){
                glVertex3f( p.x-sz, p.y   , p.z    ); glVertex3f( p.x+sz, p.y,    p.z    );
                glVertex3f( p.x   , p.y-sz, p.z    ); glVertex3f( p.x,    p.y+sz, p.z    );
                glVertex3f( p.x   , p.y   , p.z-sz ); glVertex3f( p.x,    p.y,    p.z+sz );
            }
            glVertex3f( p.x, p.y, p.z );
            Vec3d v = func(p); //printf( "(%f,%f,%f) (%f,%f,%f)\n", p.x, p.y, p.z,  v.x, v.y, v.z );
            p.add_mul( v, dt);
            glVertex3f( p.x, p.y, p.z );
        }
    }
    glEnd();
}

void plotStreamLine( int n, double dt, Vec3d p, VecFieldFunc func ){
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<n; i++ ){
        Vec3d v = func(p);
        p.add_mul(v,dt);
        glVertex3f( p.x, p.y, p.z );
        //printf( "(%f,%f,%f) (%f,%f,%f)\n", p.x, p.y, p.z,  v.x, v.y, v.z );
    }
    glEnd();
    //exit(0);
}

void plotStreamLinePlane( Vec2i n, int ns, Vec3d p0, Vec3d a, Vec3d b, double dt, VecFieldFunc func ){
    for(int ia=0; ia<n.a; ia++ ){
        for(int ib=0; ib<n.b; ib++ ){
            Vec3d p = p0 + a*ia + b*ib;
            //printf( "(%i,%i) (%f,%f,%f) \n", ia, ib, p.x, p.y, p.z );
            plotStreamLine( ns, dt, p, func );
        }
    }
}

#endif
