
#ifndef TankHelpers_h
#define TankHelpers_h

void armorColorScale( float f ){
    float mf = 1 - f;
    //float r = f*f; float g = 4*mf*f; float b = mf*mf;

    float r = 1-mf*mf; float g = f*f*f; float b = mf*mf*mf*mf;
    //float r = 1-mf*mf; float g = f*f*f; float b = f*f*f*f*f*f + 0.5*mf*mf*mf*mf;

    //float renorm = 3/(r+g+b); r*=renorm; g*=renorm; b*=renorm;
    glColor3f( r, g, b );
}

void renderArmor( const VehicleBlock& block, double maxThickness ){
    int i=0;
    for( Polygon* pl : block.polygons ){
        //printf( " pl %i npoints %i \n", i, pl->ipoints.size() );
        armorColorScale( block.armor[i].thickness/maxThickness );
        Draw3D::drawPlanarPolygon( pl->ipoints.size(), &pl->ipoints.front(), &block.points.front() );
        //Vec3d c = block.faceCog( i );
        //glColor3f(1.0f,1.0f,1.0f); Draw3D::drawVecInPos( block.armor[i].normal, c );
        i++;
    }

}

void renderArmorCaptions( const VehicleBlock& block, float sz ){
    char str[64];
    int i=0;
    for( Polygon* pl : block.polygons ){
        sprintf(str,"%i:%3.0fmm%2.1fton\0",i+1, block.armor[i].thickness, block.armor[i].mass*1e-3 );
        Vec3d c = block.faceCog( i );
        Draw3D::drawText(str, c, fontTex, sz, 0 );
        i++;
    }
}

void drawTankWheels(Tank * tank){
    Vec3d gpos;
    for(int i=0; i<tank->nwheel; i++){
        tank->getWheelPos( i, gpos );
        //printf( "%i (%g,%g,%g)\n", i, gpos.x, gpos.y, gpos.z );
        Draw3D::drawPointCross( gpos, 0.5 );
    }
}

Terrain25D *  prepareTerrain(){
//    Terrain25D * terrain = new Terrain25D();

    Terrain25D_bicubic * terrain = new Terrain25D_bicubic();
    terrain->ruler.setup( (Vec2d){10.0,10.0}*-16, (Vec2d){10.0d,10.0d} );
    terrain->allocate( {32,32} );
    terrain->makeRandom( -2.0, 2.0 );

    terrain->shape = glGenLists(1);
    glNewList( terrain->shape , GL_COMPILE );
    int na=100,nb=100;
    float da=1.0,db=1.0;
    float x0=-0.5*da*na,y0=-0.5*db*nb;
    glEnable(GL_LIGHTING);
    glColor3f(0.5f,0.5f,0.5f);
    glNormal3f(0.0f,1.0f,0.0f);
    /*
    glBegin(GL_QUADS);
        glVertex3f( na*da,     0, 0 );
        glVertex3f( 0,         0, 0 );
        glVertex3f( 0,     nb*db, 0 );
        glVertex3f( na*da, nb*db, 0 );
    glEnd();
    */
    float * oldvals = new float[na*3];
    for(int ia=0; ia<na; ia++){
        glBegin(GL_TRIANGLE_STRIP);
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d dv1,dv2;
            Vec2d p1; p1.set( (ia  )*da+x0, ib*db+y0 );
            Vec2d p2; p2.set( (ia+1)*da+x0, ib*db+y0 );
            float v1,v2;
            if( ia == 0 ){
                v1 = (float)terrain->eval( p1, dv1 );
            }else{
                v1 = oldvals[i3]; dv1.x=oldvals[i3+1]; dv1.y=oldvals[i3+2];
            }
            v2 = (float)terrain->eval( p2, dv2 );
            oldvals[i3] = v2; oldvals[i3+1] = dv2.x; oldvals[i3+2] = dv2.y;
            //glNormal3f(-dv1.x,1.0,-dv1.y); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            //glNormal3f(-dv2.x,1.0,-dv2.y); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            glNormal3f(dv1.x,-1.0,dv1.y); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            glNormal3f(dv2.x,-1.0,dv2.y); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            //glColor3f(v1,0.5,-v1); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            //glColor3f(v2,0.5,-v2); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            //printf( " %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", p1.x, p1.y, v1 ,  p2.x, p2.y, v2  );
        }
        glEnd();
    }

    glBegin(GL_LINES);
    for(int ia=0; ia<na; ia++){
        for(int ib=0; ib<nb; ib++){
            int i3 = 3*ib;
            Vec2d p,dv; p.set( ia*da+x0, ib*db+y0 );
            double v = (float)terrain->eval( p, dv );
            glVertex3f( (float)p.x,         v, (float)p.y );
            glVertex3f( (float)(p.x-dv.x),  v+1.0, (float)(p.y-dv.y) );
        }

    }

    glEnd();
    glEndList();
    return terrain;
}

#endif
