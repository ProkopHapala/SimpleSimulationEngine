
#ifndef AeroDraw_h
#define AeroDraw_h

#include "AeroCraft.h"
#include "AeroTest.h"

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

void renderAeroCraft ( const AeroCraft& craft, bool bPos ){
	glPushMatrix();
        Mat3d camMat;
        camMat.set(craft.rotMat);
        float glmat[16];

        Vec3d pos;
        if   (bPos){ pos = craft.pos;            }
        else       { pos = (Vec3d){0.0,0.0,0.0}; }

        Draw3D::toGLMat( pos, camMat, glmat);
        //printf( "%g %g %g\n", pos.x, pos.y, pos.z);
        //printf( "%g %g %g\n", camMat.ax, camMat.ay, camMat.az);
        //printf( "%g %g %g\n", camMat.bx, camMat.by, camMat.bz);
        //printf( "%g %g %g\n", camMat.cx, camMat.cy, camMat.cz);
        glMultMatrixf( glmat );
        //glDisable (GL_LIGHTING);
        //glColor3f( 0.3f,0.3f,0.3f );

        Draw3D::drawSphere_oct(1,0.2,{0.0,0.0,0.0});

        //Draw3D::drawLine( {0,0,0},wingLeft .lpos);
        //Draw3D::drawLine( {0,0,0},wingRight.lpos);
        //Draw3D::drawLine( {0,0,0},elevator .lpos);
        //wingLeft .render();
        //wingRight.render();
        //rudder   .render();
        //elevator .render();


		for( int i=0; i<craft.nPanels; i++ ){
            Draw3D::drawLine( {0,0,0}, craft.panels[i].lpos);
            //printf( "%g %g %g\n", panels[i].lpos.x, panels[i].lpos.y, panels[i].lpos.z);
            //panels[i].render( );
            //Draw3D::drawKite( craft.panels[i].lpos, craft.panels[i].lrot, sqrt(craft.panels[i].area) );
            double sz = sqrt(craft.panels[i].area);
            Draw3D:: drawPanel( craft.panels[i].lpos, craft.panels[i].lrot, {sz,sz*0.25} );
		}

	glPopMatrix();
};

void renderAeroCraftDebug ( const AeroCraft& craft, bool bPos ){
	glPushMatrix();
        Mat3d camMat;
        camMat.setT(craft.rotMat);
        float glmat[16];

        Vec3d pos;
        if   (bPos){ pos = craft.pos;            }
        else       { pos = (Vec3d){0.0,0.0,0.0}; }

        Draw3D::toGLMat( pos, camMat, glmat);
        //printf( "%g %g %g\n", pos.x, pos.y, pos.z);
        //printf( "%g %g %g\n", camMat.ax, camMat.ay, camMat.az);
        //printf( "%g %g %g\n", camMat.bx, camMat.by, camMat.bz);
        //printf( "%g %g %g\n", camMat.cx, camMat.cy, camMat.cz);
        glMultMatrixf( glmat );
        //glDisable (GL_LIGHTING);
        //glColor3f( 0.3f,0.3f,0.3f );

        Draw3D::drawSphere_oct(1,0.2,{0.0,0.0,0.0});

        //Draw3D::drawLine( {0,0,0},wingLeft .lpos);
        //Draw3D::drawLine( {0,0,0},wingRight.lpos);
        //Draw3D::drawLine( {0,0,0},elevator .lpos);
        //wingLeft .render();
        //wingRight.render();
        //rudder   .render();
        //elevator .render();


		for( int i=0; i<craft.nPanels; i++ ){
            Draw3D::drawLine( {0,0,0}, craft.panels[i].lpos);
            //printf( "%g %g %g\n", panels[i].lpos.x, panels[i].lpos.y, panels[i].lpos.z);
            //panels[i].render( );
            //Draw3D::drawKite( craft.panels[i].lpos, craft.panels[i].lrot, sqrt(craft.panels[i].area) );
            double sz = sqrt(craft.panels[i].area);
            Draw3D:: drawPanel( craft.panels[i].lpos, craft.panels[i].lrot, {sz,sz*0.25} );
		}

	glPopMatrix();
};




Terrain25D * prepareTerrain( int nsz, int nsub, double step, double hmax ){
    Terrain25D_bicubic * terrain = new Terrain25D_bicubic();
    terrain->ruler.setup( (Vec2d){nsz*0.5,nsz*0.5}*-step, (Vec2d){step,step} );
    terrain->allocate( {nsz,nsz} );
    terrain->makeRandom( 0.0, hmax );

    terrain->shape = glGenLists(1);
    glNewList( terrain->shape , GL_COMPILE );
    //int na=100,nb=100;

    int na = (terrain->ruler.n.a - 3)*nsub;
    int nb = (terrain->ruler.n.b - 3)*nsub;
    float da=terrain->ruler.step.a/float(nsub);
    float db=terrain->ruler.step.b/float(nsub);
    float x0=terrain->ruler.pos0.x;
    float y0=terrain->ruler.pos0.y;

    glEnable(GL_LIGHTING);
    glColor3f (0.5f,0.5f,0.5f);
    glNormal3f(0.0f,1.0f,0.0f);

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
            glNormal3f(-dv1.x,1.0,-dv1.y); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            glNormal3f(-dv2.x,1.0,-dv2.y); glVertex3f( (float)p2.x,  v2, (float)p2.y );

            //glColor3f(v1,0.5,-v1); glVertex3f( (float)p1.x,  v1, (float)p1.y );
            //glColor3f(v2,0.5,-v2); glVertex3f( (float)p2.x,  v2, (float)p2.y );
            //printf( " %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f)\n", p1.x, p1.y, v1 ,  p2.x, p2.y, v2  );
        }
        glEnd();
    }

    // Normals
    /*
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
    */
    glEndList();
    delete [] oldvals;
    return terrain;
}

void renderSkyBox( float x0, float y0, float z0, float skysz ){
	glDepthMask(0);
	glDisable (GL_LIGHTING);
	glDisable (GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	//float skysz = VIEW_DEPTH*0.25;
	float R0=0.1,G0=0.1,B0=0.5;
	float R1=0.7,G1=0.8,B1=0.8;
	glBegin(GL_QUADS);
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, -skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, -skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, -skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, -skysz+z0 );

		glColor3f( R1, G1, B1 );  glVertex3f( -skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R1, G1, B1 );  glVertex3f( +skysz+x0,     0+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( +skysz+x0, skysz+y0, +skysz+z0 );
		glColor3f( R0, G0, B0 );  glVertex3f( -skysz+x0, skysz+y0, +skysz+z0 );
	glEnd();
	glDepthMask(1);
}

int makeBuildingsGrid( int nx, int ny, float sx, float sy, float cx, float cy,  float min_height, float max_height ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
	for (int ix=-nx; ix<nx; ix++){
		float x = ix*sx;
		for (int iy=-ny; iy<ny; iy++){
			float height = randf() * (max_height-min_height) + min_height;
			float y = iy*sy;
			Draw3D::drawBox( x, x + sx*cx, 0, height, y, y + sy*cy, 0.75f, 0.75f, 0.75f );
		}
	}
	glEndList();
	return( ilist );
}

int makeBuildingsClusters( int nclustest, int nmin, int nmax, float minx, float maxx, float miny, float maxy, float min_dist, float max_dist, float min_size, float max_size, float min_height, float max_height ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
	int nboxes = 0;
	for (int icluster=0; icluster<nclustest; icluster++){
		float x0 = randf()*(maxx-minx) + minx;
		float y0 = randf()*(maxy-miny) + miny;
		float nb = round(randf()*(nmax - nmin)) + nmin;
		for (int ib=0; ib<nb; ib++){
			float height = randf() * (max_height-min_height) + min_height;
			float x  = x0 + randf()*(max_dist-min_dist) + min_dist;
			float y  = y0 + randf()*(max_dist-min_dist) + min_dist;
			float dx = randf()*(max_size-min_size) + min_size;
			float dy = randf()*(max_size-min_size) + min_size;
			Draw3D::drawBox( x-dx, x+dx, 0, height, y-dy, y+dy, 0.75f, 0.75f, 0.75f );
			nboxes++;
		};
	};
	printf(" %i buildings \n", nboxes );
	glEndList();
	return( ilist );
}

void drawStaticTest2D( const AeroTester& data, int fontTex, float WIDTH, float HEIGHT ){
    char str[256];
    /*
        Vec3d pos   = data->trjPos  [i];
        Vec3d vel   = data->trjVel  [i];
        Vec3d Force = data->trjForce[i];
        Vec3d Fw    = data->trjFw   [i];
        Vec3d Up    = data->trjUp   [i];
    */

    int n = data.ntrj;
     double t;
    Vec3d pos,vel,force,Up,Fw;
    // --- attitude
    double attitude;
    //glPushMatrix();
    //glScalef(1.0,0.5,1.0);
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f,0.0f,0.0f);
    for(int i=0; i<n; i++){
        t   = data.trjT    [i];
        pos = data.trjPos[i];
        glVertex3f(t,pos.y, 1.0);
    }
    glEnd();
    glColor4f(1.0f,1.0f,1.0f,0.9f);
    t        = data.trjT  [0]; attitude = data.trjPos[0].y;
    sprintf(str, "attitude=%4.3f m\0", attitude );
    Draw2D::drawText( str, 0, {t,attitude},0.0,fontTex, 10 );
    t        = data.trjT[n-1]; attitude = data.trjPos[n-1].y;
    sprintf(str, "attitude=%4.3f m\0", attitude );
    Draw2D::drawText( str, 0, {t,attitude},0.0,fontTex, 10 );
    //glPopMatrix();

    // --- pos.zy
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f,0.8f,0.8f);
    for(int i=0; i<n; i++){
        pos = data.trjPos[i];
        glVertex3f(pos.z*0.02+WIDTH*0.5,pos.y*0.02+HEIGHT*0.5, 1.0);
    }
    glEnd();

    // --- pos.xy
    glBegin(GL_LINE_STRIP);
    glColor3f(0.8f,0.8f,0.0f);
    for(int i=0; i<n; i++){
        pos = data.trjPos[i];
        glVertex3f(pos.x*0.02+WIDTH*0.5,pos.y*0.02+HEIGHT*0.5, 1.0);
    }
    glEnd();

    // --- pos.xz
    int iturn_old=0;
    int iturn    =0;
    double xold  =0.0;
    double xmin=1e+8,xmax=-1e+8,zmin=1e+8,zmax=-1e+8;
    glBegin(GL_LINE_STRIP);
    glColor3f(0.8f,0.8f,0.8f);
    for(int i=0; i<n; i++){
        pos = data.trjPos[i];
        vel = data.trjVel[i];
        if(pos.x<xmin) xmin=pos.x;
        if(pos.x>xmax) xmax=pos.x;
        if(pos.z<zmin) zmin=pos.z;
        if(pos.z>zmax) zmax=pos.z;
        if(vel.x*xold <0 ){ iturn_old=iturn; iturn=i; } xold=vel.x;
        glVertex3f(pos.x*0.02+WIDTH*0.5,pos.z*0.02+HEIGHT*0.5, 1.0);
    }
    glEnd();
    //sprintf(str, "xspan=%4.3f[m] zspan=%4.3f[m] T=%4.3f [s]\0", xmax-xmin, zmax-zmin, data->trjT[iturn]-data->trjT[iturn_old] );
    sprintf(str, "xspan=%4.3f[m] T=%4.3f [s]\0", xmax-xmin, data.trjT[iturn]-data.trjT[iturn_old] );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawText( str, 0, {WIDTH*0.5,HEIGHT*0.5},0.0,fontTex, 10 );
	pos=data.trjPos[iturn    ]; Draw2D::drawPointCross({pos.x*0.02+WIDTH*0.5,pos.z*0.02+HEIGHT*0.5},10);
	pos=data.trjPos[iturn_old]; Draw2D::drawPointCross({pos.x*0.02+WIDTH*0.5,pos.z*0.02+HEIGHT*0.5},10);


    // --- v.y
    float d0 = HEIGHT*0.5f;
    glBegin(GL_LINE_STRIP);
    glColor3f(1.0f,0.0f,1.0f);
    for(int i=0; i<n; i++){
        t     = data.trjT    [i];
        vel   = data.trjVel  [i];
        glVertex3f(t,vel.y+d0, 1.0);
    }
    glEnd();
    glBegin(GL_LINES);
        glVertex3f(0,d0, 1.0); glVertex3f(WIDTH,d0, 1.0);
    glEnd();


    // --- speed
    glBegin(GL_LINE_STRIP);
    glColor3f(0.0f,0.0f,1.0f);
    double speed;
    for(int i=0; i<n; i++){
        t     = data.trjT    [i];
        vel   = data.trjVel  [i];
        speed = vel.norm();
        glVertex3f(t,speed, 1.0);
    }
    glEnd();
    sprintf(str, "v=%4.3fm/s(%4.3fkm/h)\0", speed, speed*3.6 );
	glColor4f(1.0f,1.0f,1.0f,0.9f); Draw2D::drawText( str, 0, {t,speed},0.0,fontTex, 10 );

}

#endif
