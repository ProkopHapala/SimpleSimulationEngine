
//#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
//#include <GLU.h>

#include "Vec2.h"
#include "Draw.h"

#include "Draw3D.h" // THE HEADER


namespace Draw3D{

//void vertex(const Vec3f& v ){ glVertex3f(v.x,v.y,v.z); }
//void vertex(const Vec3d& v ){ glVertex3f(v.x,v.y,v.z); }
//void color (const Vec3f& v ){ glColor3f (v.x,v.y,v.z); }
//void color (const Vec3d& v ){ glColor3f (v.x,v.y,v.z); }
//void normal(const Vec3f& v ){ glNormal3f(v.x,v.y,v.z); }
//void normal(const Vec3d& v ){ glNormal3f(v.x,v.y,v.z); }

void drawPoint( const Vec3f& vec ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_POINTS);
		glVertex3d( vec.x, vec.y, vec.z );
	glEnd();
};

void drawPointCross_bare( const Vec3f& vec, float sz ){
    glVertex3f( vec.x-sz, vec.y, vec.z ); glVertex3f( vec.x+sz, vec.y, vec.z );
    glVertex3f( vec.x, vec.y-sz, vec.z ); glVertex3f( vec.x, vec.y+sz, vec.z );
    glVertex3f( vec.x, vec.y, vec.z-sz ); glVertex3f( vec.x, vec.y, vec.z+sz );
}

void drawPointCross( const Vec3f& vec, float sz ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
	drawPointCross_bare(vec,sz );
    glEnd();
};

void drawPointCross( const Vec3f& vec, double sz ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3d( vec.x-sz, vec.y, vec.z ); glVertex3d( vec.x+sz, vec.y, vec.z );
		glVertex3d( vec.x, vec.y-sz, vec.z ); glVertex3d( vec.x, vec.y+sz, vec.z );
		glVertex3d( vec.x, vec.y, vec.z-sz ); glVertex3d( vec.x, vec.y, vec.z+sz );
	glEnd();
};

void drawVec( const Vec3f& vec ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3d( 0, 0, 0 ); glVertex3d( vec.x, vec.y, vec.z );
	glEnd();
};

void drawVecInPos( const Vec3f& v, const Vec3f& pos ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3d( pos.x, pos.y, pos.z ); glVertex3d( pos.x+v.x, pos.y+v.y, pos.z+v.z );
	glEnd();
};

void drawLine( const Vec3f& p1, const Vec3f& p2 ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glVertex3d( p1.x, p1.y, p1.z ); glVertex3d( p2.x, p2.y, p2.z );
	glEnd();
};

void drawArrow( const Vec3f& p1, const Vec3f& p2, float sz ){
	//glDisable (GL_LIGHTING);
	Vec3f up,lf,p;
    Vec3f fw = p2-p1; fw.normalize();
    fw.getSomeOrtho(up,lf);
    fw.mul(sz); lf.mul(sz); up.mul(sz);
	glBegin   (GL_LINES);
		glVertex3d( p1.x, p1.y, p1.z ); glVertex3d( p2.x, p2.y, p2.z );
		p = p2 - fw + up; glVertex3d( p.x, p.y, p.z ); glVertex3d( p2.x, p2.y, p2.z );
        p = p2 - fw - up; glVertex3d( p.x, p.y, p.z ); glVertex3d( p2.x, p2.y, p2.z );
        p = p2 - fw + lf; glVertex3d( p.x, p.y, p.z ); glVertex3d( p2.x, p2.y, p2.z );
        p = p2 - fw - lf; glVertex3d( p.x, p.y, p.z ); glVertex3d( p2.x, p2.y, p2.z );
	glEnd();
};

void vecsInPoss( int n, const Vec3d* vs, const Vec3d* ps, float sc ){
    //printf("%i %i\n", n, closed );
    glBegin(GL_LINES);
    for(int i=0; i<n; i++){
        //printf("%i (%3.3f,%3.3f,%3.3f)\n", i, ps[i].x, ps[i].y, ps[i].z );
        glVertex3d( ps[i].x,            ps[i].y,            ps[i].z            );
        glVertex3d( ps[i].x+vs[i].x*sc, ps[i].y+vs[i].y*sc, ps[i].z+vs[i].z*sc );
    };
    glEnd();
};

void drawPolyLine( int n, Vec3d * ps, bool closed ){   // closed=false
    //printf("%i %i\n", n, closed );
    if(closed){ glBegin(GL_LINE_LOOP); }else{ glBegin(GL_LINE_STRIP); }
    for(int i=0; i<n; i++){
        //printf("%i (%3.3f,%3.3f,%3.3f)\n", i, ps[i].x, ps[i].y, ps[i].z );
        glVertex3d( ps[i].x, ps[i].y, ps[i].z );
    };
    glEnd();
};

void drawScale( const Vec3f& p1, const Vec3f& p2, const Vec3f& up, float tick, float sza, float szb ){
	//glDisable (GL_LIGHTING);
	Vec3f d,a,b,p;
	d.set_sub( p2, p1 );
	float L = d.norm();
	int n = (L+0.0001)/tick;
	d.mul( 1/L );
	a.set(up);
	a.add_mul( d, -d.dot( a ) );
	b.set_cross( d, a  );
	glBegin   (GL_LINES);
	p.set(p1);
	d.mul( L/n );
	a.mul(sza);
	b.mul(szb);
    for( int i=0; i<=n; i++){
        glVertex3d( p.x-a.x, p.y-a.y, p.z-a.z );
        glVertex3d( p.x+a.x, p.y+a.y, p.z+a.z );
        glVertex3d( p.x-b.x, p.y-b.y, p.z-a.z );
        glVertex3d( p.x+b.x, p.y+b.y, p.z+b.z );
        p.add( d );
    }
    glVertex3d( p1.x, p1.y, p1.z ); glVertex3d( p2.x, p2.y, p2.z );
	glEnd();
};

void drawTriangle_bare( const Vec3f& p1, const Vec3f& p2, const Vec3f& p3 ){
    //glNormal3d( normal.x, normal.y, normal.z );
    Vec3f d1,d2,nr;
	d1.set( p2 - p1 );
	d2.set( p3 - p1 );
	nr.set_cross(d1,d2);
	nr.normalize();
	glNormal3d( nr.x, nr.y, nr.z );
    glVertex3d( p1.x, p1.y, p1.z );
    glVertex3d( p2.x, p2.y, p2.z );
    glVertex3d( p3.x, p3.y, p3.z );
    /*
    glColor3f(0.0,0.0,0.0);
    glVertex3d( p1.x, p1.y, p1.z ); glVertex3d( p2.x, p2.y, p2.z );
    glVertex3d( p2.x, p2.y, p2.z ); glVertex3d( p3.x, p3.y, p3.z );
    glVertex3d( p3.x, p3.y, p3.z ); glVertex3d( p1.x, p1.y, p1.z );
    glColor3f(1.0,1.0,1.0);
    nr.mul(0.1);
    glVertex3d( p1.x, p1.y, p1.z ); glVertex3d( p1.x+nr.x, p1.y+nr.y, p1.z+nr.z );
    glVertex3d( p2.x, p2.y, p2.z ); glVertex3d( p2.x+nr.x, p2.y+nr.y, p2.z+nr.z );
    glVertex3d( p3.x, p3.y, p3.z ); glVertex3d( p3.x+nr.x, p3.y+nr.y, p3.z+nr.z );
    */
}

void drawTriangle( const Vec3f& p1, const Vec3f& p2, const Vec3f& p3 ){
	glBegin   (GL_TRIANGLES);
        drawTriangle_bare( p1, p2, p3 );
	glEnd();
}

void drawTriangle ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3, bool filled ){
    int primitive;
    if(filled){ primitive=GL_TRIANGLE_FAN; }else{ primitive=GL_LINE_LOOP; }
    glBegin(primitive);
    vertex(p1); vertex(p2); vertex(p3);
    glEnd();
}


void drawQuad_bare( const Vec3f& p1, const Vec3f& p2, const Vec3f& p3, const Vec3f& p4 ){
    //glNormal3d( normal.x, normal.y, normal.z );
    double r13=(p3-p1).norm2();
    double r24=(p4-p2).norm2();
    if(r13>r24){ drawTriangle_bare( p1, p2, p4 ); drawTriangle_bare( p2, p3, p4 ); }
    else       { drawTriangle_bare( p1, p2, p3 ); drawTriangle_bare( p3, p4, p1 ); }
}


void drawQuad     ( const Vec3f& p1,  const Vec3f& p2, const Vec3f& p3, const Vec3f& p4, bool filled ){
    int primitive;
    if(filled){ primitive=GL_TRIANGLE_FAN; }else{ primitive=GL_LINE_LOOP; }
    glBegin(primitive);
    //printf( " drawQuad (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) (%g,%g,%g) \n", p1.x,p1.y,p1.z,  p2.x,p2.y,p2.z,  p3.x,p3.y,p3.z, p4.x,p4.y,p4.z  );
    vertex(p1); vertex(p2); vertex(p3); vertex(p4);
    glEnd();
}


Vec3f lincomb(const Vec3f& p1, const Vec3f& p2, double v1, double v2 ){
    double f  = v1/(v1-v2);
    Vec3f p;
    p.set_lincomb( 1-f, p1, f, p2 );
    return p;
}

void drawTetraIso( Vec3f** ps, Quat4d vals ){
    bool b0 = vals.x>0;
    bool b1 = vals.y>0;
    bool b2 = vals.z>0;
    bool b3 = vals.w>0;
    int n01 = b0+b1;
    int n23 = b2+b3;
    int n   = n01+n23;
    if( n==0 || n==4 ) return;
    int i1,i2;
    int j1,j2,j3;
    if(n==1){        // triangle
        if(n01){
            if(b0){ i1=0; j1=1; j2=2; j3=3; } // b1
            else  { i1=1; j1=0; j2=3; j3=2; } // b2
        }else{
            if(b2){ i1=2; j1=3; j2=0; j3=1; } // b3
            else  { i1=3; j1=2; j2=1; j3=0; } // b4
        }
        drawTriangle_bare(
            lincomb(*ps[i1],*ps[j1],vals.array[i1],vals.array[j1]),
            lincomb(*ps[i1],*ps[j2],vals.array[i1],vals.array[j2]),
            lincomb(*ps[i1],*ps[j3],vals.array[i1],vals.array[j3])
        );
    }else if (n==3){ // triangle
        if(n01==1){
            if(!b0){ i1=0; j1=1; j2=3; j3=2; } // b1
            else   { i1=1; j1=0; j2=2; j3=3; } // b2
        }else{
            if(!b2){ i1=2; j1=1; j2=0; j3=3; } // b3
            else   { i1=3; j1=1; j2=2; j3=0; } // b4
        }
        //printf("n=3 %i | %i %i %i   (%i%i%i%i)\n", i1, j1, j2, j3, b0,b1,b2,b3 );
        drawTriangle_bare(
            lincomb(*ps[i1],*ps[j1],-vals.array[i1],-vals.array[j1]),
            lincomb(*ps[i1],*ps[j2],-vals.array[i1],-vals.array[j2]),
            lincomb(*ps[i1],*ps[j3],-vals.array[i1],-vals.array[j3])
        );
    }else if (n==2){ // quad
        if(n01==1){
            if(b0){
                if(b2){ i1=0; i2=2; j1=3; j2=1; } // b1-b3 | b2-b4
                else  { i1=0; i2=3; j1=1; j2=2; } // b1-b4 | b2-b3
            }else{
                if(b3){ j1=0; j2=2; i1=3; i2=1; } // b2-b3 | b1-b4
                else  { j1=0; j2=3; i1=1; i2=2; } // b2-b4 | b1-b3
            }
        }else{
            if(n01==2){ i1=0; i2=1; j1=2; j2=3; }
            else      { j1=0; j2=1; i1=2; i2=3; }
        }
        //printf("n=%i  %i,%i|%i,%i   (%i%i%i%i)\n", n, i1, i2, j1, j2, b0,b1,b2,b3 )
        drawQuad_bare(
            lincomb(*ps[i1],*ps[j1],vals.array[i1],vals.array[j1]),
            lincomb(*ps[i1],*ps[j2],vals.array[i1],vals.array[j2]),
            lincomb(*ps[i2],*ps[j2],vals.array[i2],vals.array[j2]),
            lincomb(*ps[i2],*ps[j1],vals.array[i2],vals.array[j1])
        );
    }
};

void drawSimplexLines( Vec3f** ps ){
    vertex( *ps[0] ); vertex( *ps[1] );
    vertex( *ps[0] ); vertex( *ps[2] );
    vertex( *ps[0] ); vertex( *ps[3] );
    vertex( *ps[1] ); vertex( *ps[2] );
    vertex( *ps[1] ); vertex( *ps[3] );
    vertex( *ps[2] ); vertex( *ps[3] );
}

void drawMatInPos( const Mat3f& mat, const Vec3f& pos, const Vec3f& sc ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glColor3f( 1, 0, 0 ); glVertex3d( pos.x, pos.y, pos.z ); glVertex3d( pos.x+mat.xx*sc.x, pos.y+mat.xy*sc.x, pos.z+mat.xz*sc.x );
		glColor3f( 0, 1, 0 ); glVertex3d( pos.x, pos.y, pos.z ); glVertex3d( pos.x+mat.yx*sc.y, pos.y+mat.yy*sc.y, pos.z+mat.yz*sc.y );
		glColor3f( 0, 0, 1 ); glVertex3d( pos.x, pos.y, pos.z ); glVertex3d( pos.x+mat.zx*sc.z, pos.y+mat.zy*sc.z, pos.z+mat.zz*sc.z );
	glEnd();
};

void drawShape( int shape, const Vec3f& pos, const Mat3f& rot, bool trasposed ){
	glPushMatrix();
	float glMat[16];
	if( trasposed ){
        toGLMatT ( pos, rot, glMat );
	}else{
        toGLMat( pos, rot, glMat );
	}
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};

void drawShape    ( int shape, const Vec3f& pos, const Quat4f& qrot ){
	glPushMatrix();
	float glMat[16];
	toGLMat ( pos, qrot, glMat );
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};

void drawShape    ( int shape, const Vec3f& pos, const Quat4f& qrot, const Vec3f& scale ){
	glPushMatrix();
	float glMat[16];
	toGLMat ( pos, qrot, scale, glMat );
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};

void shapeInPoss(  int shape, int n, const Vec3d* pos, const double* sizes, const Mat3d& rot, bool transposed ){
    Mat3f mat = (Mat3f)rot;
    for(int i=0; i<n; i++){
        if(sizes) mat.mul((float)sizes[i]);
        drawShape( shape, (Vec3f)pos[i], mat, transposed );
    }
}

/*
void drawShapeT( const Vec3f& pos, const Mat3f& rot, int shape ){
	glPushMatrix();
	float glMat[16];
	toGLMatT( pos, rot, glMat );
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};

void drawShapeT( const Vec3d& pos, const Mat3d& rot, int shape ){
	glPushMatrix();
	float glMat[16];
	toGLMatT( pos, rot, glMat );
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};
void drawShapeT    ( const Vec3d& pos, const Quat4d& qrot, int shape ){
	glPushMatrix();
	float glMat[16];
	Quat4d qrotT;
	qrotT.setInverseUnitary(qrot);
	toGLMat( pos, qrotT, glMat );
	glMultMatrixf( glMat );
	glCallList( shape );
	glPopMatrix();
};
*/

int drawConeFan( int n, float r, const Vec3f& base, const Vec3f& tip ){
	int nvert=0;
	Vec3f a,b,c,c_hat;
	c.set_sub( tip, base );
	c_hat.set_mul( c, 1/c.norm() );
	c_hat.getSomeOrtho( a, b );
	a.normalize();
	b.normalize();
    //float alfa = 2*M_PI/n;
    float alfa = 2*M_PI/n;
    Vec2f rot,drot;
    rot .set(1.0f,0.0f);
    drot.set( cos( alfa ), sin( alfa ) );

	Vec3f q; q.set(c); q.add_mul( a, -r );
	float pnab =  c_hat.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	glBegin   ( GL_TRIANGLE_FAN );
	//glBegin   ( GL_LINES );
	//glBegin   ( GL_LINE_STRIP );

	glNormal3f( c_hat.x, c_hat.z, c_hat.z );
	//printf( "pn0 %f %f %f \n", c_hat.x, c_hat.z, c_hat.z );
	glVertex3f( tip.x, tip.y, tip.z ); nvert++;
	for(int i=0; i<=n; i++ ){
        Vec3f p,pn;
		p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
		pn.set( pnab*p.x + pnc*c_hat.x, pnab*p.y + pnc*c_hat.y, pnab*p.z + pnc*c_hat.z  );

        glNormal3f( pn.x, pn.y, pn.z );
		glVertex3f( base.x + r*p.x, base.y + r*p.y, base.z + r*p.z ); nvert++;
        rot.mul_cmplx( drot );
	}
	glEnd();
	return nvert;
};

int drawCylinderStrip( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip ){
	int nvert=0;

	Vec3f a,b,c,c_hat;
	c.set_sub( tip, base );
	c_hat.set_mul( c, 1/c.norm() );
	c_hat.getSomeOrtho( a, b );
	a.normalize();
	b.normalize();

    float alfa = 2*M_PI/n;
    Vec2f rot,drot;
    rot .set(1.0f,0.0f);
    drot.set( cos( alfa ), sin( alfa ) );

	Vec3f q; q.set(c); q.add_mul( a, -(r1-r2) );
	float pnab =  c_hat.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	glBegin   ( GL_TRIANGLE_STRIP );
	//glBegin   ( GL_LINES );
	for(int i=0; i<=n; i++ ){
		Vec3f p,pn;
		p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
		pn.set( pnab*p.x + pnc*c_hat.x, pnab*p.y + pnc*c_hat.y, pnab*p.z + pnc*c_hat.z  );
		//printf( "p %f %f %f   pn %f %f %f |pn| %f \n", p.x, p.y, p.z,   pn.x, pn.y, pn.z, pn.norm() );
		glNormal3f( pn.x, pn.y, pn.z );
		glVertex3f( base.x + r1*p.x, base.y + r1*p.y, base.z + r1*p.z ); nvert++;
		glVertex3f( tip .x + r2*p.x, tip .y + r2*p.y, tip .z + r2*p.z ); nvert++;
        rot.mul_cmplx( drot );
	}
	glEnd();
	return nvert;
};

int drawCylinderStrip_wire( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip ){
	int nvert=0;

	Vec3f a,b,c,c_hat;
	c.set_sub( tip, base );
	c_hat.set_mul( c, 1/c.norm() );
	c_hat.getSomeOrtho( a, b );
	a.normalize();
	b.normalize();

    float alfa = 2*M_PI/n;
    Vec2f rot,drot;
    rot .set(1.0f,0.0f);
    drot.set( cos( alfa ), sin( alfa ) );

	Vec3f q; q.set(c); q.add_mul( a, -(r1-r2) );
	float pnab =  c_hat.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	glBegin   ( GL_LINE_LOOP );
	//glBegin   ( GL_LINES );
	for(int i=0; i<n; i++ ){
		Vec3f p;
		p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
		glVertex3f( base.x + r1*p.x, base.y + r1*p.y, base.z + r1*p.z ); nvert++;
		glVertex3f( tip .x + r2*p.x, tip .y + r2*p.y, tip .z + r2*p.z ); nvert++;
        rot.mul_cmplx( drot );
	}
    for(int i=0; i<n; i++ ){
		Vec3f p;
		p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
		glVertex3f( base.x + r1*p.x, base.y + r1*p.y, base.z + r1*p.z ); nvert++;
        rot.mul_cmplx( drot );
	}
    for(int i=0; i<n; i++ ){
		Vec3f p;
		p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
		glVertex3f( tip .x + r2*p.x, tip .y + r2*p.y, tip .z + r2*p.z ); nvert++;
        rot.mul_cmplx( drot );
	}
	glEnd();
	return nvert;
};

int drawCone( int n, float phi0, float phi1, float r1, float r2, const Vec3f& base, const Vec3f& tip, bool smooth ){
	int nvert=0;

	Vec3f a,b,c,c_hat;
	c.set_sub( tip, base );
	c_hat.set_mul( c, 1/c.norm() );
	c_hat.getSomeOrtho( a, b );
	a.normalize();
	b.normalize();

    //float alfa = 2*M_PI/n;
    float alfa = (phi1-phi0)/n;
    Vec2f rot,drot;
    //rot .set(1.0f,0.0f);
    rot.set( cos( phi0 ), sin( phi0 ) );
    drot.set( cos( alfa ), sin( alfa ) );

	Vec3f q; q.set(c); q.add_mul( a, -(r1-r2) );
	float pnab =  c_hat.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	glBegin   ( GL_QUADS );
	//glBegin   ( GL_LINES );
	Vec3f p,pn,op,opn;
    op .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
	opn.set( pnab*op.x + pnc*c_hat.x, pnab*op.y + pnc*c_hat.y, pnab*op.z + pnc*c_hat.z  );
	if( smooth ){
        for(int i=0; i<n; i++ ){

            glNormal3f( opn.x, opn.y, opn.z );		glVertex3f( base.x + r1*op.x, base.y + r1*op.y, base.z + r1*op.z ); nvert++;
            glNormal3f( opn.x, opn.y, opn.z );		glVertex3f( tip .x + r2*op.x, tip .y + r2*op.y, tip .z + r2*op.z ); nvert++;

            rot.mul_cmplx( drot );
            p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
            pn.set( pnab*p.x + pnc*c_hat.x, pnab*p.y + pnc*c_hat.y, pnab*p.z + pnc*c_hat.z  );
            pn.normalize();

            glNormal3f( pn.x, pn.y, pn.z );
            glVertex3f( tip .x + r2*p.x, tip .y + r2*p.y, tip .z + r2*p.z ); nvert++;
            glVertex3f( base.x + r1*p.x, base.y + r1*p.y, base.z + r1*p.z ); nvert++;

            op.set(p);
            opn.set(pn);

        }
	}else{
        for(int i=0; i<n; i++ ){

            //printf( " %i (%3.3f,%3.3f) \n", i, rot.x, rot.y );

            rot.mul_cmplx( drot );

            p .set( rot.x*a.x +  rot.y*b.x, rot.x*a.y + rot.y*b.y, rot.x*a.z + rot.y*b.z    );
            pn.set( pnab*p.x + pnc*c_hat.x, pnab*p.y + pnc*c_hat.y, pnab*p.z + pnc*c_hat.z  );

            Vec3f normal; normal.set_add( opn, pn ); normal.normalize();

            glNormal3f( normal.x, normal.y, normal.z );
            glVertex3f( base.x + r1*op.x, base.y + r1*op.y, base.z + r1*op.z ); nvert++;
            glVertex3f( tip .x + r2*op.x, tip .y + r2*op.y, tip .z + r2*op.z ); nvert++;

            glVertex3f( tip .x + r2*p.x, tip .y + r2*p.y, tip .z + r2*p.z ); nvert++;
            glVertex3f( base.x + r1*p.x, base.y + r1*p.y, base.z + r1*p.z ); nvert++;

            op.set(p);
            opn.set(pn);

        }
	}
	glEnd();
	return nvert;
};

int drawSphereTriangle( int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c ){
	int nvert=0;
	float d = 1.0f/n;
	Vec3f da,db;
	da.set_sub( a, c ); da.mul( d );
	db.set_sub( b, c ); db.mul( d );
	for( int ia=0; ia<n; ia++ ){
		Vec3f p0,p; p0.set( c );
		p0.add_mul( da, ia );
		p.set_mul( p0, 1.0f/p0.norm() );
		glBegin   (GL_TRIANGLE_STRIP);
		//glBegin   (GL_LINES);
		//glColor3f( d*ia, 0, 0 );
		glNormal3f( p.x, p.y, p.z );
		glVertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
		//glVertex3f( r*p.x+pos.x+p.x, r*p.y+pos.y+p.y, r*p.z+pos.z+p.z );
		for( int ib=0; ib<(n-ia); ib++ ){
			Vec3f p;
			p.set_add( p0, da );
			p.normalize();
			//glColor3f( 0, 1, 0 );
			glNormal3f( p.x, p.y, p.z );
			glVertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
			//glVertex3f( r*p.x+pos.x+p.x, r*p.y+pos.y+p.y, r*p.z+pos.z+p.z );
			p.set_add( p0, db );
			p.normalize();
			//glColor3f( 0, 0, 1 );
			glNormal3f( p.x, p.y, p.z );
			glVertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
			//glVertex3f( r*p.x+pos.x+p.x, r*p.y+pos.y+p.y, r*p.z+pos.z+p.z );
			p0.add( db );
			//printf(" %f %f %f %f \n", p.x, p.y, p.z, p.norm() );
		}
		glEnd();
	}
	return nvert;
};

int drawSphereTriangle_wire( int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c ){
	int nvert=0;
	float d = 1.0f/n;
	Vec3f da,db;
	da.set_sub( a, c ); da.mul( d );
	db.set_sub( b, c ); db.mul( d );
	for( int ia=0; ia<n; ia++ ){
		Vec3f p0,p; p0.set( c );
		p0.add_mul( da, ia );
        glBegin   (GL_LINE_STRIP); //glColor3f(0.0,0.0,1.0);
        p.set(p0); p.normalize();
        glVertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
		for( int ib=0; ib<(n-ia); ib++ ){
			p.set_add( p0, da ); p.normalize();
			glVertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
			p.set_add( p0, db ); p.normalize();
			glVertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
			p0.add( db );
		}
        glEnd();
		glBegin   (GL_LINE_STRIP);         //glColor3f(1.0,0.0,0.0);
        for( int ib=0; ib<=(n-ia); ib++ ){
			//p.set_add( p0, da );
			p.set(p0); p.normalize();
			glVertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
			p0.sub( db );
		}
		glEnd();
	}
	return nvert;
};

int drawSphere_oct( int n, float r, const Vec3f& pos, bool wire ){
	int nvert=0;
	Vec3f px,mx,py,my,pz,mz;
	px.set( 1,0,0); py.set(0, 1,0); pz.set(0,0, 1);
	mx.set(-1,0,0); my.set(0,-1,0); mz.set(0,0,-1);
	if(wire){
        nvert += drawSphereTriangle_wire( n, r, pos, mz, mx, my );
        nvert += drawSphereTriangle_wire( n, r, pos, mz, my, px );
        nvert += drawSphereTriangle_wire( n, r, pos, mz, px, py );
        nvert += drawSphereTriangle_wire( n, r, pos, mz, py, mx );
        nvert += drawSphereTriangle_wire( n, r, pos, pz, mx, my );
        nvert += drawSphereTriangle_wire( n, r, pos, pz, my, px );
        nvert += drawSphereTriangle_wire( n, r, pos, pz, px, py );
        nvert += drawSphereTriangle_wire( n, r, pos, pz, py, mx );
	}else{
        nvert += drawSphereTriangle( n, r, pos, mz, mx, my );
        nvert += drawSphereTriangle( n, r, pos, mz, my, px );
        nvert += drawSphereTriangle( n, r, pos, mz, px, py );
        nvert += drawSphereTriangle( n, r, pos, mz, py, mx );
        nvert += drawSphereTriangle( n, r, pos, pz, mx, my );
        nvert += drawSphereTriangle( n, r, pos, pz, my, px );
        nvert += drawSphereTriangle( n, r, pos, pz, px, py );
        nvert += drawSphereTriangle( n, r, pos, pz, py, mx );
	}
	return nvert;
};

/*
int drawSphereStrip( Vec3f p, Vec3f ax,  int nPhi. int nThet, float R, float sinTheta1, float sinTheta1 ){
    Vec2f cph, dph, cth,dph;
    dph.fromAngle( );
    dph.fromAngle( );
    cos();
    for(int ith=0; ith<nTheta; ith++){
        glBegin(GL_TRIANGLE_STRIP);
        for(int iph=0; iph<nPhi; iph++){
            glVertex2d();
            glVertex2d();
        }
        glEnd();
    }
};
*/


int  drawCapsula( Vec3f p0, Vec3f p1, float r1, float r2, float theta1, float theta2, float dTheta, int nPhi, bool capped ){
    int nvert=0;
    Vec3f ax   = p1-p0;  float L = ax.normalize();
    Vec3f up,left;       ax.getSomeOrtho(up,left);
    Vec2f cph=Vec2fX, dph;
    dph.fromAngle( 2*M_PI/nPhi );
    // Cylinder
    Vec2f cth,dth;
    float dr = (r2-r1);
    float cv = sqrt(L*L+dr*dr);
    cth.set( L/cv, -dr/cv );
    glBegin(GL_TRIANGLE_STRIP);
    for(int iph=0; iph<(nPhi+1); iph++){
        Vec3f pa = p0 + left*(cph.x*r1) + up*(cph.y*r1);
        Vec3f pb = p1 + left*(cph.x*r2) + up*(cph.y*r2);
        Vec3f na = (left*(cph.x) + up*(cph.y)*cth.x + ax*cth.y)*-1.0;
        Vec3f nb = (left*(cph.x) + up*(cph.y)*cth.x + ax*cth.y)*-1.0;
        glNormal3f(na.x,na.y,na.z); glVertex3f(pa.x,pa.y,pa.z);
        glNormal3f(nb.x,nb.y,nb.z); glVertex3f(pb.x,pb.y,pb.z);
        cph.mul_cmplx(dph);
        nvert+=2;
    }
    glEnd();

    float DTh,h;
    int nTheta;
    // Spherical Cap
    cph=Vec2fX;
    cth.set( L/cv, -dr/cv );
    //dth.fromSin(v1/r1);
    DTh = (-theta1 - asin(cth.y));
    nTheta = (int)(fabs(DTh)/dTheta);
    dth.fromAngle( DTh/nTheta );
    //printf( " cth (%f,%f)  dth (%f,%f) \n", cth.x, cth.y,  dth.x, dth.y );
    r1/=cth.x;
    h  =-cth.y*r1;
    // Left
    for(int ith=0; ith<(nTheta+1); ith++){
        Vec2f cth_ = Vec2f::mul_cmplx(cth,dth);
        glBegin(GL_TRIANGLE_STRIP);
        //glBegin(GL_LINES);
        for(int iph=0; iph<(nPhi+1); iph++){
            Vec3f pa = p0 + (left*(cph.x*r1) + up*(cph.y*r1))*cth.x  + ax*(h+cth.y*r1);
            Vec3f pb = p0 + (left*(cph.x*r1) + up*(cph.y*r1))*cth_.x + ax*(h+cth_.y*r1);
            Vec3f na = (left*(cph.x) + up*(cph.y)*cth.x  + ax*cth.y)*1.0;
            Vec3f nb = (left*(cph.x) + up*(cph.y)*cth_.x + ax*cth_.y)*1.0;
            glNormal3f(na.x,na.y,na.z); glVertex3f(pa.x,pa.y,pa.z);
            glNormal3f(nb.x,nb.y,nb.z); glVertex3f(pb.x,pb.y,pb.z);
            nvert+=2;
            //na.mul(0.2);
            //glVertex3f(pa.x,pa.y,pa.z);   glVertex3f(pa.x+na.x,pa.y+na.y,pa.z+na.z);
            cph.mul_cmplx(dph);
        }
        glEnd();
        //printf( "%i cth (%f,%f)  cth_ (%f,%f) \n", ith, cth.x, cth.y,  cth_.x, cth_.y );
        cth=cth_;
    }
    //return 0;
    cph=Vec2fX;
    cth.set( L/cv, -dr/cv );
    //cth = Vec2fX;
    //cth.set( dr/cv, L/cv);
    //dth.fromAngle( asin(v2/r2)/nTheta );
    DTh    = (theta2-asin(cth.y));
    nTheta = (int)(fabs(DTh)/dTheta);
    dth.fromAngle(DTh/nTheta );
    r2/= cth.x;
    h  =-cth.y*r2;
    // Right
    for(int ith=0; ith<(nTheta+1); ith++){
        Vec2f cth_ = Vec2f::mul_cmplx(cth,dth);
        glBegin(GL_TRIANGLE_STRIP);
        for(int iph=0; iph<(nPhi+1); iph++){
            Vec3f pa = p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth.x  + ax*(h+cth.y*r2);
            Vec3f pb = p1 + (left*(cph.x*r2) + up*(cph.y*r2))*cth_.x + ax*(h+cth_.y*r2);
            Vec3f na = (left*(cph.x) + up*(cph.y)*cth.x  + ax*cth.y)*-1.0;
            Vec3f nb = (left*(cph.x) + up*(cph.y)*cth_.x + ax*cth_.y)*-1.0;
            glNormal3f(na.x,na.y,na.z); glVertex3f(pa.x,pa.y,pa.z);
            glNormal3f(nb.x,nb.y,nb.z); glVertex3f(pb.x,pb.y,pb.z);
            nvert+=2;
            cph.mul_cmplx(dph);
        }
        glEnd();
        cth=cth_;
    }
    return nvert;
}


int drawParaboloid     ( Vec3f p0, Vec3f ax, float r, float l, float nR, int nPhi, bool capped ){
    int nvert=0;
    float L = ax.normalize();
    Vec3f up,left;       ax.getSomeOrtho(up,left);
    Vec2f cph=Vec2fX, dph;
    dph.fromAngle( 2*M_PI/nPhi );
    float dr = r/nR;
    float a  = 1.0; // TODO
    for(int ir=0; ir<nR; ir++){
        glBegin(GL_TRIANGLE_STRIP);
        //glBegin(GL_LINES);
        float r1 =(ir-1)*dr;
        float r2 =(ir  )*dr;
        for(int iph=0; iph<(nPhi+1); iph++){
            float h1 = a*r1*r1;
            float h2 = a*r2*r2;
            Vec3f pa = p0 + left*(cph.x*r1) + up*(cph.y*r1) + ax*h1;
            Vec3f pb = p0 + left*(cph.x*r2) + up*(cph.y*r2) + ax*h2;
            //Vec3f na = left*(cph.x) + up*(cph.y)*cth.x  + ax*cth.;
            //Vec3f nb = left*(cph.x) + up*(cph.y)*cth_.x + ax*cth_.y;
            //glNormal3f(na.x,na.y,na.z);
            glVertex3f(pa.x,pa.y,pa.z);
            //glNormal3f(nb.x,nb.y,nb.z);
            glVertex3f(pb.x,pb.y,pb.z);
            nvert+=2;
            cph.mul_cmplx(dph);
        }
        glEnd();
    }
    return nvert;
}


int drawCircleAxis( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R, float dca, float dsa ){
    int nvert=0;
    Vec3f v; v.set(v0);
    glBegin( GL_LINE_LOOP );
    for( int i=0; i<n; i++ ){
        glVertex3f( pos.x+v.x*R, pos.y+v.y*R, pos.z+v.z*R ); nvert++;
        //printf( " drawCircleAxis %i (%3.3f,%3.3f,%3.3f) \n", i, v.x, v.y, v.z );
        v.rotate_csa( dca, dsa, uaxis );
    }
    glEnd();
    return nvert;
}

int drawCircleAxis( int n, const Vec3f& pos, const Vec3f& v0, const Vec3f& uaxis, float R ){
    float dphi = 2*M_PI/n;
    float dca  = cos( dphi );
    float dsa  = sin( dphi );
    return drawCircleAxis( n, pos, v0, uaxis, R, dca, dsa );
}

/*
int drawSphereOctLines( int n, float R, const Vec3f& pos ){
	int nvert=0;
    float dphi = 2*M_PI/n;
    float dca  = cos( dphi );
    float dsa  = sin( dphi );
    nvert += drawCircleAxis( n, pos, {0,1,0}, {1.0d,0.0d,0.0d}, R, dca, dsa );
    nvert += drawCircleAxis( n, pos, {0,0,1}, {0.0d,1.0d,0.0d}, R, dca, dsa );
    nvert += drawCircleAxis( n, pos, {1,0,0}, {0.0d,0.0d,1.0d}, R, dca, dsa );
	return nvert;
}
*/

int drawSphereOctLines( int n, float R, const Vec3f& pos, const Mat3f& rot, bool bRGB ){
	int nvert=0;
    float dphi = 2*M_PI/n;
    float dca  = cos( dphi );
    float dsa  = sin( dphi );
    if(bRGB)glColor3f(1,0,0);
    nvert += drawCircleAxis( n, pos, rot.b, rot.a, R, dca, dsa );
    if(bRGB)glColor3f(0,1,0);
    nvert += drawCircleAxis( n, pos, rot.c, rot.b, R, dca, dsa );
    if(bRGB)glColor3f(0,0,1);
    nvert += drawCircleAxis( n, pos, rot.a, rot.c, R, dca, dsa );
	return nvert;
}

void drawPlanarPolygon( int n, const int * inds, const Vec3d * points ){
    if( n < 3 ) return;

    Vec3f a,b,c,normal;
    a = (Vec3f) points[inds[0]];
    b = (Vec3f) points[inds[1]];
    c = (Vec3f) points[inds[2]];
    normal.set_cross( a-b, b-c );
    normal.normalize( );

    glBegin( GL_TRIANGLE_FAN );
    glNormal3f( normal.x, normal.y, normal.z );
    glVertex3f( a.x, a.y, a.z );
    glVertex3f( b.x, b.y, b.z );
    glVertex3f( c.x, c.y, c.z );
    for( int i=3; i<n; i++ ){
        convert( points[inds[i]], a );
        glVertex3f( a.x, a.y, a.z );
        //average.add( a );
    }
    glEnd();
}

void drawPolygonNormal( int n, const int * inds, const Vec3d * points ){
    if( n < 3 ) return;

    Vec3f a,b,c,normal;
    a = (Vec3f) points[inds[0]];
    b = (Vec3f) points[inds[1]];
    c = (Vec3f) points[inds[2]];
    normal.set_cross( a-b, b-c );
    normal.normalize( );

    glBegin( GL_LINES );
        glVertex3f( a.x, a.y, a.z );
        glVertex3f( a.x+normal.x, a.y+normal.y, a.z+normal.z );
    glEnd();
}

void drawPolygonBorder( int n, const int * inds, const Vec3d * points ){
    glBegin( GL_LINE_LOOP );
    Vec3f a;
    for( int i=0; i<n; i++ ){
        convert( points[inds[i]], a );
        glVertex3f( a.x, a.y, a.z );
    }
    glEnd();
}

void drawPlanarPolygon( int ipl, Mesh& mesh ){
    Polygon * pl = mesh.polygons[ipl];
    Draw3D:: drawPlanarPolygon( pl->ipoints.size(), &pl->ipoints.front(), &mesh.points.front() );
}

void drawPolygonNormal( int ipl, Mesh& mesh ){
    Polygon * pl = mesh.polygons[ipl];
    Draw3D:: drawPolygonNormal( pl->ipoints.size(), &pl->ipoints.front(), &mesh.points.front() );
}

void drawPolygonBorder( int ipl, Mesh& mesh ){
    Polygon * pl = mesh.polygons[ipl];
    Draw3D:: drawPolygonBorder( pl->ipoints.size(), &pl->ipoints.front(), &mesh.points.front() );
}


void drawPoints( int n, const  Vec3d * points, float sz ){
    if(sz<=0){
        glBegin( GL_POINTS );
        for( int i=0; i<n; i++ ){
            Vec3f a;
            convert( points[i], a );
            glVertex3f( a.x, a.y, a.z );
        }
        glEnd();
	}else{
        glBegin( GL_LINES );
        for( int i=0; i<n; i++ ){
            Vec3f vec;
            convert( points[i], vec );
            glVertex3f( vec.x-sz, vec.y, vec.z ); glVertex3f( vec.x+sz, vec.y, vec.z );
            glVertex3f( vec.x, vec.y-sz, vec.z ); glVertex3f( vec.x, vec.y+sz, vec.z );
            glVertex3f( vec.x, vec.y, vec.z-sz ); glVertex3f( vec.x, vec.y, vec.z+sz );
        }
        glEnd();
	}
}

void drawLines( int nlinks, const  int * links, const  Vec3d * points ){
	int n2 = nlinks<<1;
	glBegin( GL_LINES );
	for( int i=0; i<n2; i+=2 ){
		//drawLine( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
		Vec3f a,b;
		convert( points[links[i  ]], a );
        convert( points[links[i+1]], b );
        glVertex3f( a.x, a.y, a.z );
        glVertex3f( b.x, b.y, b.z );
	}
	glEnd();
}

void drawMeshWireframe(const CMesh& msh){ drawLines( msh.nedge, (int*)msh.edges, msh.verts ); }

    void drawTriangles( int nlinks, const int * links, const Vec3d * points, int mode ){
        int n2 = nlinks*3;
        if((mode==2)||(mode==1)){ glBegin( GL_LINES ); }else{ glBegin( GL_TRIANGLES ); };
        for( int i=0; i<n2; i+=3 ){
            //drawTriangle( points[links[i]], points[links[i+1]], points[links[i+2]] );
            //printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
            Vec3f a,b,c,nor;
            convert( points[links[i  ]], a );
            convert( points[links[i+1]], b );
            convert( points[links[i+2]], c );
            //printf( " %i (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) (%3.3f,%3.3f,%3.3f) \n", i, a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z  );
            nor.set_cross( a-b, b-c );
            nor.normalize( );
            if(mode==2){
                Vec3f cog = (a+b+c)*(1./3.);
                vertex( cog );
                vertex( cog + nor );
            }else if(mode==1){
                vertex( a );vertex( b );
                vertex( b );vertex( c );
                vertex( c );vertex( a );
            }else{
                //glNormal3f( normal.x, normal.y, normal.z );
                //glVertex3f( a.x, a.y, a.z );
                //glVertex3f( b.x, b.y, b.z );
                //glVertex3f( c.x, c.y, c.z );
                normal(nor);
                vertex( a ); vertex( b ); vertex( c );

            }
        }
        glEnd();
    }


    void drawPolygons( int nlinks, const int * ns, const int * links, const Vec3d * points ){
        const int * inds = links;
        for( int i=0; i<nlinks; i++ ){
            int ni = ns[i];
            drawPlanarPolygon( ni, inds, points );
            inds += ni;
            /*
            glBegin( GL_TRIANGLE_FAN );
            Vec3f a,b,c,normal;
            convert( points[links[i  ]], a );
            convert( points[links[i+1]], b );
            convert( points[links[i+2]], c );
            normal.set_cross( a-b, b-c );
            normal.normalize( );
            glNormal3f( normal.x, normal.y, normal.z );
            glVertex3f( a.x, a.y, a.z );
            glVertex3f( b.x, b.y, b.z );
            glVertex3f( c.x, c.y, c.z );
            int ni = ns[i];
            if( ni > 3 ){
                convert( points[links[i  ]], a );
                glVertex3f( a.x, a.y, a.z );
            }
            glEnd();
            */
        }

    }

    void drawMeshPolygons( const CMesh& msh ){ drawPolygons( msh.nfaces,  msh.ngons, msh.faces, msh.verts ); };
    void drawMesh( const CMesh& msh, uint32_t cpoly, uint32_t cwire ){
        if( cpoly>0 ){ Draw::setRGBA(cpoly); drawPolygons( msh.nfaces,  msh.ngons, msh.faces, msh.verts );  };
        if( cwire>0 ){ Draw::setRGBA(cwire); drawLines   ( msh.nedge, (int*)msh.edges, msh.verts );         };
    }

    void drawKite( const Vec3f& pos, const Mat3f& rot, double sz ){
	    //drawLine( const Vec3d& p1, const Vec3d& p2 );
	    //double sz = sqrt( area );
	    //glEnable (GL_LIGHTING);
	    //glColor3f( 0.5f,0.5f,0.5f );
	    glBegin  (GL_QUADS);
		    glNormal3d( rot.b.x, rot.b.y, rot.b.z );
		    glVertex3d( pos.x-sz*rot.a.x, pos.y-sz*rot.a.y, pos.z-sz*rot.a.z );
		    glVertex3d( pos.x-sz*rot.c.x, pos.y-sz*rot.c.y, pos.z-sz*rot.c.z );
		    glVertex3d( pos.x+sz*rot.a.x, pos.y+sz*rot.a.y, pos.z+sz*rot.a.z );
		    glVertex3d( pos.x+sz*rot.c.x, pos.y+sz*rot.c.y, pos.z+sz*rot.c.z );
	    glEnd();
    }

    void drawPanel( const Vec3f& pos, const Mat3f& rot, const Vec2f& sz ){
	    //drawLine( const Vec3d& p1, const Vec3d& p2 );
	    //double sz = sqrt( area );
	    //glEnable (GL_LIGHTING);
	    //glColor3f( 0.5f,0.5f,0.5f );
	    glBegin  (GL_QUADS);
		    glNormal3f( rot.b.x, rot.b.y, rot.b.z );
		    Vec3f p;
		    glTexCoord2f(0.0,1.0); p=pos-rot.a*sz.a + rot.c*sz.b; glVertex3f( p.x, p.y, p.z );
		    glTexCoord2f(0.0,0.0); p=pos-rot.a*sz.a - rot.c*sz.b; glVertex3f( p.x, p.y, p.z );
		    glTexCoord2f(1.0,0.0); p=pos+rot.a*sz.a - rot.c*sz.b; glVertex3f( p.x, p.y, p.z );
		    glTexCoord2f(1.0,1.0); p=pos+rot.a*sz.a + rot.c*sz.b; glVertex3f( p.x, p.y, p.z );
	    glEnd();
    }

    void drawVectorArray(int n,const  Vec3d* ps,const  Vec3d* vs, double sc, double lmax ){
        glBegin(GL_LINES);
        double l2max=sq(lmax/sc);
        for(int i=0; i<n; i++){
            if(lmax>0){ if(vs[i].norm2()>l2max ) continue; }
            Vec3d p=ps[i];        glVertex3f(p.x,p.y,p.z);
            p.add_mul( vs[i], sc); glVertex3f(p.x,p.y,p.z);
        }
        glEnd();
    }


    void drawScalarArray(int n,const Vec3d* ps,const double* vs, double vmin, double vmax, const uint32_t * colors, int ncol ){
        glBegin(GL_POINTS);
        double sc = 1/(vmax-vmin);
        for(int i=0; i<n; i++){
            Vec3d p=ps[i];
            double c = (vs[i]-vmin)*sc;
            if(colors){
                Draw::colorScale(c,ncol,colors);
            }else{
                glColor3f(c,c,c);
            }
            glVertex3f(p.x,p.y,p.z);
            //printf( "i %i p(%g,%g,%g) v: %g c: %g\n", i, p.x,p.y,p.z, vs[i], c );
        }
        glEnd();
    }

    void drawScalarField( Vec2i ns, const Vec3d* ps,const  double* data,  double vmin, double vmax, const uint32_t * colors, int ncol ){
        //printf( " debug_draw_GridFF \n" );
        double z0  = 1.5;
        double dz0 = 0.1;
        double clsc = 1/(vmax-vmin);
        glShadeModel(GL_SMOOTH);
        //glEnable( GL_POLYGON_SMOOTH);
        for(int iy=1;iy<ns.y;iy++){
            glBegin( GL_TRIANGLE_STRIP );
            for(int ix=0;ix<ns.x;ix++){
                Vec3d p;
                int i = (iy-1)*ns.x + ix;
                //glColor3f ( data[i].x+0.5, data[i].y+0.5, 0.5 );
                double c = clamp( clsc*(data[i]-vmin), 0, 1 );
                if(colors){ Draw::colorScale( c,ncol,colors); }else{ glColor3f(c,c,c); }
                //p = (gsh.dCell.a*(ix + (gsh.n.x*-0.5))) + (gsh.dCell.b*(iy-1 + (gsh.n.y*-0.5) ));
                p = ps[i];
                glVertex3f(p.x,p.y,p.z);

                i += ns.x;
                //glColor3f ( data[i].x+0.5, data[i].y+0.5, 0.5 );
                c = clamp(  clsc*(data[i]-vmin), 0, 1 );
                if(colors){ Draw::colorScale( c,ncol,colors); }else{ glColor3f(c,c,c); }
                p = ps[i];
                glVertex3f(p.x,p.y,p.z);
            }
            glEnd();
        }
    }

    void drawScalarGrid(Vec2i ns, const Vec3d& p0, const Vec3d& a, const Vec3d& b,const double* data,  double vmin, double vmax, const uint32_t * colors, int ncol ){
        //printf( " debug_draw_GridFF \n" );
        double z0  = 1.5;
        double dz0 = 0.1;
        double clsc = 1/(vmax-vmin);
        glShadeModel(GL_SMOOTH);
        //glEnable( GL_POLYGON_SMOOTH);
        for(int iy=1;iy<ns.y;iy++){
            glBegin( GL_TRIANGLE_STRIP );
            for(int ix=0;ix<ns.x;ix++){
                Vec3d p;
                int i = (iy-1)*ns.x + ix;
                //glColor3f ( data[i].x+0.5, data[i].y+0.5, 0.5 );
                double c = clamp( clsc*(data[i]-vmin), 0, 1 );
                if(colors){ Draw::colorScale( c,ncol,colors); }else{ glColor3f(c,c,c); }
                //p = (gsh.dCell.a*(ix + (gsh.n.x*-0.5))) + (gsh.dCell.b*(iy-1 + (gsh.n.y*-0.5) ));
                p = a*ix + b*(iy-1) + p0;
                glVertex3f(p.x,p.y,p.z);

                i += ns.x;
                //glColor3f ( data[i].x+0.5, data[i].y+0.5, 0.5 );
                c = clamp(  clsc*(data[i]-vmin), 0, 1 );
                if(colors){ Draw::colorScale( c,ncol,colors); }else{ glColor3f(c,c,c); }
                p.add(b);
                glVertex3f(p.x,p.y,p.z);
            }
            glEnd();
        }
    }

    void drawColorScale( int n, const Vec3d& p0, const Vec3d& fw, const Vec3d& up, const uint32_t * colors, int ncol ){
        //printf( " debug_draw_GridFF \n" );
        double step = 1./(n-1);
        //glBegin( GL_LINE_STRIP );
        glShadeModel(GL_SMOOTH);
        glBegin( GL_TRIANGLE_STRIP );
        for(int iy=0;iy<n;iy++){
            double c = iy*step;
            Draw::colorScale( c,ncol,colors);
            //glColor3f(1.,1.,1.);
            Vec3d p = fw*c + p0;
            //printf( "%i %g (%g,%g,%g)\n", iy, c, p.x,p.y,p.z );
            glVertex3f(p.x,p.y,p.z);
            p.add(up);
            //printf( "%i %g (%g,%g,%g)\n", iy, c, p.x,p.y,p.z );
            glVertex3f(p.x,p.y,p.z);
        }
        glEnd();
    }



    inline void simplex_deriv(
        const Vec2d& da, const Vec2d& db,
        double p7, double p8, double p9, double p4, double p5, double p6, double p2, double p3,
        Vec2d& deriv
    ){
        deriv.x = da.x*(p6-p4) + db.x*(p8-p2) + (da.x-db.x)*(p3-p7);
        deriv.y = da.y*(p6-p4) + db.y*(p8-p2) + (da.y-db.y)*(p3-p7);
    }

    void drawSimplexGrid( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs, const double * clrs, int ncolors, const uint32_t * cscale ){
        //const double * heights
        //const double * colors
        //const Vec2d  * normals

        Vec2d pa; pa.set(0.0d);
        if( !cscale ){ cscale=&Draw::colors_rainbow[0]; ncolors=Draw::ncolors; }
        int ii=0;
        glNormal3f(0.0f,1.0f,0.0f);
        for (int ia=0; ia<(na-1); ia++){
            glBegin( GL_TRIANGLE_STRIP );
            Vec2d p; p.set(pa);
            for (int ib=0; ib<nb; ib++){
                double h=0.0d;
                //printf( " %i %i %i (%3.3f,%3.3f) %f %f \n", ia, ib, ii, p.x, p.y, hs[ii], clrs[ii] );
                if(clrs) Draw::colorScale( clrs[ii], ncolors, cscale );
                //if(hs){ simplex_deriv(); glNormal3f(0.0f,1.0f,0.0f); }
                if(hs){ h=hs[ii]; }
                glVertex3f( (float)(p.x), (float)(p.y), (float)h );

                if(clrs) Draw::colorScale( clrs[ii+nb], ncolors, cscale );
                //if(hs){ simplex_deriv(); glNormal3f(0.0f,1.0f,0.0f); }
                if(hs){ h=hs[ii+nb]; }
                glVertex3f( (float)(p.x+da.x), (float)(p.y+da.y), (float)h );
                p.add(db);
                ii++;
            }
            pa.add(da);
        glEnd();
        }

    }

    void drawSimplexGridLines( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs ){
        Vec2d p,pa; pa.set(0.0d);
        for (int ia=0; ia<(na-1); ia++){
            glBegin( GL_LINE_STRIP );
            p.set(pa);
            for (int ib=0; ib<nb; ib++){
                glVertex3f( (float)(p.x),      (float)(p.y),      (float)hs[ia*nb+ib] );
                p.add(db);
            }
            glEnd();
            p.set(pa);
            glBegin( GL_LINE_STRIP );
            for (int ib=0; ib<nb; ib++){
                int ii=ia*nb+ib;
                glVertex3f( (float)(p.x),      (float)(p.y),      (float)hs[ii   ] );
                glVertex3f( (float)(p.x+da.x), (float)(p.y+da.y), (float)hs[ii+nb] );
                p.add(db);
                ii++;
            }
            glEnd();
            pa.add(da);
        }
        p.set(pa);
        glBegin( GL_LINE_STRIP );
        for (int ib=0; ib<nb; ib++){
            glVertex3f( (float)(p.x),  (float)(p.y), (float)hs[(na-1)*nb+ib] );
            p.add(db);
        }
        glEnd();
    }

    void drawSimplexGridLinesToned( int na, int nb, const Vec2d& da, const Vec2d& db,  const double * hs ){
        Vec2d p,pa; pa.set(0.0d);
        float h;
        for (int ia=0; ia<(na-1); ia++){
            glBegin( GL_LINE_STRIP );
            p.set(pa);
            for (int ib=0; ib<nb; ib++){
                h = (float)hs[ia*nb+ib];
                glColor3f( h,h*4,h*16 ); glVertex3f( (float)(p.x),      (float)(p.y),      h );
                p.add(db);
            }
            glEnd();
            p.set(pa);
            glBegin( GL_LINE_STRIP );
            for (int ib=0; ib<nb; ib++){
                int ii=ia*nb+ib;
                h=(float)hs[ii   ]; glColor3f( h,h*4,h*16 ); glVertex3f( (float)(p.x),      (float)(p.y),      h );
                h=(float)hs[ii+nb]; glColor3f( h,h*4,h*16 ); glVertex3f( (float)(p.x+da.x), (float)(p.y+da.y), h );
                p.add(db);
                ii++;
            }
            glEnd();
            pa.add(da);
        }
        p.set(pa);
        glBegin( GL_LINE_STRIP );
        for (int ib=0; ib<nb; ib++){
            h=(float)hs[(na-1)*nb+ib]; glColor3f( h,h*4,h*16 ); glVertex3f( (float)(p.x),  (float)(p.y), h );
            p.add(db);
        }
        glEnd();
    }


    void drawRectGridLines( Vec2i n, const Vec3d& p0, const Vec3d& da, const Vec3d& db ){
        glBegin( GL_LINES );
        Vec3d p  = p0;
        Vec3d dn = db*n.b;
        for (int ia=0; ia<n.a; ia++){
            glVertex3f( (float)(p .x), (float)(p .y), (float)(p .z) );  Vec3d p_ = p+dn;
            glVertex3f( (float)(p_.x), (float)(p_.y), (float)(p_.z) );
            //printf( "ia (%g,%g,%g) (%g,%g,%g)\n", p.x,p.y,p.z,   p_.x,p_.y,p_.z );
            p.add(da);
        }
        p   = p0;
        dn  = da*n.a;
        for (int ib=0; ib<n.b; ib++){
            glVertex3f( (float)(p .x), (float)(p .y), (float)(p .z) );  Vec3d p_ = p+dn;
            glVertex3f( (float)(p_.x), (float)(p_.y), (float)(p_.z) );
            //printf( "ib (%g,%g,%g) (%g,%g,%g)\n", p.x,p.y,p.z,   p_.x,p_.y,p_.z );
            p.add(db);
        }
        glEnd();
    }

    int drawMesh( const Mesh& mesh  ){
        int nvert=0;
        for( Polygon* pl : mesh.polygons ){
            Draw3D::drawPlanarPolygon( pl->ipoints.size(), &pl->ipoints.front(), &mesh.points.front() );
            nvert +=pl->ipoints.size();
        }
        return nvert;
    }

    void drawText( const char * str, const Vec3f& pos, int fontTex, float textSize, int iend ){
        glDisable    ( GL_LIGHTING   );
        glDisable    ( GL_DEPTH_TEST );
        glShadeModel ( GL_FLAT       );
        glPushMatrix();
            //glMatrixMode(GL_MODELVIEW);
            //glMatrixMode(GL_PROJECTION);
            //printf("-- txt p (%f,%f,%f)\n", pos.x, pos.y, pos.z);
            glTranslatef( pos.x, pos.y, pos.z );
            //Draw::billboardCam( );
            Draw::billboardCamProj( );
            //Draw2D::drawString( inputText.c_str(), 0, 0, textSize, fontTex );
            Draw::drawText( str, fontTex, textSize, iend );
        glPopMatrix();
	}

    void drawText3D( const char * str, const Vec3f& pos, const Vec3f& fw, const Vec3f& up, int fontTex, float textSize, int iend ){
        glDisable    ( GL_LIGHTING   );
        glDisable    ( GL_DEPTH_TEST );
        glShadeModel ( GL_FLAT       );
        glPushMatrix();
            glTranslatef( pos.x, pos.y, pos.z );
            Mat3f rot;
            //rot.fromDirUp(fw,up);
            rot.fromSideUp(fw,up);
            float glmat[16];
            //toGLMatCam( {0.0,0.0,0.0},rot, glmat );
            toGLMat( {0.0,0.0,0.0},rot, glmat );
            glMatrixMode(GL_MODELVIEW);
            //glLoadMatrixf( glmat );
            glMultMatrixf( glmat );
            Draw::drawText( str, fontTex, textSize, iend );
        glPopMatrix();
	}

	void drawInt( const Vec3d& pos, int i, int fontTex, float sz ){
        char str[16];
        sprintf(str,"%i\0",i);
        Draw3D::drawText(str, pos, fontTex, sz, 0);
    }
    void drawDouble( const Vec3d& pos, double f, int fontTex, float sz ){
        char str[24];
        sprintf(str,"%g\0",f);
        Draw3D::drawText(str, pos, fontTex, sz, 0);
    }

	void drawCurve( float tmin, float tmax, int n, Func1d3 func ){
        glBegin(GL_LINE_STRIP);
        float dt = (tmax-tmin)/n;
        for( float t=tmin; t<=tmax; t+=dt ){
            double x,y,z;
            func( t, x, y, z );
            glVertex3f( (float)x, (float)y, (float)z );
        }
        glEnd();
    }

    void drawColorScale( int n, Vec3d pos, Vec3d dir, Vec3d up, void (_colorFunc_)(float f) ){
        glBegin(GL_TRIANGLE_STRIP);
        double d = 1.0/(n-1);
        for(int i=0; i<n; i++){
            double f = i*d;
            _colorFunc_( f );
            //glColor3f(1.0,1.0,1.0);
            Vec3d p = pos + dir*f;
            glVertex3f( (float)(p.x     ),(float)( p.y     ),(float)( p.z     ) );
            glVertex3f( (float)(p.x+up.x),(float)( p.y+up.y),(float)( p.z+up.z) );
            //printf( "(%g,%g,%g) (%g,%g,%g) \n", p.x, p.y, p.z, (float)(pos.x+up.x),(float)( pos.y+up.y),(float)( pos.z+up.z)  );
        }
        glEnd();
    }

// =================
// from drawUtils.h
// =================

void drawBox( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b ){
	glBegin(GL_QUADS);
		glColor3f( r, g, b );
		glNormal3f(0,0,-1); glVertex3f( x0, y0, z0 ); glVertex3f( x1, y0, z0 ); glVertex3f( x1, y1, z0 ); glVertex3f( x0, y1, z0 );
		glNormal3f(0,-1,0); glVertex3f( x0, y0, z0 ); glVertex3f( x1, y0, z0 ); glVertex3f( x1, y0, z1 ); glVertex3f( x0, y0, z1 );
		glNormal3f(-1,0,0); glVertex3f( x0, y0, z0 ); glVertex3f( x0, y1, z0 ); glVertex3f( x0, y1, z1 ); glVertex3f( x0, y0, z1 );
		glNormal3f(0,0,+1); glVertex3f( x1, y1, z1 ); glVertex3f( x0, y1, z1 ); glVertex3f( x0, y0, z1 ); glVertex3f( x1, y0, z1 );
		glNormal3f(0,+1,1); glVertex3f( x1, y1, z1 ); glVertex3f( x0, y1, z1 ); glVertex3f( x0, y1, z0 ); glVertex3f( x1, y1, z0 );
		glNormal3f(+1,0,0); glVertex3f( x1, y1, z1 ); glVertex3f( x1, y0, z1 ); glVertex3f( x1, y0, z0 ); glVertex3f( x1, y1, z0 );
	glEnd();
}

void drawBBox( const Vec3f& p0, const Vec3f& p1 ){
	glBegin(GL_LINES);
		glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p0.y, p0.z );
		glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p0.x, p1.y, p0.z );
		glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p0.x, p0.y, p1.z );
        glVertex3f( p1.x, p1.y, p1.z ); glVertex3f( p0.x, p1.y, p1.z );
		glVertex3f( p1.x, p1.y, p1.z ); glVertex3f( p1.x, p0.y, p1.z );
		glVertex3f( p1.x, p1.y, p1.z ); glVertex3f( p1.x, p1.y, p0.z );
		glVertex3f( p1.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p0.z );
		glVertex3f( p1.x, p0.y, p0.z ); glVertex3f( p1.x, p0.y, p1.z );
		glVertex3f( p0.x, p1.y, p0.z ); glVertex3f( p1.x, p1.y, p0.z );
		glVertex3f( p0.x, p1.y, p0.z ); glVertex3f( p0.x, p1.y, p1.z );
		glVertex3f( p0.x, p0.y, p1.z ); glVertex3f( p1.x, p0.y, p1.z );
		glVertex3f( p0.x, p0.y, p1.z ); glVertex3f(p0.x, p1.y, p1.z );
	glEnd();
}

void drawBBox( const Vec3f& p, float r ){ drawBBox( (Vec3f){p.x-r,p.y-r,p.z-r}, (Vec3f){p.x+r,p.y+r,p.z+r} ); };

void drawTriclinicBox( const Mat3f& lvec, const Vec3f& c0, const Vec3f& c1 ){
    Vec3f p0,p1;
	glBegin(GL_LINES);
        lvec.dot_to({c0.x,c0.y,c0.z},p0);
        lvec.dot_to({c0.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to({c0.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to({c1.x,c0.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
		lvec.dot_to({c1.x,c1.y,c1.z},p0);
        lvec.dot_to({c0.x,c1.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to({c1.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to({c1.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to({c1.x,c0.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to({c1.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to({c0.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to({c0.x,c1.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to({c0.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to({c1.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
	glEnd();
}

void drawTriclinicBoxT( const Mat3f& lvec, const Vec3f& c0, const Vec3f& c1 ){
    Vec3f p0,p1;
	glBegin(GL_LINES);
        lvec.dot_to_T({c0.x,c0.y,c0.z},p0);
        lvec.dot_to_T({c0.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to_T({c0.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to_T({c1.x,c0.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
		lvec.dot_to_T({c1.x,c1.y,c1.z},p0);
        lvec.dot_to_T({c0.x,c1.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to_T({c1.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        lvec.dot_to_T({c1.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to_T({c1.x,c0.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to_T({c1.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to_T({c0.x,c0.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to_T({c0.x,c1.y,c1.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to_T({c0.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
        p0=p1; lvec.dot_to_T({c1.x,c1.y,c0.z},p1); glVertex3f( p0.x, p0.y, p0.z ); glVertex3f( p1.x, p1.y, p1.z );
	glEnd();
}


int makeBoxList( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b  ){
	int ilist=glGenLists(1);
	glNewList( ilist, GL_COMPILE );
		drawBox( x0, x1, y0, y1, z0, z1, r, g, b );
	glEndList();
	return( ilist );
	// don't forget use glDeleteLists( ilist ,1); later
}

void drawAxis( float sc ){
	//glDisable (GL_LIGHTING);
	glBegin   (GL_LINES);
		glColor3f( 1, 0, 0 ); glVertex3f( 0, 0, 0 ); glVertex3f( 1*sc, 0, 0 );
		glColor3f( 0, 1, 0 ); glVertex3f( 0, 0, 0 ); glVertex3f( 0, 1*sc, 0 );
		glColor3f( 0, 0, 1 ); glVertex3f( 0, 0, 0 ); glVertex3f( 0, 0, 1*sc );
	glEnd();
}


}; // namespace Draw3D











