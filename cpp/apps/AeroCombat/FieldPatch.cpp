
#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include <SDL2/SDL_opengl.h>

#include "fastmath.h"

#include "FieldPatch.h" // THE HEADER

// ===================
// ==== Quat
// ===================

Quad::Quad(  double value_, const Vec3d& p1_, const Vec3d& p2_, const Vec3d& p3_, const Vec3d& p4_ ){
	// printVec(p1_); printVec(p2_); printVec(p3_); printVec(p4_); printf("ps_\n");
	value = value_;
	p1.set(p1_); p2.set(p2_); p3.set(p3_); p4.set(p4_);
	// printVec(p1); printVec(p2); printVec(p3); printVec(p4); printf(" ps \n");
}

void Quad::draw() const {

	glBegin(GL_QUADS);
		glNormal3f( 0, 1, 0 );
        float red0 =0.65;
		float gree0=0.55;
		float blue0=0.35;
		glColor3f ( red0+z1*0.01, gree0, blue0 );	 glVertex3d( x1, z1, y1 );
		glColor3f ( red0+z2*0.01, gree0, blue0 );	 glVertex3d( x2, z2, y2 );
		glColor3f ( red0+z4*0.01, gree0, blue0 );	 glVertex3d( x4, z4, y4 );
		glColor3f ( red0+z3*0.01, gree0, blue0 );	 glVertex3d( x3, z3, y3 );
	glEnd();

/*
    //printf("QUAD (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f) (%3.3f,%3.3f)\n", x1,y1,  x2,y2, x3,y3, x4,y4);
	glBegin(GL_LINE_LOOP);
		glNormal3f( 0, 1, 0 );
		glVertex3d( x1, z1, y1 );
		glVertex3d( x2, z2, y2 );
        glVertex3d( x3, z3, y3 );
        glVertex3d( x4, z4, y4 );
	glEnd();
*/

};

// ===================
// ==== FieldPatch
// ===================

void FieldPatch::divide_1( int nlevels, int level, const Quad& rc ){

	int dlevel = nlevels - level;
	double fx = 0.5+randf(-dmax,dmax); double mfx = 1-fx;
	double fy = 0.5+randf(-dmax,dmax); double mfy = 1-fy;
	Vec3d top,bot,lft,rgt,ctr;
	top.set( (rc.p1+rc.p2)*0.5 );
	bot.set( (rc.p3+rc.p4)*0.5 );
	lft.set( (rc.p1+rc.p3)*0.5 );
	rgt.set( (rc.p2+rc.p4)*0.5 );
	ctr.set(  (rc.p1*mfx+rc.p2*fx)*mfy + (rc.p3*mfx+rc.p4*fx)*fy );
	double szscale = ( 4.0/( dlevel + 0.5 ) );
	ctr.z += randf( -dmheight, dheight )*szscale;

	divide( nlevels, level-1, Quad( rc.value + randf(-dval,dval),   rc.p1,  top,   lft,    ctr   ) );
	divide( nlevels, level-1, Quad( rc.value + randf(-dval,dval),   top,    rc.p2, ctr,    rgt   ) );
	divide( nlevels, level-1, Quad( rc.value + randf(-dval,dval),   lft,    ctr,   rc.p3,  bot   ) );
	divide( nlevels, level-1, Quad( rc.value + randf(-dval,dval),   ctr,    rgt,   bot,    rc.p4 ) );

};

void FieldPatch::divide( int nlevels, int level, const Quad& rc ){
	//printVec(rc.p1); printVec(rc.p2); printVec(rc.p3); printVec(rc.p4); printf(" rc.ps \n");
	double rnd  = randf();
	double area = rc.area();
	printf("%f %f\n", area, minarea);
	if((level>0)&&( area>minarea )){
        printf("here\n");
		if      ( rnd < thresh[0] ){ divide_1( nlevels, level, rc ); }
		else{
			if( level>(nlevels-startkill) ){ divide_1( nlevels, level, rc );  }
			else                           { rc.draw(); quadCount++;    };
		};
	}else{ rc.draw(); quadCount++; };
};


int FieldPatch::makeList( int nlevels, const Quad& quad ){
	int ilist=glGenLists(1);
	quadCount = 0;
	glNewList( ilist, GL_COMPILE );
		divide( nlevels, nlevels, quad );
	glEndList();
	//exit(0);
	//printf(" FieldPatch made %i quads \n", quadCount );
	return( ilist );
}




