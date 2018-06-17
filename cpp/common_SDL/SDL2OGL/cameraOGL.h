
#ifndef  cameraOGL_h
#define  cameraOGL_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Camera.h"
#include "Draw3D.h"

namespace Cam{

void ortho( const Camera& cam, bool zsym ){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    float zmin = cam.zmin; if(zsym) zmin=-cam.zmax;
	glOrtho( -cam.zoom*cam.aspect, cam.zoom*cam.aspect, -cam.zoom, cam.zoom, zmin, cam.zmax );
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, cam.rot, glMat );
	glMultMatrixf( glMat );

	glMatrixMode ( GL_MODELVIEW );
	glLoadIdentity();
	glTranslatef(-cam.pos.x,-cam.pos.y,-cam.pos.z);
}

void perspective( const Camera& cam ){
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    //glFrustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, camDist/zoom, VIEW_DEPTH );
    glFrustum( -cam.aspect*cam.zoom, cam.aspect*cam.zoom, -cam.zoom, cam.zoom, cam.zmin, cam.zmax );
    //glFrustum( -cam.zoom*cam.aspect, cam.zoom*cam.aspect, -cam.zoom, cam.zoom, cam.zmin, cam.zmax );
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, cam.rot, glMat );
	glMultMatrixf( glMat );

	glMatrixMode ( GL_MODELVIEW );
	glLoadIdentity();
	glTranslatef(-cam.pos.x,-cam.pos.y,-cam.pos.z);
    //glTranslatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
    //glTranslatef ( -cam.pos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
}

} // namespace Cam

#endif

