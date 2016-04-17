
#ifndef  SceneNode_h
#define  SceneNode_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <SDL2/SDL.h>

#include "Shader.h"
#include "GLObject.h"

class SceneNode3D{
	public:
	Vec3d  pos;
	Mat3d  rot;
	double scale;
	double radius;
	bool visible;

    std::vector<GLObject* >  objects;
	std::vector<SceneNode3D*> subnodes;

    void render( const Vec3d& camPos_, const Mat3d& camRot_, const Vec2d& tgFrustrum, double scale );
};

#endif
