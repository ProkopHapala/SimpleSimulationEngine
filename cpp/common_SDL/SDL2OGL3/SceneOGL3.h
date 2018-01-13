
#ifndef  SceneOGL3_h
#define  SceneOGL3_h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"

#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
//#include <SDL2/SDL.h>

#include "Shader.h"
#include "GLObject.h"

class Camera{ public:
    Vec3f  pos   =(Vec3f){0.0f,0.0f,-50.0f};
    Mat3f  rot   =(Mat3f){1.0f,0.0f,0.0f, 0.0f,1.0f,0.0f,  0.0f,0.0f,1.0f };
    float  zoom  =10.0f;
    float  aspect=1.0;
    float  zmin  =1.0;
    float  zmax  =1000.0;

    bool pointInFrustrum( Vec3f p ){
        p.sub(pos);
        Vec3f c;
        rot.dot_to( p, c );
        float y = c.z/(zoom);
        float x = c.z/(zoom*aspect);
        return  (c.x>-x)&&(c.x<x) && (c.y>-y)&&(c.x<y) && (c.z>zmin)&&(c.z<zmax);
    }

};

inline void setCamera(Shader& sh, const Camera& cam){
    Mat4f camMat,mRot,mPersp;
    mPersp.setPerspective( cam.aspect*cam.zoom, cam.zoom, cam.zmin, cam.zmax );
    mRot.setOne(); mRot.setRot(cam.rot);
    //mPersp.setPerspective( fov, fov*ASPECT_RATIO, 1.0, 1000.0 );
    camMat.set_mmul_TN( mRot, mPersp );
    //camMat.set_mmul( mRot, mPersp );
    //Mat4f camMat; camMat.setPerspective( cam.tg, cam.tg*cam.aspect, cam.zmin, cam.zmax );
    //Mat4f camMat; camMat.setPerspective( 20.0, 20.0, 2.0, 1000.0 );
    sh.set_camPos  ( (GLfloat*)&cam.pos );
    sh.set_camMat  ( (GLfloat*)&camMat );
}

inline void setCameraOrtho(Shader& sh, const Camera& cam){
    //Mat4f camMat,mRot,mPersp;
    Mat4f camMat; camMat.setOrthographic( cam.zoom, cam.zoom*cam.aspect, cam.zmin, cam.zmax );
    sh.set_camPos( (GLfloat*)&cam.pos );
    sh.set_camMat( (GLfloat*)&camMat  );
    //sh.set_modelMat();
}

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

class SceneOGL3{ public:
    GLuint vao;

    //functions

	SceneOGL3(){
        glGenVertexArrays(1, &vao);  				// Allocate and assign a Vertex Array Object to our handle
        glBindVertexArray(vao); 					// Bind our Vertex Array Object as the current used object
	}
	~SceneOGL3(){
        glDeleteVertexArrays(1, &vao);
	}

    virtual void draw( Camera& cam ){};
};

#endif
