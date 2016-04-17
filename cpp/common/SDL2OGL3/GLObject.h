
#ifndef  GLObject_h
#define  GLObject_h

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>

#include "Shader.h"

class GLObject{
	public:
    Shader * shader;
	GLuint vbo[2];  // vertex buffer object
	GLenum draw_mode = GL_TRIANGLE_STRIP;

	int nVert   = 4;
	int vertDim = 2;
	int clrDim  = 3;

	GLfloat * vertexes;
	GLfloat * colors;

	virtual void init   ();
	virtual void draw   ();
	virtual void destroy();

	void draw_default();

};

#endif
