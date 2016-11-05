
#ifndef  GLObject_h
#define  GLObject_h

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>

#include "Shader.h"

// offset and stride
// http://stackoverflow.com/questions/16380005/opengl-3-4-glvertexattribpointer-stride-and-offset-miscalculation

class GLBuff{
    public:
    GLuint  id    = 0;
    GLuint  dim   = 3;
    bool normalized = GL_FALSE;
    GLuint dtype  = GL_FLOAT; // currently this defualt is not modified
    GLuint stride = 0;        // currently this defualt is not modified
    GLuint offset = 0;        // currently this defualt is not modified
    GLuint target = GL_ARRAY_BUFFER;
    GLuint usage  = GL_STATIC_DRAW; // currently this defualt is not modified

    char kind       = 'v';
    GLuint     vbo  = 0;
    GLfloat * cbuff = NULL;

    inline void setup( int id_, int dim_, bool normalized_, GLfloat * cbuff_, char kind_ ){ id = id_; dim = dim_; normalized = normalized_; cbuff = cbuff_; kind=kind_; }

    inline void toGPU( int n ){
        glGenBuffers( 1, &vbo );
        glBindBuffer( target, vbo );
        glBufferData( target, n*dim*sizeof(GLfloat), cbuff, usage );
    }

    inline void activate(){
        glEnableVertexAttribArray( id );
        glBindBuffer( target, vbo );
        glVertexAttribPointer( id, dim, dtype, normalized, stride, (void*)offset );
    }

};

class GLObject{
	public:
	static constexpr int nbuffs = 3;

    Shader * shader;
    int nVert     = 4;
    GLBuff buffs[3];
	//GLuint vbo   [3] = {0,0,0};  // vertex buffer object
	//GLuint vboDim[3] = {3,3,3};
	GLenum draw_mode = GL_TRIANGLE_STRIP;

	virtual void init   ();
	virtual void draw   ();
	virtual void destroy();

	void draw_default();

};

#endif
