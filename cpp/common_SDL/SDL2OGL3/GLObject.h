
#ifndef  GLObject_h
#define  GLObject_h

#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
//#include <GL/glext.h>

#include "Shader.h"

#include <stdio.h>

// offset and stride
// http://stackoverflow.com/questions/16380005/opengl-3-4-glvertexattribpointer-stride-and-offset-miscalculation

// Vertex Array Object ? - does the same job as GLBuff ?
// https://www.opengl.org/discussion_boards/showthread.php/185088-Vertex-Array-Object-vs-Vertex-Buffer-Object

class GLBuff{ public:
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

    inline void setup( int id_, int dim_, bool normalized_, const void * cbuff_, char kind_ ){ id = id_; dim = dim_; normalized = normalized_; cbuff = (GLfloat*)cbuff_; kind=kind_; }

    inline void setup( int id_, int dim_, bool normalized_, const double * cbuff_, int n, char kind_ ){
        id = id_; dim = dim_; normalized = normalized_; kind=kind_;
        int nd=n*dim;
        cbuff = new GLfloat[nd];
        for(int i=0; i<nd; i++){
            cbuff[i] = (GLfloat)cbuff_[i];
            printf("%i %f \n", i, cbuff[i] );
        }
    }

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

class GLObject{ public:
	static constexpr int nbuffs = 4;

	// ======= variables

    Shader * shader;
    int nVert     = 0;
    int nInd      = 0;
    GLBuff buffs[nbuffs];
    GLuint * index_cbuff = NULL;
    GLuint   index_vbo;
	GLenum draw_mode = GL_TRIANGLES;

    // ======= function def

	virtual void init   ();
	virtual void draw   ();
	virtual void destroy();

	void draw_default();
    void preDraw();
	void afterDraw();
	void draw_instance();

	// ======= functions inline

	void setup( int nVert_ ){
        nVert  = nVert_;
        buffs[0].setup(0,3,GL_FALSE,new GLfloat[nVert*3],'v'); // vertexes
        buffs[1].setup(1,3,GL_FALSE,new GLfloat[nVert*3],'n'); // normals
	}

	inline void setIndexes( int nInd_, int * cbuff_ ){
        nInd = nInd_;
        //int nd=nInd*3;
        index_cbuff = new GLuint[nInd];
        for(int i=0; i<nInd; i++){
            index_cbuff[i] = (GLuint)cbuff_[i];
            printf("%i %i %i \n", i, index_cbuff[i], cbuff_[i] );
        }
    }

    /*
    inline void init( int nVert_ ){
        nVert  = nVert_;
        buffs[0].setup(0,3,GL_FALSE,new GLfloat[nVert*3],'v'); // vertexes
        buffs[1].setup(1,3,GL_FALSE,new GLfloat[nVert*3],'n'); // normals
        init();
    }
    */

    /*
    GLObject(){};
    GLObject( int nVert_ ){
        nVert  = nVert_;
        buffs[0].setup(0,3,GL_FALSE,new GLfloat[nVert*3],'v'); // vertexes
        buffs[1].setup(1,3,GL_FALSE,new GLfloat[nVert*3],'n'); // normals
        init();
    }
    */

};

#endif
