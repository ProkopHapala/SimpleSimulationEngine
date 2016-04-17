
#ifndef  Shader_h
#define  Shader_h

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
//#include "IO_utils.h"

class Shader{
	public:
	GLchar *vertexsource   = NULL;
	GLchar *fragmentsource = NULL;
	GLuint vertexshader;
	GLuint fragmentshader;
	GLuint shaderprogram;

	// ==== function declarations

	void compileShader       ( GLenum shaderType, char* sourceCode, GLuint& shader, char*& errLog        );
	void compileShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint& shader, char*& errLog );
	void init                ( char const * vertName, char const * fragName                              );
	void destory             (                                                                           );

};

#endif
