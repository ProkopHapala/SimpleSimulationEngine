
#ifndef  Shader_h
#define  Shader_h

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
//#include "IO_utils.h"

class Shader{
	public:
	GLchar *vertexsource   = NULL;
	GLchar *fragmentsource = NULL;
	GLuint vertexshader   =0;
	GLuint fragmentshader =0;
	GLuint shaderprogram  =0;

	GLuint uloc_camMat=0, uloc_modelMat=0, uloc_modelPos=0, uloc_baseColor=0;

	// ==== function declarations

	void compileShader       ( GLenum shaderType, char* sourceCode, GLuint& shader, char*& errLog        );
	void compileShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint& shader, char*& errLog );
	void init                ( char const * vertName, char const * fragName                              );
	void destory             (                                                                           );

	void getDefaultUniformLocation();

	inline void set_baseColor( GLfloat * baseColor ){ glUniform4fv      (uloc_baseColor, 1, baseColor );         };
	inline void set_modelPos ( GLfloat * modelPos  ){ glUniform3fv      (uloc_modelPos,  1, modelPos );           };
    inline void set_modelMat ( GLfloat * modelMat  ){ glUniformMatrix3fv(uloc_modelMat,  1, GL_FALSE, modelMat ); };
    inline void set_camMat   ( GLfloat * camMat    ){ glUniformMatrix4fv(uloc_camMat  ,  1, GL_FALSE, camMat   ); };
    inline void setPose      ( GLfloat * modelPos, GLfloat * modelMat, GLfloat * camMat  ){ set_modelPos(modelPos); set_modelMat(modelMat); set_camMat(camMat); };

};

#endif
