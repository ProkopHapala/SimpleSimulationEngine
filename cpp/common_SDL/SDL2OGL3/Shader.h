
#ifndef  Shader_h
#define  Shader_h

#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
//#include "IO_utils.h"

#include "Vec3.h"
#include "Mat3.h"
#include "Mat4.h"

class Shader{
	public:
	GLchar *vertexsource   = NULL;
	GLchar *fragmentsource = NULL;
	GLuint vertexshader   =0;
	GLuint fragmentshader =0;
	GLuint shaderprogram  =0;

	GLuint uloc_camMat=0, uloc_modelMat=0, uloc_modelPos=0, uloc_camPos=0, uloc_baseColor=0;

	// ==== function declarations

	void compileShader       ( GLenum shaderType, char* sourceCode, GLuint& shader, char*& errLog        );
	void compileShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint& shader, char*& errLog );
	void init                ( char const * vertName, char const * fragName                              );
	void destory             (                                                                           );

	void getDefaultUniformLocation();

	void setUniformi    ( char * name, int    i);
	void setUniformf    ( char * name, float  f);
	void setUniformVec2f( char * name, Vec2f  v);
	void setUniformVec4f( char * name, Quat4f v);
	void setUniformVec3f( char * name, Vec3f  v);//{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform3fv      (u, 1, modelPos           ); };
	void setUniformMat3f( char * name, Mat3f  m);//{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix3fv(u, 1, GL_FALSE, modelMat ); };
	void setUniformMat4f( char * name, Mat4f  m);//{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix4fv(u, 1, GL_FALSE, camMat   ); };

	inline void set_baseColor( GLfloat * baseColor ){ glUniform4fv      (uloc_baseColor, 1, baseColor );          };
	inline void set_modelPos ( GLfloat * modelPos  ){ glUniform3fv      (uloc_modelPos,  1, modelPos );           };
    inline void set_modelMat ( GLfloat * modelMat  ){ glUniformMatrix3fv(uloc_modelMat,  1, GL_FALSE, modelMat ); };
    inline void set_camMat   ( GLfloat * camMat    ){ glUniformMatrix4fv(uloc_camMat  ,  1, GL_FALSE, camMat   ); };
    inline void set_camPos   ( GLfloat * camPos    ){ glUniform3fv      (uloc_camPos  ,  1, camPos   );           };
    inline void setPose      ( GLfloat * modelPos, GLfloat * modelMat, GLfloat * camMat  ){ set_modelPos(modelPos); set_modelMat(modelMat); set_camMat(camMat); };

};

#endif
