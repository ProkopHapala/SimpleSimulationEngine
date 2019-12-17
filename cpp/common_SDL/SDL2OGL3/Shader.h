
#ifndef  Shader_h
#define  Shader_h

#include <GL/glew.h>
//#define GL_GLEXT_PROTOTYPES
//#include <GL/gl.h>
//#include "IO_utils.h"

#include "Vec3.h"
#include "Mat3.h"
#include "Mat4.h"

#define GLSL(version, shader)  "#version " #version "\n" #shader
// https://stackoverflow.com/questions/13872544/gcc-stringification-and-inline-glsl
// Multiline commet R"( .... )"   // https://stackoverflow.com/questions/20508534/c-multiline-string-raw-literal

const char DEFAULT_vertex_shader_code[]= GLSL(330,
layout(location = 0) in vec3 vpos;
uniform vec3 modelPos;
uniform mat3 modelMat;
uniform vec3 camPos;
uniform mat4 camMat;
void main(){
	vec3 position_world = modelPos + modelMat * vpos;
	gl_Position         = camMat   * vec4( position_world-camPos, 1 );
};
);

const char DEFAULT_fragment_shader_code[] = GLSL(330,
smooth in vec4 gl_FragCoord;
out     vec3 color;
uniform vec4 baseColor;
void main(){
	color = baseColor.xyz;
	gl_FragDepth = gl_FragCoord.z;
};
);

class Shader{
	public:
	//GLchar *vertexsource   = NULL;
	//GLchar *fragmentsource = NULL;
	GLuint vertexshader   =0;
	GLuint fragmentshader =0;
	GLuint geomShader     =0;
	GLuint shaderprogram  =0;

	GLuint uloc_camMat=0, uloc_modelMat=0, uloc_modelPos=0, uloc_camPos=0, uloc_baseColor=0;

	// ==== function declarations

	void compileShader       ( GLenum shaderType, const char* sourceCode,                     GLuint& shader, char*& errLog  );
	void compileShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint geomShader, GLuint& shader, char*& errLog );
	int  init_str            ( const char * vertexsource , const char * fragmentsource, const char * GSsrc       );
	int  init                ( const char * vertName,      const char * fragName,       const char * GSName=NULL );
	void destory             (                                                                           );

	inline int    init_default (){ return init_str(DEFAULT_vertex_shader_code, DEFAULT_fragment_shader_code, NULL);     };
	inline void   use          (){ glUseProgram(shaderprogram); }
	inline GLuint getUloc( char * name ){ return glGetUniformLocation(shaderprogram,name); };

	void getDefaultUniformLocation();

	void setUniformi    ( const char * name, const int    i)const;
	void setUniformf    ( const char * name, const float  f)const;
	void setUniformVec2f( const char * name, const Vec2f  v)const;
	void setUniformVec4f( const char * name, const Quat4f v)const;
	void setUniformVec3f( const char * name, const Vec3f  v)const;//{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform3fv      (u, 1, modelPos           ); };
	void setUniformMat3f( const char * name, const Mat3f  m)const;//{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix3fv(u, 1, GL_FALSE, modelMat ); };
	void setUniformMat4f( const char * name, const Mat4f  m)const;//{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix4fv(u, 1, GL_FALSE, camMat   ); };

	inline void set_baseColor  ( const GLfloat * baseColor )const{ glUniform4fv      (uloc_baseColor, 1, baseColor );          };
	inline void set_modelPos   ( const GLfloat * modelPos  )const{ glUniform3fv      (uloc_modelPos,  1, modelPos );           };
    inline void set_modelMat   ( const GLfloat * modelMat  )const{ glUniformMatrix3fv(uloc_modelMat,  1, GL_FALSE, modelMat ); };
    inline void set_camMat     ( const GLfloat * camMat    )const{ glUniformMatrix4fv(uloc_camMat  ,  1, GL_FALSE, camMat   ); };
    inline void set_camPos     ( const GLfloat * camPos    )const{ glUniform3fv      (uloc_camPos  ,  1, camPos   );           };
    inline void setModelPose   ( const GLfloat * modelPos, GLfloat * modelMat )const{ set_modelPos(modelPos); set_modelMat(modelMat); };
    inline void setModelAndCam ( const GLfloat * modelPos, GLfloat * modelMat, GLfloat * camMat  )const{ set_modelPos(modelPos); set_modelMat(modelMat); set_camMat(camMat); };

    inline void set_modelPos   ( const Vec3d& pos  )const{ Vec3f pos_; convert(pos,pos_); set_modelPos((float*)&pos_); };
    inline void set_modelMat   ( const Mat3d& rot  )const{ Mat3f rot_; convert(rot,rot_); set_modelMat((float*)&rot_); };
    inline void set_modelMatT  ( const Mat3d& rot  )const{ Mat3f rot_; rot_.setT((Mat3f)rot); set_modelMat((float*)&rot_); };
    inline void setModelPose   ( const Vec3d& modelPos, const Mat3d& modelMat )const{ set_modelPos(modelPos); set_modelMat(modelMat); };
    inline void setModelPoseT  ( const Vec3d& modelPos, const Mat3d& modelMat )const{ set_modelPos(modelPos); set_modelMatT(modelMat); };

    Shader(){};
	Shader( const char * vertexsource , const char * fragmentsource, bool bFile ){
        //printf( "bFile %i \n >>%s<< ==\n==\n==\n  >>%s<< ==\n==\n==\n \n", bFile, vertexsource, fragmentsource );
        if(bFile){ init    (vertexsource, fragmentsource, NULL); }
        else     { init_str(vertexsource, fragmentsource, NULL); }
        getDefaultUniformLocation();
    };

};

#endif
