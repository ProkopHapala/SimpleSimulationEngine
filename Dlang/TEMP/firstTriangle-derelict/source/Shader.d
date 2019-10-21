module Shader;

import derelict.opengl3.gl3;

import core.stdc.stdlib;
import std.stdio;
import std.string;
//import std.conv : to;

// ============================
// ====== Free Functions
// ============================

GLuint makeShader( GLenum shaderType, const string sourceCode, ref char* errLog ){
    GLuint shader = glCreateShader( shaderType );
    errLog = null;
    int isCompiled;
    auto vertexZeroTerminated = toStringz(sourceCode);
    glShaderSource( shader, 1, &vertexZeroTerminated, null);
    //glShaderSource ( shader, 1, cast(const(GLchar*)*) &(sourceCode.ptr), null );
    glCompileShader( shader );
    glGetShaderiv( shader, GL_COMPILE_STATUS, &isCompiled );
    if( isCompiled == false)    {
        int maxLength;
        glGetShaderiv( shader, GL_INFO_LOG_LENGTH, &maxLength );
        errLog = cast(char *)malloc(maxLength); 	// The maxLength includes the NULL character
        glGetShaderInfoLog( shader, maxLength, &maxLength, errLog );
    }
    return shader;
}

GLuint makeShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint geomShader, ref char* errLog ){
    GLuint shaderprogram = glCreateProgram();
    glAttachShader(shaderprogram, vertexshader );
    glAttachShader(shaderprogram, fragmentshader );
    if(geomShader) glAttachShader(shaderprogram, geomShader );
    glBindAttribLocation(shaderprogram, 0, "in_Position");
    //glBindAttribLocation(shaderprogram, 1, "in_Color");
    glLinkProgram(shaderprogram);
    int isLinked;
    glGetProgramiv(shaderprogram, GL_LINK_STATUS, cast(int *)&isLinked);
    if( isLinked == false)    {
        int maxLength;
        glGetProgramiv(shaderprogram, GL_INFO_LOG_LENGTH, &maxLength);  			 			// Noticed that glGetProgramiv is used to get the length for a shader program, not glGetShaderiv.
        errLog = cast(char *)malloc(maxLength);   									// The maxLength includes the NULL character
        glGetProgramInfoLog(shaderprogram, maxLength, &maxLength, errLog );		// Notice that glGetProgramInfoLog, not glGetShaderInfoLog.
        //printf( " Error linking shaderProgram : \n" );
        //printf( " %s \n", errLog );
    }
    return shaderprogram;
}

// ============================
// ====== class Shader
// ============================

class Shader{ public:

	//GLchar *vertexsource   = null;
	//GLchar *fragmentsource = null;
	GLuint vertexshader   =0;
	GLuint fragmentshader =0;
	GLuint geomShader     =0;
	GLuint shaderprogram  =0;

	GLuint uloc_camMat=0, uloc_modelMat=0, uloc_modelPos=0, uloc_camPos=0, uloc_baseColor=0;

void checkError( char * errLog, const char * comment ){
    if( errLog != null ){
        writeln( comment, errLog );
        free( errLog );
        exit(-1);
    }
}

this( const string vertexsource , const string fragmentsource, const string GSsrc ){
    char * errLog = null;
    geomShader = 0;
    vertexshader   = makeShader( GL_VERTEX_SHADER,   vertexsource,    errLog ); checkError( errLog, " Error in Vertex Shader: \n" );
    fragmentshader = makeShader( GL_FRAGMENT_SHADER, fragmentsource,  errLog ); checkError( errLog, " Error in Fragment Shader: \n" );
    if( GSsrc  ){ geomShader = makeShader( GL_GEOMETRY_SHADER, GSsrc,  errLog  ); checkError( errLog, " Error in Geometry Shader: \n" ); }
    shaderprogram = makeShaderProgram( vertexshader, fragmentshader, geomShader, errLog ); checkError( errLog, " Error in linking of Shader Program: \n" );
}

/*
int this( const char * vertName, const char * fragName, const char * GSName ){
    if( GSName ){ printf( "loading shaders: Frag: %s Vert: %s Geom: %s \n", vertName, fragName, GSName ); }
    else        { printf( "loading shaders: Frag: %s Vert: %s\n",           vertName, fragName ); }
    int ret=0;
    char * vertexsource   = filetobuf( vertName );  if( vertexsource   == null ){ printf( "cannot load %s \n", vertName ); exit(1); }
    char * fragmentsource = filetobuf( fragName );  if( fragmentsource == null ){ printf( "cannot load %s \n", fragName ); exit(1); }
    char * GSsrc = null;
    if( GSName ){ GSsrc = filetobuf( GSName );  if( GSsrc == null ){ printf( "cannot load %s \n", GSName ); exit(1); } }
    ret = init_str( vertexsource , fragmentsource, GSsrc );
    //delete [] vertexsource;
    //delete [] fragmentsource;
    return ret;
}
*/

~this(){
    // Cleanup all the things we bound and allocated
    glUseProgram(0);
    glDisableVertexAttribArray(0);
    glDisableVertexAttribArray(1);
    glDetachShader(shaderprogram, vertexshader);
    glDetachShader(shaderprogram, fragmentshader);
    glDeleteProgram(shaderprogram);
    glDeleteShader(vertexshader);
    glDeleteShader(fragmentshader);
}

void getDefaultUniformLocation(){
    uloc_camMat    = glGetUniformLocation( shaderprogram, "camMat"    );
    uloc_camPos    = glGetUniformLocation( shaderprogram, "camPos"    );
    uloc_modelMat  = glGetUniformLocation( shaderprogram, "modelMat"  );
    uloc_modelPos  = glGetUniformLocation( shaderprogram, "modelPos"  ); // glUniform3fv      (uloc, 1, modelPos );
    uloc_baseColor = glGetUniformLocation( shaderprogram, "baseColor" );
}

void use(){ glUseProgram(shaderprogram); }

/*
void setUniformi    ( const char * name, const int    i)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform1iv      (u, 1, (GLint*  )&i           ); };
void setUniformf    ( const char * name, const float  f)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform1fv      (u, 1, (GLfloat*)&f           ); };
void setUniformVec2f( const char * name, const Vec2f  v)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform2fv      (u, 1, (GLfloat*)&v           ); };
void setUniformVec3f( const char * name, const Vec3f  v)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform3fv      (u, 1, (GLfloat*)&v           ); };
void setUniformVec4f( const char * name, const Quat4f v)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform4fv      (u, 1, (GLfloat*)&v           ); };
void setUniformMat3f( const char * name, const Mat3f  m)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix3fv(u, 1, GL_FALSE, (GLfloat*)&m ); };
void setUniformMat4f( const char * name, const Mat4f  m)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix4fv(u, 1, GL_FALSE, (GLfloat*)&m ); };}
*/

};
