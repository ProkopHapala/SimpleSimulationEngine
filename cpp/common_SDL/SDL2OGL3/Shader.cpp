
#include "IO_utils.h"

#include "Shader.h" // THE HEADER


void checkError( char * errLog, const char * comment ){
    if( errLog != NULL ){
        //printf( " Error in Vertex Shader: \n" );
        printf( "%s %s \n", comment, errLog );
        free( errLog );
        //return -1;
        exit(-1);
    }
}

void Shader::compileShader( GLenum shaderType, const char* sourceCode, GLuint& shader, char*& errLog ){
    shader = glCreateShader( shaderType );
    errLog = NULL;
    int isCompiled;
    glShaderSource( shader, 1, (const GLchar**)&sourceCode, 0);
    glCompileShader( shader );
    glGetShaderiv( shader, GL_COMPILE_STATUS, &isCompiled );
    if( isCompiled == false)    {
        int maxLength;
        glGetShaderiv( shader, GL_INFO_LOG_LENGTH, &maxLength );
        errLog = (char *)malloc(maxLength); 								// The maxLength includes the NULL character
        glGetShaderInfoLog( shader, maxLength, &maxLength, errLog );
        //printf( " Error in compilation of shader %s : \n",  );
        //printf( " %s \n", errLog );
    }
}

void Shader::compileShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint geomShader, GLuint& shader, char*& errLog ){
    shaderprogram = glCreateProgram();
    glAttachShader(shaderprogram, vertexshader );
    glAttachShader(shaderprogram, fragmentshader );
    if(geomShader) glAttachShader(shaderprogram, geomShader );
    glBindAttribLocation(shaderprogram, 0, "in_Position");
    //glBindAttribLocation(shaderprogram, 1, "in_Color");
    glLinkProgram(shaderprogram);
    int isLinked;
    glGetProgramiv(shaderprogram, GL_LINK_STATUS, (int *)&isLinked);
    if( isLinked == false)    {
        int maxLength;
        glGetProgramiv(shaderprogram, GL_INFO_LOG_LENGTH, &maxLength);  			 			// Noticed that glGetProgramiv is used to get the length for a shader program, not glGetShaderiv.
        errLog = (char *)malloc(maxLength);   									// The maxLength includes the NULL character
        glGetProgramInfoLog(shaderprogram, maxLength, &maxLength, errLog );		// Notice that glGetProgramInfoLog, not glGetShaderInfoLog.
        //printf( " Error linking shaderProgram : \n" );
        //printf( " %s \n", errLog );
    }
}

int Shader::init_str( const char * vertexsource , const char * fragmentsource, const char * GSsrc ){
    char * errLog = NULL;
    geomShader = 0;
    compileShader( GL_VERTEX_SHADER,   vertexsource,   vertexshader, errLog   );  checkError( errLog, " Error in Vertex Shader: \n" );
    compileShader( GL_FRAGMENT_SHADER, fragmentsource, fragmentshader, errLog  ); checkError( errLog, " Error in Fragment Shader: \n" );
    if( GSsrc  ){ compileShader( GL_GEOMETRY_SHADER, GSsrc, geomShader, errLog  ); checkError( errLog, " Error in Geometry Shader: \n" ); }
    compileShaderProgram( vertexshader, fragmentshader, geomShader, shaderprogram, errLog ); checkError( errLog, " Error in linking of Shader Program: \n" );
    return 0;
}

int Shader::init( const char * vertName, const char * fragName, const char * GSName ){
    if( GSName ){ printf( "loading shaders: Frag: %s Vert: %s Geom: %s \n", vertName, fragName, GSName ); }
    else        { printf( "loading shaders: Frag: %s Vert: %s\n",           vertName, fragName ); }
    int ret=0;
    char * vertexsource   = filetobuf( vertName );  if( vertexsource   == NULL ){ printf( "cannot load %s \n", vertName ); exit(1); }
    char * fragmentsource = filetobuf( fragName );  if( fragmentsource == NULL ){ printf( "cannot load %s \n", fragName ); exit(1); }
    char * GSsrc = NULL;
    if( GSName ){ GSsrc = filetobuf( GSName );  if( GSsrc == NULL ){ printf( "cannot load %s \n", GSName ); exit(1); } }
    ret = init_str( vertexsource , fragmentsource, GSsrc );
    delete [] vertexsource;
    delete [] fragmentsource;
    return ret;
}

void Shader::destory(){
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

void Shader::getDefaultUniformLocation(){
    uloc_camMat    = glGetUniformLocation( shaderprogram, "camMat"    );
    uloc_camPos    = glGetUniformLocation( shaderprogram, "camPos"    );
    uloc_modelMat  = glGetUniformLocation( shaderprogram, "modelMat"  );
    uloc_modelPos  = glGetUniformLocation( shaderprogram, "modelPos"  ); // glUniform3fv      (uloc, 1, modelPos );
    uloc_baseColor = glGetUniformLocation( shaderprogram, "baseColor" );
}

//void Shader::setUniformi    ( const char * name, const int    i)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform1iv      (u, 1, (GLint*  )&i           ); };
void Shader::setUniformi    ( const char * name, const int    i)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform1i      (u,  i ); };
void Shader::setUniformf    ( const char * name, const float  f)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform1fv      (u, 1, (GLfloat*)&f           ); };
void Shader::setUniformVec2f( const char * name, const Vec2f  v)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform2fv      (u, 1, (GLfloat*)&v           ); };
void Shader::setUniformVec3f( const char * name, const Vec3f  v)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform3fv      (u, 1, (GLfloat*)&v           ); };
void Shader::setUniformVec4f( const char * name, const Quat4f v)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform4fv      (u, 1, (GLfloat*)&v           ); };
void Shader::setUniformMat3f( const char * name, const Mat3f  m)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix3fv(u, 1, GL_FALSE, (GLfloat*)&m ); };
void Shader::setUniformMat4f( const char * name, const Mat4f  m)const{ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix4fv(u, 1, GL_FALSE, (GLfloat*)&m ); };

