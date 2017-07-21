
#include "IO_utils.h"

#include "Shader.h" // THE HEADER

void Shader::compileShader( GLenum shaderType, const char* sourceCode, GLuint& shader, char*& errLog ){
    shader = glCreateShader( shaderType );   // Create an empty vertex shader handle
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

void Shader::compileShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint& shader, char*& errLog ){
    shaderprogram = glCreateProgram();
    glAttachShader(shaderprogram, vertexshader );
    glAttachShader(shaderprogram, fragmentshader );
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

int Shader::init_str( const char * vertexsource , const char * fragmentsource ){
    //printf( " DEBUG 1 \n" );
    char * errLog = NULL;
    compileShader( GL_VERTEX_SHADER,   vertexsource,   vertexshader, errLog   );
    if( errLog != NULL ){
        printf( " Error in Vertex Shader %s : \n" );
        printf( " %s \n", errLog );
        free( errLog );
        return -1;
    }
    compileShader( GL_FRAGMENT_SHADER, fragmentsource, fragmentshader, errLog  );
    if( errLog != NULL ){
        printf( " Error in Fragment Shader %s : \n" );
        printf( " %s \n", errLog );
        free( errLog );
        return -2;
    }
    compileShaderProgram( vertexshader, fragmentshader, shaderprogram, errLog );
    if( errLog != NULL ){
        printf( " Error in linking of Shader Program : \n" );
        printf( " %s \n", errLog );
        free( errLog );
        return -3;
    }
}

int Shader::init( const char * vertName, const char * fragName ){
    printf( "loading shaders: Frag: %s Vert: %s\n", vertName, fragName );
    vertexsource   = filetobuf( vertName );  if( vertexsource   == NULL ){ printf( "cannot load %s \n", vertName ); exit(1); }
    fragmentsource = filetobuf( fragName );  if( fragmentsource == NULL ){ printf( "cannot load %s \n", fragName ); exit(1); }
    return init_str( vertexsource , fragmentsource );
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

void Shader::setUniformi    ( char * name, int    i){ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform1iv      (u, 1, (GLint*  )&i           ); };
void Shader::setUniformf    ( char * name, float  f){ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform1fv      (u, 1, (GLfloat*)&f           ); };
void Shader::setUniformVec2f( char * name, Vec2f  v){ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform2fv      (u, 1, (GLfloat*)&v           ); };
void Shader::setUniformVec3f( char * name, Vec3f  v){ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform3fv      (u, 1, (GLfloat*)&v           ); };
void Shader::setUniformVec4f( char * name, Quat4f v){ GLuint u = glGetUniformLocation(shaderprogram,name); glUniform4fv      (u, 1, (GLfloat*)&v           ); };
void Shader::setUniformMat3f( char * name, Mat3f  m){ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix3fv(u, 1, GL_FALSE, (GLfloat*)&m ); };
void Shader::setUniformMat4f( char * name, Mat4f  m){ GLuint u = glGetUniformLocation(shaderprogram,name); glUniformMatrix4fv(u, 1, GL_FALSE, (GLfloat*)&m ); };

