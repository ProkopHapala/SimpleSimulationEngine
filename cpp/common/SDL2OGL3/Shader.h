
class Shader{ 
	public:
	GLchar *vertexsource;
	GLchar *fragmentsource;
	GLuint vertexshader;
	GLuint fragmentshader; 
	GLuint shaderprogram;

	//Shader(){}

	void compileShader( GLenum shaderType, char* sourceCode, GLuint& shader, char*& errLog ){
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

	void compileShaderProgram( GLuint vertexshader, GLuint fragmentshader, GLuint& shader, char*& errLog ){
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


	void init( ){

		// Read our shaders into the appropriate buffers 
		char const * vertName = "shaders/vert_1.c";
		char const * fragName = "shaders/frag_1.c";
		vertexsource   = filetobuf( vertName );
		fragmentsource = filetobuf( fragName );
	 
		char * errLog = NULL;

		compileShader( GL_VERTEX_SHADER,   vertexsource,   vertexshader, errLog   );
		if( errLog != NULL ){
			printf( " Error in compilation of shader %s : \n", vertName );
			printf( " %s \n", errLog );
			free( errLog );
			return;
		}

		compileShader( GL_FRAGMENT_SHADER, fragmentsource, fragmentshader, errLog  );
		if( errLog != NULL ){
			printf( " Error in compilation of shader %s : \n", fragName );
			printf( " %s \n", errLog );
			free( errLog );
			return;
		}

		compileShaderProgram( vertexshader, fragmentshader, shaderprogram, errLog );
		if( errLog != NULL ){
			printf( " Error in linking of shader program : \n" );
			printf( " %s \n", errLog );
			free( errLog );
			return;
		}

		//free(vertexsource);
		//free(fragmentsource);
	}

	void destory(){
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

};
