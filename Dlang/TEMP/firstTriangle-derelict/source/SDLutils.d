module glUtils;

public import std.stdio;
public import std.format;

public import derelict.sdl2.sdl;
public import derelict.opengl3.gl3;
public import derelict.util.exception;

public import glFunctions;
public import glObjects;
public import Shader;

public import Vector;

import core.stdc.stdlib;
//import std.conv;

//void double2float( int n, const double[] ds, float[] fs ){ for(int i=0; i<n; i++){ fs[i]=cast(float)ds[i]; }; }

// ======== Compile Shader

class WindowSDL{
    SDL_Window*   window;
    SDL_GLContext context;

    this( const char* wname, int sw, int sh ){
        // Create a centered 640x480 OpenGL window named "Triangle"
        //writeln( "DEBUG 3.1 ");
        //writeln( format("wname: %s", wname ) );
        window = SDL_CreateWindow( wname, SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, sw, sh, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN );
        //SDL_Window*   window = SDL_CreateWindow( "Triangle", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, sw, sh, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN );
        //SDL_Window* window = SDL_CreateWindow("Triangle",SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,640, 480,SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
        //writeln( "DEBUG 3.2 ");
        if(null is window){ writeln("Failed to create the application window"); exit(1); } // Exit if window creation fails.
        //writeln( "DEBUG 3.3 ");
        context = SDL_GL_CreateContext(window);  // Create an OpenGL context for our window.
        //writeln( "DEBUG 3.4 ");
        DerelictGL3.reload();   // Load all OpenGL functions and extensions supported by Derelict.
        //writeln( "DEBUG 3.5 ");
    }

    ~this(){
        SDL_GL_DeleteContext(context);
        SDL_DestroyWindow(window);
    }

};

GLuint makeVertexArrays(){
    GLuint vertexArrayID;
    glGenVertexArrays(1, &vertexArrayID);
    glBindVertexArray(vertexArrayID);
    return vertexArrayID;
}

GLuint makeVertexBuffer( Vec3f[] vertices ){
    GLuint vertexBufferID;
    glGenBuffers(1, &vertexBufferID);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBufferID);
    glBufferData(GL_ARRAY_BUFFER, vertices.length * Vec3f.sizeof, vertices.ptr, GL_STATIC_DRAW );
    return vertexBufferID;
}

void initDerelictSDL2GL3(){
    // Load SDL2 and GL.
    try{
        //DerelictSDL2.load();
        //DerelictGL3.load();
        // https://stackoverflow.com/questions/35778861/segmentation-fault-when-using-derelict-sdl
        DerelictSDL2.load(SharedLibVersion(2, 0, 2));
        DerelictGL3.load();
    }
    // Print errors, if any.
    catch(SharedLibLoadException e){ writeln("SDL2 or GL not found: " ~ e.msg); }
    catch(SymbolLoadException    e){ writeln("Missing SDL2 or GL symbol (old version installed?): " ~ e.msg); }

    // Initialize SDL Video subsystem.
    if(SDL_Init(SDL_INIT_VIDEO) < 0){
        // SDL_Init returns a negative number on error.
        writeln("SDL Video subsystem failed to initialize");
        //return 1;
        exit(1);
    }

    // We want OpenGL 3.2 (change this if you want a newer version),
    // and the core profile (i.e. no deprecated functions)
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);

    // 32bit RGBA window
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE,     8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,   8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,    8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,   8);
    // Double buffering to avoid tearing
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    // Depth buffer. Not useful when drawing a triangle, but almost always
    // useful when drawing 3D
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,   24);
}

void closeDerelictSDL2GL3(){
    SDL_Quit();
    DerelictSDL2.unload();
    DerelictGL3.unload();
}


GLuint compileShaders(string vertexSrc, string fragmentSrc){
    // Create the shaders
    GLuint vertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

    import std.string;
    // Compile the vertex shader
    writeln("we're about to attempt compiling a vertex shader");
    auto vertexZeroTerminated = toStringz(vertexSrc);
    glShaderSource(vertexShaderID, 1, &vertexZeroTerminated, null);
    glCompileShader(vertexShaderID);

    // Use this to determine how much to allocate if infoLog is too short
    // glGetShaderiv(vertexShaderID, GL_INFO_LOG_LENGTH, &infoLogLength);

    // Check for errors
    GLint compiled;
    glGetShaderiv(vertexShaderID, GL_COMPILE_STATUS, &compiled);
    char[1024 * 8] infoLog;
    glGetShaderInfoLog(vertexShaderID, infoLog.length, null, infoLog.ptr);
    import core.stdc.stdio;
    writeln("vertex shader info log:");
    puts(infoLog.ptr);
    if(!compiled){
        throw new Exception("Failed to compile vertex shader " ~ vertexSrc);
    }

    // Compile Fragment Shader
    writeln("we're about to attempt compiling a fragment shader");
    auto fragmentZeroTerminated = toStringz(fragmentSrc);
    glShaderSource(fragmentShaderID, 1, &fragmentZeroTerminated, null);
    glCompileShader(fragmentShaderID);

    // Check for errors
    glGetShaderiv(fragmentShaderID, GL_COMPILE_STATUS, &compiled);
    glGetShaderInfoLog(fragmentShaderID, infoLog.length, null, infoLog.ptr);
    writeln("fragment shader info log:");
    puts(infoLog.ptr);
    if(!compiled){
        throw new Exception("Failed to compile fragment shader " ~ fragmentSrc);
    }

    // Link the program
    writeln("we're about to attempt linking");
    GLuint programID = glCreateProgram();
    glAttachShader(programID, vertexShaderID);
    glAttachShader(programID, fragmentShaderID);
    glLinkProgram(programID);

    // Check the program
    GLint linked;
    glGetProgramiv(programID, GL_LINK_STATUS, &linked);
    glGetProgramInfoLog(programID, infoLog.length, null, infoLog.ptr);
    writeln("linking info log:");
    puts(infoLog.ptr);
    if(!linked){
        throw new Exception("Failed to link shaders " ~ vertexSrc ~ " " ~ fragmentSrc);
    }

    glDeleteShader(vertexShaderID);
    glDeleteShader(fragmentShaderID);

    return programID;
}
