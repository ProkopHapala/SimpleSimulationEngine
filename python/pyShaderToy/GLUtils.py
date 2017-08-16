#!/usr/bin/python

import OpenGL.GL as gl
import OpenGL.arrays.vbo as glvbo

# Vertex shader
DEFAULT_VERTEX_SHADER = """
#version 330
layout(location = 0) in vec2 pos;
void main(){
    gl_Position = vec4(pos.x, pos.y, 0.0, 1.0);
}
"""

# Fragment shader
DEFAULT_FRAGMENT_SHADER = """
#version 330
out vec4 out_color;
uniform vec2  iResolution;
uniform float iTime;
void main(){
    vec2 fuv  = gl_FragCoord.xy/iResolution;
    out_color = vec4(sin(iTime*fuv), sin(0.5*iTime*fuv.x)-cos(0.5*iTime*fuv.y), 1.0);
}
"""

def compile_vertex_shader(source):
    """Compile a vertex shader from source."""
    vertex_shader = gl.glCreateShader(gl.GL_VERTEX_SHADER)
    gl.glShaderSource(vertex_shader, source)
    gl.glCompileShader(vertex_shader)
    # check compilation error
    result = gl.glGetShaderiv(vertex_shader, gl.GL_COMPILE_STATUS)
    if not(result):
        raise RuntimeError(gl.glGetShaderInfoLog(vertex_shader))
        exit()
    return vertex_shader

def compile_fragment_shader(source):
    """Compile a fragment shader from source."""
    fragment_shader = gl.glCreateShader(gl.GL_FRAGMENT_SHADER)
    gl.glShaderSource(fragment_shader, source)
    gl.glCompileShader(fragment_shader)
    # check compilation error
    result = gl.glGetShaderiv(fragment_shader, gl.GL_COMPILE_STATUS)
    if not(result):
        raise RuntimeError(gl.glGetShaderInfoLog(fragment_shader))
        exit()
    return fragment_shader

def link_shader_program(vertex_shader, fragment_shader):
    """Create a shader program with from compiled shaders."""
    program = gl.glCreateProgram()
    gl.glAttachShader(program, vertex_shader)
    gl.glAttachShader(program, fragment_shader)
    gl.glLinkProgram(program)
    # check linking error
    result = gl.glGetProgramiv(program, gl.GL_LINK_STATUS)
    if not(result):
        raise RuntimeError(gl.glGetProgramInfoLog(program))
        exit()
    return program


def fromShaderToy(vertex_code):
    vertex_code_ = '''
    uniform vec3      iResolution;           // viewport resolution (in pixels)
    uniform float     iTime;                 // shader playback time (in seconds)
    uniform float     iTimeDelta;            // render time (in seconds)
    uniform int       iFrame;                // shader playback frame
    uniform float     iChannelTime[4];       // channel playback time (in seconds)
    uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
    uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
    //uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
    uniform vec4      iDate;                 // (year, month, day, time in seconds)
    uniform float     iSampleRate;           // sound sample rate (i.e., 44100)
    ''' + vertex_code + '''
    void main(){
        mainImage( gl_FragColor, gl_FragCoord );
    };
    '''
    return vertex_code_
