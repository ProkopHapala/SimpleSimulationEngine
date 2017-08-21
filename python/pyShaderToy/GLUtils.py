#!/usr/bin/python

import OpenGL.GL as gl
import OpenGL.arrays.vbo as glvbo

# Vertex shader
DEFAULT_VERTEX_SHADER = """
#version 100
//in  vec2 in_Position;
attribute vec2 pos;
void main(){
    gl_Position = vec4(pos.x, pos.y, 0.0, 1.0);
}
"""

# Fragment shader
DEFAULT_FRAGMENT_SHADER = """
#version 100
//uniform sampler2D tex0;
//varying vec2 vTexCoord;
void main() {
    //vec4 color = texture2D(tex0, vTexCoord);
    //gl_FragColor = color;
    //gl_FragColor = vec4(1.0,0.0,0.0,0.0);
    //gl_FragColor = vec4( sin(50.0*gl_FragCoord.xyz), 0.0 );
    //gl_FragColor = vec4( sin(50.0*gl_FragCoord.xyz), 0.0 );
}
"""

# Vertex shader
DEFAULT_VERTEX_SHADER_330 = """
#version 330
layout(location = 0) in vec2 pos;
void main(){
    gl_Position = vec4(pos.x, pos.y, 0.0, 1.0);
}
"""

# Fragment shader
DEFAULT_FRAGMENT_SHADER_330 = """
#version 330
out vec4 out_color;
uniform vec2  iResolution;
uniform float iTime;
void main(){
    vec2 fuv  = gl_FragCoord.xy/iResolution;
    out_color = vec4(sin(iTime*fuv), sin(0.5*iTime*fuv.x)-cos(0.5*iTime*fuv.y), 1.0);
}
"""


SHADERTOY_HEADER = '''
#version 100
precision highp float;
uniform vec3      iResolution;           // viewport resolution (in pixels)
uniform float     iTime;                 // shader playback time (in seconds)
uniform float     iTimeDelta;            // render time (in seconds)
uniform int       iFrame;                // shader playback frame
uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
uniform vec4      iDate;                 // (year, month, day, time in seconds)
//uniform float     iChannelTime[4];       // channel playback time (in seconds)
//uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
//uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
//uniform float     iSampleRate;           // sound sample rate (i.e., 44100)
'''

def getShaderToyUrl(id):
    import urllib, urllib2, json
    url = 'https://www.shadertoy.com/shadertoy'
    headers = { 'Referer' : 'https://www.shadertoy.com/' }
    values  = { 's' : json.dumps ({'shaders' : [id]}) }
    data = urllib.urlencode (values).encode ('utf-8')
    req  = urllib2.Request (url, data, headers)
    response   = urllib2.urlopen (req)
    sh_str = response.read().decode ('utf-8')
    return sh_str

def downloadFromShaderToy(id, save_json=None):
    import json
    shs = getShaderToyUrl(id)
    shj = json.loads(shs)
    if save_json:
        with open(save_json, "w") as f: f.write(json.dumps(shj, indent=2))
    assert (len(shj) == 1)
    shj = shj[0]
    name  = shj["info"]["name"]
    desc  = shj["info"]["description"]
    codes = [ p["code"] for p in shj["renderpass"] ]
    return name, desc, codes

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
    vertex_code_ =  SHADERTOY_HEADER + vertex_code + "\n\nvoid main(){mainImage( gl_FragColor, gl_FragCoord.xy );}"
    return vertex_code_
