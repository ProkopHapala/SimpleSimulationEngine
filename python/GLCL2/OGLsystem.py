import numpy as np
import os
from OpenGL.GL import *
import ctypes

def check_gl_error(msg=""):
    err = glGetError()
    if err != GL_NO_ERROR:
        print(f"OpenGL Error @ {msg}: {err}")
        return True
    return False

def compile_shader_program(vertex_src, fragment_src):
    if not vertex_src or not fragment_src:
        raise ValueError("Shader source code cannot be empty.")
    
    vertex_shader = glCreateShader(GL_VERTEX_SHADER)
    if not vertex_shader: raise RuntimeError("glCreateShader(GL_VERTEX_SHADER) failed.")
    glShaderSource(vertex_shader, vertex_src)
    glCompileShader(vertex_shader)
    if glGetShaderiv(vertex_shader, GL_COMPILE_STATUS) != GL_TRUE:
        log = glGetShaderInfoLog(vertex_shader).decode()
        glDeleteShader(vertex_shader)
        raise RuntimeError(f"Vertex shader compilation failed:\n{log}")

    fragment_shader = glCreateShader(GL_FRAGMENT_SHADER)
    if not fragment_shader: raise RuntimeError("glCreateShader(GL_FRAGMENT_SHADER) failed.")
    glShaderSource(fragment_shader, fragment_src)
    glCompileShader(fragment_shader)
    if glGetShaderiv(fragment_shader, GL_COMPILE_STATUS) != GL_TRUE:
        log = glGetShaderInfoLog(fragment_shader).decode()
        glDeleteShader(vertex_shader)
        glDeleteShader(fragment_shader)
        raise RuntimeError(f"Fragment shader compilation failed:\n{log}")

    program = glCreateProgram()
    if not program: raise RuntimeError("glCreateProgram() failed.")
    glAttachShader(program, vertex_shader)
    glAttachShader(program, fragment_shader)
    glLinkProgram(program)
    if glGetProgramiv(program, GL_LINK_STATUS) != GL_TRUE:
        log = glGetProgramInfoLog(program).decode()
        glDeleteProgram(program)
        raise RuntimeError(f"Shader program linking failed:\n{log}")

    glDeleteShader(vertex_shader)
    glDeleteShader(fragment_shader)
    return program

class GLobject:
    def __init__(self, vao=None, vbo=None, ebo=None, nelements=0, mode=GL_TRIANGLES):
        self.vao, self.vbo, self.ebo = vao, vbo, ebo
        self.nelements = nelements
        self.mode = mode

    def alloc_vao_vbo_ebo(self, components):
        self.vao = glGenVertexArrays(1)
        glBindVertexArray(self.vao)
        self.vbo = glGenBuffers(1)
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        
        total_stride = sum(components) * ctypes.sizeof(GLfloat)
        offset = 0
        for i, n_comp in enumerate(components):
            glEnableVertexAttribArray(i)
            glVertexAttribPointer(i, n_comp, GL_FLOAT, GL_FALSE, total_stride, ctypes.c_void_p(offset))
            offset += n_comp * ctypes.sizeof(GLfloat)
            
        glBindVertexArray(0)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        return self.vao, self.vbo, None

    def draw_arrays(self, n=-1):
        if self.vao is None: return
        glBindVertexArray(self.vao)
        glDrawArrays(self.mode, 0, self.nelements if n == -1 else n)
        glBindVertexArray(0)

    def upload_vbo(self, vertices):
        if self.vbo is None: return
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        # Use glBufferSubData for updates if the size is the same, glBufferData for initial/resize
        glBufferData(GL_ARRAY_BUFFER, vertices.nbytes, vertices, GL_DYNAMIC_DRAW)
        glBindBuffer(GL_ARRAY_BUFFER, 0)

class OGLSystem:
    def __init__(self):
        self.shader_programs = {}
        self.shader_configs = {}
        self.script_dir = "."

    def compile_all_shaders(self):
        print("Compiling all configured shaders...")
        for name, (vert_path, frag_path, _) in self.shader_configs.items():
            try:
                vert_full_path = os.path.join(self.script_dir, vert_path)
                frag_full_path = os.path.join(self.script_dir, frag_path)
                with open(vert_full_path, 'r') as f: vert_src = f.read()
                with open(frag_full_path, 'r') as f: frag_src = f.read()
                
                program = compile_shader_program(vert_src, frag_src)
                self.shader_programs[name] = program
                print(f"  Successfully compiled and linked shader '{name}' (Program ID: {program})")
            except (IOError, RuntimeError) as e:
                print(f"ERROR: Failed to load or compile shader '{name}':\n{e}")

    def get_shader_program(self, name):
        return self.shader_programs.get(name)
        
    def clear(self):
        print("OGLSystem::clear() Deleting all shader programs.")
        for program in self.shader_programs.values():
            glDeleteProgram(program)
        self.shader_programs.clear()