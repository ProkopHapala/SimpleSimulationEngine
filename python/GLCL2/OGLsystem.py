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
        # --- DEBUG registries
        self.textures = {}   # name -> {id:int, size:(w,h)}
        self.fbos     = {}   # name -> fbo id
        # Fullscreen quad geometry (lazy init)
        self._quad_vao = 0
        self._quad_vbo = 0

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

    # ------------------------------------------------------------------
    # Textures & FBOs
    # ------------------------------------------------------------------
    def create_texture_2d(self, name, size, *, internal_format=GL_RGBA32F, fmt=GL_RGBA, typ=GL_FLOAT, filter=GL_LINEAR, wrap=GL_CLAMP_TO_EDGE):
        if name in self.textures:
            print(f"WARNING: texture '{name}' already exists; re-creating")
        w, h = size
        tex = glGenTextures(1)
        if tex == 0: raise RuntimeError("glGenTextures failed")
        glBindTexture(GL_TEXTURE_2D, tex)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,     wrap)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,     wrap)
        glTexImage2D(GL_TEXTURE_2D, 0, internal_format, w, h, 0, fmt, typ, None)
        glBindTexture(GL_TEXTURE_2D, 0)
        if check_gl_error("create_texture_2d"): raise RuntimeError("OpenGL error in create_texture_2d")
        self.textures[name] = {"id": tex, "size": (w, h)}
        print(f"create_texture_2d('{name}') id={tex} size=({w},{h})")
        return tex

    def get_texture(self, name):
        t = self.textures.get(name)
        if t is None: raise KeyError(f"Texture '{name}' not found")
        return t["id"]

    def bind_texture_unit(self, name, unit):
        tex = self.get_texture(name)
        glActiveTexture(GL_TEXTURE0 + unit)
        glBindTexture(GL_TEXTURE_2D, tex)
        if check_gl_error(f"bind_texture_unit('{name}',{unit})"): raise RuntimeError("OpenGL error in bind_texture_unit")

    def create_fbo(self, name, color_tex_name):
        if name in self.fbos:
            print(f"WARNING: fbo '{name}' already exists; re-creating")
        tex = self.get_texture(color_tex_name)
        fbo = glGenFramebuffers(1)
        if fbo == 0: raise RuntimeError("glGenFramebuffers failed")
        glBindFramebuffer(GL_FRAMEBUFFER, fbo)
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex, 0)
        status = glCheckFramebufferStatus(GL_FRAMEBUFFER)
        if status != GL_FRAMEBUFFER_COMPLETE:
            glBindFramebuffer(GL_FRAMEBUFFER, 0)
            raise RuntimeError(f"Framebuffer incomplete: status=0x{status:04X}")
        glBindFramebuffer(GL_FRAMEBUFFER, 0)
        if check_gl_error("create_fbo"): raise RuntimeError("OpenGL error in create_fbo")
        self.fbos[name] = fbo
        print(f"create_fbo('{name}') id={fbo} color={color_tex_name}")
        return fbo

    def get_fbo(self, name):
        f = self.fbos.get(name)
        if f is None: raise KeyError(f"FBO '{name}' not found")
        return f

    # ------------------------------------------------------------------
    # Fullscreen quad utilities
    # ------------------------------------------------------------------
    def _init_fullscreen_quad(self):
        if self._quad_vao != 0: return
        # 4 vertices for triangle strip: (-1,-1),(1,-1),(-1,1),(1,1)
        verts = np.array([
            -1.0, -1.0,
             1.0, -1.0,
            -1.0,  1.0,
             1.0,  1.0,
        ], dtype=np.float32)
        self._quad_vao = glGenVertexArrays(1)
        self._quad_vbo = glGenBuffers(1)
        if self._quad_vao == 0 or self._quad_vbo == 0: raise RuntimeError("Failed to create quad VAO/VBO")
        glBindVertexArray(self._quad_vao)
        glBindBuffer(GL_ARRAY_BUFFER, self._quad_vbo)
        glBufferData(GL_ARRAY_BUFFER, verts.nbytes, verts, GL_STATIC_DRAW)
        glEnableVertexAttribArray(0)
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, ctypes.c_void_p(0))
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        glBindVertexArray(0)
        if check_gl_error("_init_fullscreen_quad"): raise RuntimeError("OpenGL error in _init_fullscreen_quad")
        print(f"Fullscreen quad initialized vao={self._quad_vao} vbo={self._quad_vbo}")

    def draw_fullscreen_quad(self):
        if self._quad_vao == 0: self._init_fullscreen_quad()
        glBindVertexArray(self._quad_vao)
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4)
        glBindVertexArray(0)
        if check_gl_error("draw_fullscreen_quad"): raise RuntimeError("OpenGL error in draw_fullscreen_quad")

    def render_fs_to_fbo(self, program, out_fbo_name, bind_textures=None):
        """Render currently bound FS program to named FBO using the fullscreen quad.
        bind_textures: list of texture names bound to texture units 0..N-1.
        """
        if bind_textures is None: bind_textures = []
        fbo = self.get_fbo(out_fbo_name)
        for i, tname in enumerate(bind_textures): self.bind_texture_unit(tname, i)
        # Preserve the currently bound draw framebuffer (QOpenGLWidget uses an internal FBO)
        prev_fb = glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING)
        glBindFramebuffer(GL_FRAMEBUFFER, fbo)
        glUseProgram(program)
        self.draw_fullscreen_quad()
        glUseProgram(0)
        # Restore the previous framebuffer to remain compatible with QOpenGLWidget
        glBindFramebuffer(GL_FRAMEBUFFER, int(prev_fb))
        if check_gl_error("render_fs_to_fbo"): raise RuntimeError("OpenGL error in render_fs_to_fbo")


# ----------------------------------------------------------------------
# OpenGL Compute Shader Helper
# ----------------------------------------------------------------------
def compile_compute_program(compute_src):
    if not compute_src: raise ValueError("Compute shader source code cannot be empty.")
    cs = glCreateShader(GL_COMPUTE_SHADER)
    if not cs: raise RuntimeError("glCreateShader(GL_COMPUTE_SHADER) failed.")
    glShaderSource(cs, compute_src)
    glCompileShader(cs)
    if glGetShaderiv(cs, GL_COMPILE_STATUS) != GL_TRUE:
        log = glGetShaderInfoLog(cs).decode()
        glDeleteShader(cs)
        raise RuntimeError(f"Compute shader compilation failed:\n{log}")
    prog = glCreateProgram()
    if not prog: raise RuntimeError("glCreateProgram() failed for compute shader.")
    glAttachShader(prog, cs)
    glLinkProgram(prog)
    if glGetProgramiv(prog, GL_LINK_STATUS) != GL_TRUE:
        log = glGetProgramInfoLog(prog).decode()
        glDeleteProgram(prog)
        raise RuntimeError(f"Compute program linking failed:\n{log}")
    glDeleteShader(cs)
    print(f"compile_compute_program() OK program={prog}")
    return prog


class GLComputeProgram:
    def __init__(self, src=None):
        self.program = None
        if src is not None:
            self.compile(src)

    def compile(self, src):
        self.program = compile_compute_program(src)
        return self.program

    def use(self):
        if self.program is None: raise RuntimeError("GLComputeProgram.use(): program is None")
        glUseProgram(self.program)

    def set1f(self, name, v):
        if self.program is None: raise RuntimeError("set1f(): program is None")
        loc = glGetUniformLocation(self.program, name)
        if loc != -1: glUniform1f(loc, v)

    def set1i(self, name, v):
        if self.program is None: raise RuntimeError("set1i(): program is None")
        loc = glGetUniformLocation(self.program, name)
        if loc != -1: glUniform1i(loc, v)

    def set3f(self, name, x, y, z):
        if self.program is None: raise RuntimeError("set3f(): program is None")
        loc = glGetUniformLocation(self.program, name)
        if loc != -1: glUniform3f(loc, x, y, z)

    def bind_ssbo(self, binding, buffer_id):
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, binding, buffer_id)

    def dispatch(self, gx, gy=1, gz=1, barrier_bits=GL_SHADER_STORAGE_BARRIER_BIT):
        if self.program is None: raise RuntimeError("dispatch(): program is None")
        glUseProgram(self.program)
        glDispatchCompute(int(gx), int(gy), int(gz))
        glMemoryBarrier(barrier_bits)
        glUseProgram(0)
        if check_gl_error("GLComputeProgram.dispatch"): raise RuntimeError("OpenGL error in GLComputeProgram.dispatch")