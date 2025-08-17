import sys
import numpy as np
import os

from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D, QSurfaceFormat
from OpenGL.GL import *

from .OGLsystem import GLobject, compile_shader_program

class GLCLWidget(QOpenGLWidget):
    @staticmethod
    def apply_default_surface_format(swap_interval=None):
        """Set the default QSurfaceFormat before any GL context is created.
        Call this early (CLI entry) to control vsync via swap interval.
        """
        fmt = QSurfaceFormat()
        if swap_interval is not None:
            fmt.setSwapInterval(int(swap_interval))
        QSurfaceFormat.setDefaultFormat(fmt)

    def __init__(self, parent=None, enable_opengl_debug=False):
        super().__init__(parent)
        self.browser = parent
        self.ogl_system = None
        self.ocl_system = None
        self.enable_opengl_debug = enable_opengl_debug

        self.view_matrix = QMatrix4x4()
        self.projection_matrix = QMatrix4x4()
        
        self.camera_pos = QVector3D(0, 0, 5)
        self.last_mouse_pos = None
        self.cam_rot = QMatrix4x4()

        self.gl_objects = {}
        self.render_pipeline_info = []
        self.buffer_data = {}
        # FS pipeline state
        self.fs_textures_cfg = {}
        self.fs_fbos_cfg = {}
        self.fs_pipeline = []
        # Frame counter for FS shaders expecting iFrame
        self.frame_counter = 0
        # Viewer overlay state
        self.display_tex_name = None  # name of texture to show on default framebuffer
        self._viewer_program = 0      # lazily compiled simple blit shader
        # Cache for uniform locations per program ID to avoid per-frame glGetUniformLocation
        self._uniform_cache = {}      # { program_id: { 'projection':loc, 'view':loc } }

        self.setFocusPolicy(Qt.StrongFocus)

    def set_systems(self, ogl_system, ocl_system, browser):
        self.ogl_system = ogl_system
        self.ocl_system = ocl_system
        self.browser = browser

    def set_render_config(self, buffer_data, render_pipeline_info):
        """Receive rendering configuration from the browser."""
        self.buffer_data = buffer_data
        self.render_pipeline_info = render_pipeline_info
        # Defer resource creation until initializeGL is called
        if self.isValid():
            self.bake_render_objects()
            # Apply geometry uniforms for new/changed programs
            try:
                self.makeCurrent()
                self.apply_viewproj_for_geometry()
            finally:
                self.doneCurrent()
            self.update()

    def set_fs_config(self, textures_cfg, fbos_cfg, fs_pipeline):
        """Optional FS pipeline config.
        textures_cfg: {name:(w,h)}; fbos_cfg: {fbo_name: color_tex_name}; fs_pipeline: [(shader_name, out_fbo, [in_textures])]
        """
        self.fs_textures_cfg = dict(textures_cfg or {})
        self.fs_fbos_cfg = dict(fbos_cfg or {})
        self.fs_pipeline = list(fs_pipeline or [])
        if self.isValid():
            self._create_fs_resources()
            # Apply FS static uniforms for new/changed passes
            try:
                self.makeCurrent()
                self.apply_fs_static_uniforms()
            finally:
                self.doneCurrent()
            self.update()

    def upload_fs_texture_data(self, textures_dict):
        """Upload pixel data arrays into FS textures. textures_dict: {name: np.ndarray(H,W,[1|3|4])}
        Safe to call anytime; will no-op if context not yet valid.
        """
        if not textures_dict: return
        if not self.isValid():
            # Defer: resources not yet created; caller should call again after rebuild
            return
        self.makeCurrent()
        try:
            for name, arr in textures_dict.items():
                if name not in self.fs_textures_cfg:
                    print(f"upload_fs_texture_data: WARNING texture '{name}' not configured; skipping")
                    continue
                try:
                    self.ogl_system.upload_texture_2d_data(name, arr)
                except Exception as e:
                    print(f"upload_fs_texture_data: ERROR uploading '{name}': {e}")
        finally:
            self.doneCurrent()

    def initializeGL(self):
        try:
            print("GLCLWidget::initializeGL() - OpenGL Context is now valid.")
            glClearColor(0.1, 0.1, 0.15, 1.0)
            glEnable(GL_DEPTH_TEST)
            glEnable(GL_PROGRAM_POINT_SIZE)
            
            # Now it's safe to compile shaders and create GL resources
            self.ogl_system.compile_all_shaders()
            self.bake_render_objects()
            self._create_fs_resources()
            # Apply static uniforms once (projection/view, samplers, resolutions)
            self.apply_static_uniforms_all()
            # Upload any initial FS texture data provided by the browser (from script init())
            try:
                if self.browser is not None and getattr(self.browser, "initial_textures", None):
                    self.upload_fs_texture_data(self.browser.initial_textures)
            except Exception as e:
                print(f"initializeGL: WARNING could not upload initial FS textures: {e}")
        except Exception as e:
            if self.browser is not None:
                self.browser.on_exception(e)
            else:
                import traceback, os
                traceback.print_exc()
                os._exit(1)

    def rebuild_gl_resources(self):
        """Rebuild GL pipeline after OGLSystem.clear() on script reload.
        Requires a valid context; does nothing if context not yet initialized.
        """
        if not self.isValid():
            return
        self.makeCurrent()
        # NOTE: crash on error (no try/except); ensure context is released
        try:
            self.ogl_system.compile_all_shaders()
            self.bake_render_objects()
            self._create_fs_resources()
            # Re-apply static uniforms after relink/recreation
            self.apply_static_uniforms_all()
        finally:
            self.doneCurrent()
        # Clear uniform cache because program IDs/locations may have changed after rebuild
        self._uniform_cache.clear()
        self.update()

    def _create_fs_resources(self):
        if not self.ogl_system: return
        # Create textures
        for tname, (w, h) in self.fs_textures_cfg.items():
            self.ogl_system.create_texture_2d(tname, (int(w), int(h)))
        # Create FBOs
        for fname, color_tex in self.fs_fbos_cfg.items():
            self.ogl_system.create_fbo(fname, color_tex)

    def bake_render_objects(self):
        """Create and configure GLobject for each buffer in the render pipeline."""
        print("Baking render objects...")
        self.gl_objects.clear()

        # Create a unique set of vertex buffers that need a GLobject
        vertex_buffers_to_bake = {pass_info[2] for pass_info in self.render_pipeline_info if len(pass_info) > 2}

        # Map per-buffer geometry options collected from render pipeline optional dict (5th element)
        # options: { 'mode': 'POINTS'|'LINES'|'LINE_STRIP'|'TRIANGLES', 'attribs': [n0,n1,...] }
        geom_opts = {}
        for pass_info in self.render_pipeline_info:
            if len(pass_info) >= 5 and isinstance(pass_info[4], dict):
                buf = pass_info[2]
                # First occurrence wins to avoid conflicting configs silently; users should keep consistent
                if buf not in geom_opts:
                    geom_opts[buf] = pass_info[4]

        for buffer_name in vertex_buffers_to_bake:
            data = self.buffer_data.get(buffer_name)
            if data is None:
                print(f"Warning: No initial data for buffer '{buffer_name}'. Cannot bake GLobject.")
                continue

            nelements = len(data)
            components = data.shape[1] if data.ndim > 1 else 1

            # Resolve geometry options for this buffer
            opts = geom_opts.get(buffer_name, {})
            mode_str = str(opts.get('mode', 'POINTS')).upper() if opts else 'POINTS'
            # Map string to GL constant (defaults to GL_POINTS)
            mode_map = {
                'POINTS': GL_POINTS,
                'LINES': GL_LINES,
                'LINE_STRIP': GL_LINE_STRIP,
                'TRIANGLES': GL_TRIANGLES,
                'TRIANGLE_STRIP': GL_TRIANGLE_STRIP,
            }
            mode_gl = mode_map.get(mode_str, GL_POINTS)

            attribs = opts.get('attribs') if isinstance(opts, dict) else None
            if attribs is not None:
                try:
                    # Basic sanity: sum of attribs should match data stride
                    if int(sum(int(a) for a in attribs)) != int(components):
                        print(f"bake_render_objects: WARNING attribs {attribs} sum != components {components} for buffer '{buffer_name}'; falling back to single attrib")
                        attribs = None
                except Exception:
                    attribs = None

            gl_obj = GLobject(nelements=nelements, mode=mode_gl)
            gl_obj.alloc_vao_vbo_ebo(attribs if attribs is not None else [components])
            gl_obj.upload_vbo(data)
            
            self.gl_objects[buffer_name] = gl_obj
            print(f"  Baked GLobject for buffer '{buffer_name}' with {nelements} elements.")

    def update_buffer_data(self, buffer_name, new_data):
        """Update VBO data for a specific GLobject."""
        gl_obj = self.gl_objects.get(buffer_name)
        if gl_obj and new_data is not None:
            gl_obj.upload_vbo(new_data)

    def paintGL(self):
        try:
            # DEBUG: vivid clear to verify that we see any GL content at all
            if self.frame_counter < 60:
                glClearColor(1.0, 0.0, 1.0, 1.0)
            else:
                glClearColor(0.1, 0.1, 0.15, 1.0)
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
            # DEBUG: advance frame counter per paint; we rely on browser timer driving updates
            self.frame_counter += 1
            # Execute optional FS pipeline (offscreen)
            self._execute_fs_pipeline()
            
            # Iterate through the user-defined render pipeline
            for pass_info in self.render_pipeline_info:
                # Support optional 5th element (dict of options) in render pipeline entries
                # Expected layouts:
                #  (shader_name, elem_count, vertex_buffer, index_buffer)
                #  (shader_name, elem_count, vertex_buffer, index_buffer, {options})
                try:
                    shader_name = pass_info[0]
                    vertex_buffer_name = pass_info[2] if len(pass_info) > 2 else None
                except Exception:
                    # Malformed entry; skip
                    continue

                shader_program = self.ogl_system.get_shader_program(shader_name)
                gl_obj = self.gl_objects.get(vertex_buffer_name) if vertex_buffer_name else None

                if not shader_program or not gl_obj:
                    continue

                glUseProgram(shader_program)
                # Draw the object (arrays only; indices not yet supported)
                gl_obj.draw_arrays()

            # Optional: draw selected texture to default framebuffer as an overlay
            if self.display_tex_name:
                self._draw_display_texture()

            glUseProgram(0)
        except Exception as e:
            if self.browser is not None:
                self.browser.on_exception(e)
            else:
                import traceback, os
                traceback.print_exc()
                os._exit(1)

    def _ensure_viewer_program(self):
        if self._viewer_program:
            return
        # Minimal passthrough VS and sampler FS
        vs_src = """
        #version 330 core
        layout(location=0) in vec2 pos;
        out vec2 uv;
        void main(){
            uv = pos*0.5 + 0.5;          // map [-1,1] -> [0,1]
            gl_Position = vec4(pos,0.0,1.0);
        }
        """
        fs_src = """
        #version 330 core
        in vec2 uv;
        uniform sampler2D uTex;
        out vec4 FragColor;
        void main(){
            FragColor = texture(uTex, uv);
        }
        """
        self._viewer_program = compile_shader_program(vs_src, fs_src)
        # Initialize viewer sampler uniform once
        glUseProgram(self._viewer_program)
        loc = self._get_uniform_location(self._viewer_program, "uTex")
        if loc != -1: glUniform1i(loc, 0)
        glUseProgram(0)

    def _draw_display_texture(self):
        if self.ogl_system is None:
            return
        if self.display_tex_name not in self.ogl_system.textures:
            # Silent ignore if texture not known (e.g., before FS resources exist)
            return
        # Ensure shader exists
        self._ensure_viewer_program()
        # Preserve state
        depth_was_enabled = glIsEnabled(GL_DEPTH_TEST)
        if depth_was_enabled: glDisable(GL_DEPTH_TEST)
        # Set viewport to widget size
        glViewport(0, 0, int(self.width()), int(self.height()))
        # Bind texture to unit 0 and draw
        self.ogl_system.bind_texture_unit(self.display_tex_name, 0)
        glUseProgram(self._viewer_program)
        self.ogl_system.draw_fullscreen_quad()
        glUseProgram(0)
        # Restore state
        if depth_was_enabled: glEnable(GL_DEPTH_TEST)

    def _execute_fs_pipeline(self):
        if not self.fs_pipeline: return
        if self.ogl_system is None: return
        # Cache and restore viewport
        vp = glGetIntegerv(GL_VIEWPORT)
        # DEBUG: temporarily disable depth & scissor test for FS passes
        depth_was_enabled = glIsEnabled(GL_DEPTH_TEST)
        scissor_was_enabled = glIsEnabled(GL_SCISSOR_TEST)
        if depth_was_enabled: glDisable(GL_DEPTH_TEST)
        if scissor_was_enabled: glDisable(GL_SCISSOR_TEST)
        for (shader_name, out_fbo, bind_textures) in self.fs_pipeline:
            # if self.frame_counter < 5:
            #     try:
            #         print(f"[FS DEBUG] frame={self.frame_counter} pass='{shader_name}' out='{out_fbo}' inputs={bind_textures}")
            #     except Exception:
            #         pass
            program = self.ogl_system.get_shader_program(shader_name)
            if not program:
                raise KeyError(f"FS shader '{shader_name}' not compiled or registered")
            # Pre-compute pass resolution for viewport, but do not set uniforms here (done on events)
            pass_w = self.width()
            pass_h = self.height()
            if out_fbo and out_fbo not in ("default", "", None):
                if out_fbo not in self.fs_fbos_cfg:
                    raise KeyError(f"FS FBO '{out_fbo}' not configured")
                tex_name = self.fs_fbos_cfg.get(out_fbo)
                if tex_name and tex_name in self.fs_textures_cfg:
                    w, h = self.fs_textures_cfg[tex_name]
                    pass_w, pass_h = int(w), int(h)
            # Set only dynamic per-frame uniform
            glUseProgram(program)
            loc_frame = self._get_uniform_location(program, "iFrame")
            if loc_frame != -1:
                glUniform1i(loc_frame, int(self.frame_counter))
            glUseProgram(0)
            # Set viewport to target texture size if known
            tex_name = self.fs_fbos_cfg.get(out_fbo) if out_fbo else None
            if out_fbo and out_fbo not in ("default", "", None):
                if out_fbo not in self.fs_fbos_cfg:
                    raise KeyError(f"FS FBO '{out_fbo}' not configured")
                if tex_name and tex_name in self.fs_textures_cfg:
                    w, h = self.fs_textures_cfg[tex_name]
                    glViewport(0, 0, int(w), int(h))
            if out_fbo in ("default", "", None):
                # Render directly to default framebuffer
                glViewport(0, 0, int(self.width()), int(self.height()))
                # if self.frame_counter < 5:
                #     try:
                #         fb = glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING)
                #         vp2 = glGetIntegerv(GL_VIEWPORT)
                #         print(f"[FS DEBUG] default pass fbo={int(fb)} viewport={tuple(int(x) for x in vp2)} widgetWH=({int(self.width())},{int(self.height())})")
                #     except Exception:
                #         pass
                # DEBUG: clear to bright color for first frames to ensure we see something even if shader fails
                if self.frame_counter < 3:
                    glClearColor(1.0, 0.0, 1.0, 1.0)
                    glClear(GL_COLOR_BUFFER_BIT)
                for i, tname in enumerate(bind_textures or []):
                    self.ogl_system.bind_texture_unit(tname, i)
                glUseProgram(program)
                self.ogl_system.draw_fullscreen_quad()
                glUseProgram(0)
            else:
                self.ogl_system.render_fs_to_fbo(program, out_fbo, bind_textures or [])
        # Restore viewport and depth/scissor tests
        if vp is not None and len(vp) == 4:
            glViewport(int(vp[0]), int(vp[1]), int(vp[2]), int(vp[3]))
        if depth_was_enabled: glEnable(GL_DEPTH_TEST)
        if scissor_was_enabled: glEnable(GL_SCISSOR_TEST)

    def resizeGL(self, w, h):
        glViewport(0, 0, w, h)
        self.projection_matrix.setToIdentity()
        if h > 0:
            aspect = w / h
            self.projection_matrix.perspective(45.0, aspect, 0.1, 1000.0)
        # Re-apply static uniforms affected by size change
        try:
            self.makeCurrent()
            self.apply_viewproj_for_geometry()
            self.apply_fs_static_uniforms()
        finally:
            self.doneCurrent()

    def apply_static_uniforms_all(self):
        """Set all non-per-frame uniforms once (projection/view, color, FS samplers, iResolution, driver)."""
        self.apply_viewproj_for_geometry()
        self.apply_fs_static_uniforms()

    def apply_viewproj_for_geometry(self):
        if not self.ogl_system: return
        # Build view matrix from current camera
        vm = QMatrix4x4()
        vm.setToIdentity()
        vm.translate(-self.camera_pos)
        vm = self.cam_rot * vm
        # Find unique geometry shader programs used by render pipeline
        shader_names = [p[0] for p in self.render_pipeline_info if len(p) > 0]
        seen = set()
        for sname in shader_names:
            prog = self.ogl_system.get_shader_program(sname)
            if not prog or prog in seen: continue
            seen.add(prog)
            glUseProgram(prog)
            # projection/view
            loc = self._get_uniform_location(prog, "projection")
            if loc != -1: glUniformMatrix4fv(loc, 1, GL_FALSE, self.projection_matrix.data())
            loc = self._get_uniform_location(prog, "view")
            if loc != -1: glUniformMatrix4fv(loc, 1, GL_FALSE, vm.data())
            # static color (if used)
            loc = self._get_uniform_location(prog, "color")
            if loc != -1: glUniform4f(loc, 0.8, 0.8, 1.0, 0.5)
            glUseProgram(0)

    def apply_fs_static_uniforms(self):
        print("apply_fs_static_uniforms()")
        if not self.ogl_system: return
        # Determine for each FS program: max sampler units used, representative resolution, and declared uniforms
        # prog_info: prog -> {"max_units":N, "res":(w,h), "uniforms": set(names)}
        prog_info = {}
        for (shader_name, out_fbo, bind_textures) in self.fs_pipeline:
            prog = self.ogl_system.get_shader_program(shader_name)
            if not prog: continue
            info = prog_info.get(prog)
            if info is None:
                info = {"max_units": 0, "res": None, "uniforms": set()}
                prog_info[prog] = info
            # sampler units
            mu = max(info["max_units"], len(bind_textures or []))
            info["max_units"] = mu
            # pick resolution for this program; warn on mismatch
            pw, ph = int(self.width()), int(self.height())
            if out_fbo and out_fbo not in ("default", "", None):
                tname = self.fs_fbos_cfg.get(out_fbo)
                if tname and tname in self.fs_textures_cfg:
                    w, h = self.fs_textures_cfg[tname]
                    pw, ph = int(w), int(h)
            if info["res"] is None:
                info["res"] = (pw, ph)
            elif info["res"] != (pw, ph):
                try:
                    print(f"apply_fs_static_uniforms: WARNING program {prog} used with multiple resolutions {info['res']} vs {(pw,ph)}; using last.")
                except Exception:
                    pass
                info["res"] = (pw, ph)
            # aggregate declared uniforms from shader config
            try:
                cfg_entry = self.ogl_system.shader_configs.get(shader_name)
                if cfg_entry and len(cfg_entry) >= 3 and isinstance(cfg_entry[2], (list, tuple)):
                    info["uniforms"].update(cfg_entry[2])
            except Exception:
                pass
        # Apply uniforms per program
        params = (self.browser.current_config.get("parameters", {}) if (self.browser and self.browser.current_config) else {})
        for prog, info in prog_info.items():
            glUseProgram(prog)
            # Sampler bindings: iChannel{i}/tex{i}/s{i} -> unit i
            for i in range(info["max_units"]):
                for uname in (f"iChannel{i}", f"tex{i}", f"s{i}"):
                    loc = self._get_uniform_location(prog, uname)
                    if loc != -1: glUniform1i(loc, i)
            # Resolution (static here; iFrame remains dynamic in _execute_fs_pipeline)
            w, h = info["res"] if info["res"] is not None else (int(self.width()), int(self.height()))
            loc = self._get_uniform_location(prog, "iResolution")
            if loc != -1: glUniform3f(loc, float(w), float(h), 1.0)
            # Other declared uniforms resolved from parameters by name
            for uname in sorted(info["uniforms"]):
                if uname in ("iFrame", "iResolution"):  # handled elsewhere
                    continue
                # skip explicit sampler aliases if present in declared list; already set above
                if uname.startswith("iChannel") or uname.startswith("tex") or uname.startswith("s") and uname[1:].isdigit():
                    continue
                loc = self._get_uniform_location(prog, uname)
                if loc == -1:
                    continue
                if uname not in params:
                    # Uniform declared for shader but no value provided in parameters; skip silently
                    continue
                val, typ, _step = params.get(uname, (None, None, None))
                try:
                    if typ == "int":
                        glUniform1i(loc, int(val))
                    elif typ == "float":
                        glUniform1f(loc, float(val))
                    elif typ == "vec2" and isinstance(val, (list, tuple)) and len(val) >= 2:
                        glUniform2f(loc, float(val[0]), float(val[1]))
                    elif typ == "vec3" and isinstance(val, (list, tuple)) and len(val) >= 3:
                        glUniform3f(loc, float(val[0]), float(val[1]), float(val[2]))
                    elif typ == "vec4" and isinstance(val, (list, tuple)) and len(val) >= 4:
                        glUniform4f(loc, float(val[0]), float(val[1]), float(val[2]), float(val[3]))
                    else:
                        # Unsupported or mismatched type; ignore
                        pass
                except Exception:
                    # Be tolerant to bad user values; fail quietly for uniforms
                    pass
            glUseProgram(0)

    def _get_uniform_location(self, program, name):
        """Return cached uniform location for (program,name), querying and caching on a miss."""
        d = self._uniform_cache.get(program)
        if d is None:
            d = {}
            self._uniform_cache[program] = d
        loc = d.get(name)
        if loc is None:
            loc = glGetUniformLocation(program, name)
            d[name] = loc
        return loc

    def wheelEvent(self, event):
        delta = event.angleDelta().y() / 120
        self.camera_pos.setZ(self.camera_pos.z() - delta)
        # Apply view/projection immediately (event-driven)
        try:
            if self.isValid():
                self.makeCurrent()
                self.apply_viewproj_for_geometry()
        finally:
            if self.isValid():
                self.doneCurrent()
        self.update()

    def mousePressEvent(self, event):
        if event.buttons() == Qt.LeftButton:
            self.last_mouse_pos = event.pos()

    def mouseMoveEvent(self, event):
        if event.buttons() == Qt.LeftButton and self.last_mouse_pos:
            dx = event.x() - self.last_mouse_pos.x()
            dy = event.y() - self.last_mouse_pos.y()
            
            rotation_speed = 0.5
            rot_y = QMatrix4x4()
            rot_y.rotate(dx * rotation_speed, 0, 1, 0)
            
            rot_x = QMatrix4x4()
            rot_x.rotate(dy * rotation_speed, 1, 0, 0)
            
            self.cam_rot = rot_y * rot_x * self.cam_rot
            self.last_mouse_pos = event.pos()
            # Apply view/projection immediately (event-driven)
            try:
                if self.isValid():
                    self.makeCurrent()
                    self.apply_viewproj_for_geometry()
            finally:
                if self.isValid():
                    self.doneCurrent()
            self.update()

    def mouseReleaseEvent(self, event):
        self.last_mouse_pos = None