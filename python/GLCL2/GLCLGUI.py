import sys
import numpy as np
import os

from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QMatrix4x4, QVector3D
from OpenGL.GL import *

from .OGLsystem import GLobject

class GLCLWidget(QOpenGLWidget):
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
        finally:
            self.doneCurrent()
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

        for buffer_name in vertex_buffers_to_bake:
            data = self.buffer_data.get(buffer_name)
            if data is None:
                print(f"Warning: No initial data for buffer '{buffer_name}'. Cannot bake GLobject.")
                continue

            nelements = len(data)
            components = data.shape[1] if data.ndim > 1 else 1
            
            gl_obj = GLobject(nelements=nelements, mode=GL_POINTS)
            gl_obj.alloc_vao_vbo_ebo([components])
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
                shader_name, _, vertex_buffer_name, _ = pass_info
                
                shader_program = self.ogl_system.get_shader_program(shader_name)
                gl_obj = self.gl_objects.get(vertex_buffer_name)

                if not shader_program or not gl_obj:
                    continue

                glUseProgram(shader_program)
                
                # Set generic uniforms
                self.update_matrices(shader_program)
                color_loc = glGetUniformLocation(shader_program, "color")
                if color_loc != -1: glUniform4f(color_loc, 0.8, 0.8, 1.0, 0.5)

                # Draw the object
                gl_obj.draw_arrays()

            glUseProgram(0)
        except Exception as e:
            if self.browser is not None:
                self.browser.on_exception(e)
            else:
                import traceback, os
                traceback.print_exc()
                os._exit(1)

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
            if self.frame_counter < 5:
                try:
                    print(f"[FS DEBUG] frame={self.frame_counter} pass='{shader_name}' out='{out_fbo}' inputs={bind_textures}")
                except Exception:
                    pass
            program = self.ogl_system.get_shader_program(shader_name)
            if not program:
                raise KeyError(f"FS shader '{shader_name}' not compiled or registered")
            # Set common sampler uniforms to match bound units
            glUseProgram(program)
            max_units = 8
            for i in range(min(len(bind_textures or []), max_units)):
                for uname in (f"iChannel{i}", f"tex{i}", f"s{i}"):
                    loc = glGetUniformLocation(program, uname)
                    if loc != -1: glUniform1i(loc, i)
            # Set iFrame
            loc_frame = glGetUniformLocation(program, "iFrame")
            if loc_frame != -1:
                glUniform1i(loc_frame, int(self.frame_counter))
            # Determine resolution for this pass (target FBO if present, else widget size)
            pass_w = self.width()
            pass_h = self.height()
            if out_fbo and out_fbo not in ("default", "", None):
                if out_fbo not in self.fs_fbos_cfg:
                    raise KeyError(f"FS FBO '{out_fbo}' not configured")
                tex_name = self.fs_fbos_cfg.get(out_fbo)
                if tex_name and tex_name in self.fs_textures_cfg:
                    w, h = self.fs_textures_cfg[tex_name]
                    pass_w, pass_h = int(w), int(h)
            loc_res = glGetUniformLocation(program, "iResolution")
            if loc_res != -1:
                glUniform3f(loc_res, float(pass_w), float(pass_h), 1.0)
            # Optional custom uniform: driver (vec4)
            try:
                cfg = self.browser.current_config if self.browser else None
                if cfg is not None and "driver" in cfg.get("parameters", {}):
                    drv_val, drv_type, _ = cfg["parameters"]["driver"]
                    if drv_type == "vec4" and isinstance(drv_val, (list, tuple)) and len(drv_val) == 4:
                        loc_drv = glGetUniformLocation(program, "driver")
                        if loc_drv != -1:
                            glUniform4f(loc_drv, float(drv_val[0]), float(drv_val[1]), float(drv_val[2]), float(drv_val[3]))
            except Exception:
                # DEBUG: fail loud in browser if desired; here we keep rendering even if driver missing
                pass
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
                if self.frame_counter < 5:
                    try:
                        fb = glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING)
                        vp2 = glGetIntegerv(GL_VIEWPORT)
                        print(f"[FS DEBUG] default pass fbo={int(fb)} viewport={tuple(int(x) for x in vp2)} widgetWH=({int(self.width())},{int(self.height())})")
                    except Exception:
                        pass
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

    def update_matrices(self, shader_program):
        self.view_matrix.setToIdentity()
        self.view_matrix.translate(-self.camera_pos)
        self.view_matrix = self.cam_rot * self.view_matrix

        proj_loc = glGetUniformLocation(shader_program, "projection")
        view_loc = glGetUniformLocation(shader_program, "view")
        
        if proj_loc != -1: glUniformMatrix4fv(proj_loc, 1, GL_FALSE, self.projection_matrix.data())
        if view_loc != -1: glUniformMatrix4fv(view_loc, 1, GL_FALSE, self.view_matrix.data())

    def wheelEvent(self, event):
        delta = event.angleDelta().y() / 120
        self.camera_pos.setZ(self.camera_pos.z() - delta)
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
            self.update()

    def mouseReleaseEvent(self, event):
        self.last_mouse_pos = None