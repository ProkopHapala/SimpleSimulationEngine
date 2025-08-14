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

    def initializeGL(self):
        print("GLCLWidget::initializeGL() - OpenGL Context is now valid.")
        glClearColor(0.1, 0.1, 0.15, 1.0)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_PROGRAM_POINT_SIZE)
        
        # Now it's safe to compile shaders and create GL resources
        self.ogl_system.compile_all_shaders()
        self.bake_render_objects()

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
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        
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