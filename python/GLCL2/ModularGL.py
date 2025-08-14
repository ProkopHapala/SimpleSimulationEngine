"""
Modular OpenGL System for FireCore

This module provides a highly modular and reusable OpenGL system
that achieves maximum flexibility with minimal code duplication.
"""

import numpy as np
from OpenGL.GL import *
from .OGLsystem import GLobject, InstancedData, Mesh, OGLSystem, create_sphere_mesh

class BufferObject:
    """Generic buffer object that can handle any type of data"""
    
    def __init__(self, name, src=None, usage=GL_DYNAMIC_DRAW, elem_size=None):
        """
        Polymorphic initialization - can create buffer immediately or defer
        
        Args:
            name: Buffer name
            src: Optional data to initialize with (None for deferred)
            usage: GL usage constant
            elem_size: Optional element size for empty allocation
        """
        self.name = name
        self.usage = usage
        self.vbo = None
        self.data = None
        self.size = 0
        self.elem_size = elem_size
        
        # If src provided, create immediately (concise one-liner)
        if src is not None:
            self.create(src)
        elif elem_size is not None:
            # Pre-allocate empty buffer
            self.allocate_empty(elem_size)
            
    def create(self, data):
        """Create buffer with given data"""
        self.data = data
        self.size = data.nbytes
        
        if self.vbo is None:
            self.vbo = glGenBuffers(1)
            
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, self.size, data, self.usage)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        
    def allocate_empty(self, elem_size, count=1):
        """Allocate empty buffer of given size"""
        self.size = elem_size * count
        
        if self.vbo is None:
            self.vbo = glGenBuffers(1)
            
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glBufferData(GL_ARRAY_BUFFER, self.size, None, self.usage)
        glBindBuffer(GL_ARRAY_BUFFER, 0)
        
    def update(self, data):
        """Update buffer data"""
        if self.vbo is None:
            self.create(data)
            return
            
        if data.nbytes != self.size:
            self.create(data)  # Recreate if size changed
        else:
            glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
            glBufferSubData(GL_ARRAY_BUFFER, 0, data.nbytes, data)
            glBindBuffer(GL_ARRAY_BUFFER, 0)
            
    def bind(self, location, components=4, stride=0, offset=0):
        """Bind buffer to shader attribute location"""
        glBindBuffer(GL_ARRAY_BUFFER, self.vbo)
        glVertexAttribPointer(location, components, GL_FLOAT, GL_FALSE, stride, ctypes.c_void_p(offset))
        glEnableVertexAttribArray(location)
        
    def cleanup(self):
        """Clean up buffer"""
        if self.vbo is not None:
            glDeleteBuffers(1, [self.vbo])
            self.vbo = None


class RenderObject:
    """Generic render object that encapsulates VAO/VBO management"""
    
    def __init__(self, name, buffers=None, shader_program=None, draw_mode=GL_TRIANGLES, element_count=0):
        """
        Polymorphic initialization - multiple ways to create
        
        Args:
            name: Object name
            buffers: Dict of {name: data} to create immediately
            shader_program: Optional shader program
            draw_mode: GL draw mode
            element_count: Number of elements to draw
        """
        self.name = name
        self.vao = None
        self.buffers = {}  # name -> BufferObject mapping
        self.shader_program = shader_program
        self.draw_mode = draw_mode
        self.element_count = element_count
        
        # Create VAO immediately
        self.create_vao()
        
        # If buffers provided, create them immediately
        if buffers is not None:
            if isinstance(buffers, dict):
                # Dict of {name: data}
                for name, data in buffers.items():
                    self.add_buffer(name, data)
            elif isinstance(buffers, (list, tuple, np.ndarray)):
                # Single buffer data - use name as key
                self.add_buffer(name, buffers)
                
    def create_vao(self):
        """Create VAO for this render object"""
        if self.vao is None:
            self.vao = glGenVertexArrays(1)
        return self.vao
        
    def add_buffer(self, name, data=None, location=None, components=4, src=None, 
                   elem_size=None, count=None):
        """
        Polymorphic buffer addition - multiple ways to call
        
        Args:
            name: Buffer name (required)
            data: Buffer data (optional)
            location: Shader attribute location
            components: Number of components per element
            src: Alternative data source
            elem_size: Element size for empty allocation
            count: Number of elements for empty allocation
        """
        # Handle polymorphic data source
        buffer_data = data if data is not None else src
        
        if buffer_data is None and (elem_size is None or count is None):
            # Create empty buffer
            buffer = BufferObject(name)
        elif buffer_data is None and elem_size is not None and count is not None:
            # Pre-allocate empty buffer with specific size
            buffer = BufferObject(name, elem_size=elem_size * count)
        else:
            # Create with data
            buffer = BufferObject(name, src=buffer_data)
            
        self.buffers[name] = buffer
        
        if location is not None:
            self.bind_buffer(name, location, components)
            
    def add_buffers(self, buffers_dict, locations=None):
        """
        Add multiple buffers at once - polymorphic
        
        Args:
            buffers_dict: Dict of {name: data}
            locations: Dict of {name: location} or list of locations
        """
        for name, data in buffers_dict.items():
            location = None
            if isinstance(locations, dict):
                location = locations.get(name)
            elif isinstance(locations, (list, tuple)) and len(locations) > 0:
                location = locations[0] if len(locations) == 1 else locations[len(self.buffers)]
                
            self.add_buffer(name, data, location=location)
            
    def bind_buffer(self, name, location, components=4):
        """Bind buffer to shader attribute location"""
        if name in self.buffers:
            glBindVertexArray(self.vao)
            self.buffers[name].bind(location, components)
            glBindVertexArray(0)
            
    def update_buffer(self, name, data):
        """Update buffer data"""
        if name in self.buffers:
            self.buffers[name].update(data)
        else:
            # Create new buffer if doesn't exist
            self.add_buffer(name, data)
            
    def draw(self, count=None):
        """Draw this render object"""
        if count is None:
            count = self.element_count
            
        glBindVertexArray(self.vao)
        if self.element_count > 0:
            glDrawArrays(self.draw_mode, 0, count)
        glBindVertexArray(0)
        
    def cleanup(self):
        """Clean up all resources"""
        for buffer in self.buffers.values():
            buffer.cleanup()
        if self.vao is not None:
            glDeleteVertexArrays(1, [self.vao])
            self.vao = None


class InstancedRenderObject:
    """Instanced rendering object for efficient batch rendering"""
    
    def __init__(self, name, base_mesh=None, instance_attribs=None, instance_data=None):
        """
        Polymorphic initialization - can setup immediately or defer
        
        Args:
            name: Object name
            base_mesh: Optional base mesh
            instance_attribs: Dict of {attrib_name: components}
            instance_data: Optional instance data to load immediately
        """
        self.name = name
        self.base_mesh = base_mesh
        self.instanced_data = None
        self.instance_count = 0
        
        # Setup immediately if all components provided
        if base_mesh is not None and instance_attribs is not None:
            self.setup(base_mesh, instance_attribs)
            
        if instance_data is not None:
            self.update_instances(instance_data)
            
    def setup(self, base_mesh=None, instance_attribs=None):
        """Setup instanced rendering - polymorphic"""
        if base_mesh is not None:
            self.base_mesh = base_mesh
            
        if self.base_mesh is None:
            return
            
        self.instanced_data = InstancedData(base_attrib_location=2)
        self.instanced_data.associate_mesh(self.base_mesh)
        
        # Handle multiple input formats
        if instance_attribs is not None:
            attribs = []
            if isinstance(instance_attribs, dict):
                # Dict format: {name: components}
                for i, (name, components) in enumerate(instance_attribs.items()):
                    attribs.append((name, i, components))
            elif isinstance(instance_attribs, (list, tuple)):
                # List format: [(name, components), ...]
                for i, (name, components) in enumerate(instance_attribs):
                    attribs.append((name, i, components))
            
            self.instanced_data.setup_instance_vbos(attribs)
        
    def update_instances(self, instance_data=None, **kwargs):
        """
        Polymorphic instance data update
        
        Args:
            instance_data: Dict of instance data
            **kwargs: Alternative way to provide data
        """
        # Handle polymorphic data sources
        data = instance_data if instance_data is not None else kwargs
        
        if data and self.instanced_data is not None:
            self.instanced_data.update(data)
            if isinstance(data, dict):
                self.instance_count = len(next(iter(data.values())))
            elif isinstance(data, (list, tuple)) and len(data) > 0:
                self.instance_count = len(data[0])
        
    def draw(self):
        """Draw instanced objects"""
        if self.instanced_data:
            self.instanced_data.draw()
            
    def cleanup(self):
        """Clean up resources"""
        if self.instanced_data:
            self.instanced_data.cleanup()


class ShaderProgram:
    """Encapsulates shader program management"""
    
    def __init__(self, name, vertex_src, fragment_src):
        self.name = name
        self.vertex_src = vertex_src
        self.fragment_src = fragment_src
        self.program_id = None
        
    def compile(self):
        """Compile shader program"""
        from .OGLsystem import compile_shader_program
        self.program_id = compile_shader_program(self.vertex_src, self.fragment_src)
        return self.program_id is not None
        
    def use(self):
        """Use this shader program"""
        if self.program_id is not None:
            glUseProgram(self.program_id)
            
    def get_uniform_location(self, name):
        """Get uniform location"""
        if self.program_id is not None:
            return glGetUniformLocation(self.program_id, name)
        return -1
        
    def set_uniform_matrix4fv(self, name, matrix):
        """Set 4x4 matrix uniform"""
        location = self.get_uniform_location(name)
        if location != -1:
            glUniformMatrix4fv(location, 1, GL_FALSE, matrix)
            
    def set_uniform3fv(self, name, vector):
        """Set 3-component vector uniform"""
        location = self.get_uniform_location(name)
        if location != -1:
            glUniform3fv(location, 1, vector)
            
    def set_uniform1f(self, name, value):
        """Set float uniform"""
        location = self.get_uniform_location(name)
        if location != -1:
            glUniform1f(location, value)


class RenderSystem:
    """High-level render system that manages all render objects"""
    
    def __init__(self):
        self.render_objects = {}
        self.instanced_objects = {}
        self.shaders = {}
        self.current_shader = None
        
    def add_render_object(self, name, render_obj):
        """Add render object to system"""
        self.render_objects[name] = render_obj
        
    def add_instanced_object(self, name, instanced_obj):
        """Add instanced render object"""
        self.instanced_objects[name] = instanced_obj
        
    def add_shader(self, name, shader):
        """Add shader program"""
        self.shaders[name] = shader
        
    def use_shader(self, name):
        """Use specific shader"""
        if name in self.shaders:
            self.shaders[name].use()
            self.current_shader = self.shaders[name]
            
    def render_all(self, view_matrix=None, proj_matrix=None):
        """Render all objects with optional matrices"""
        for name, obj in self.render_objects.items():
            if self.current_shader is not None:
                if view_matrix is not None:
                    self.current_shader.set_uniform_matrix4fv("view", view_matrix)
                if proj_matrix is not None:
                    self.current_shader.set_uniform_matrix4fv("proj", proj_matrix)
            obj.draw()
            
        for name, obj in self.instanced_objects.items():
            if self.current_shader is not None:
                if view_matrix is not None:
                    self.current_shader.set_uniform_matrix4fv("view", view_matrix)
                if proj_matrix is not None:
                    self.current_shader.set_uniform_matrix4fv("proj", proj_matrix)
            obj.draw()
            
    def cleanup(self):
        """Clean up all resources"""
        for obj in self.render_objects.values():
            obj.cleanup()
        for obj in self.instanced_objects.values():
            obj.cleanup()


# Convenience factory functions for common use cases
class RenderFactory:
    """Factory for creating common render objects"""
    
    @staticmethod
    def create_particle_system(name, positions=None, colors=None, size=1.0, **kwargs):
        """
        Polymorphic particle system creation
        
        Args:
            name: System name
            positions: Particle positions (optional)
            colors: Particle colors (optional)
            size: Particle size
            **kwargs: Additional parameters
        """
        buffers = {}
        element_count = 0
        
        # Build buffers dict polymorphically
        if positions is not None:
            buffers["positions"] = positions
            element_count = len(positions)
            
        if colors is not None:
            buffers["colors"] = colors
            
        # Create render object with polymorphic initialization
        render_obj = RenderObject(name, buffers=buffers, element_count=element_count, **kwargs)
        return render_obj
        
    @staticmethod
    def create_mesh_renderer(name, vertices, normals=None, indices=None):
        """Create mesh renderer"""
        mesh = Mesh(vertices, normals, indices)
        mesh.setup_buffers()
        
        render_obj = RenderObject(name)
        render_obj.create_vao()
        render_obj.element_count = len(vertices) // 3
        return render_obj
        
    @staticmethod
    def create_sphere_instancer(name, instance_positions=None, radius=1.0, **kwargs):
        """
        Polymorphic sphere instancer creation
        
        Args:
            name: System name
            instance_positions: Instance positions (optional)
            radius: Sphere radius
            **kwargs: Additional parameters
        """
        # Create sphere mesh
        sphere_vertices, sphere_normals, sphere_indices = create_sphere_mesh(radius)
        sphere_mesh = Mesh(sphere_vertices, sphere_normals, sphere_indices)
        sphere_mesh.setup_buffers()
        
        # Polymorphic instance attributes
        instance_attribs = {"positions": 4}
        
        # Create instanced object with polymorphic initialization
        instanced_obj = InstancedRenderObject(
            name, 
            base_mesh=sphere_mesh,
            instance_attribs=instance_attribs,
            instance_data={"positions": instance_positions} if instance_positions is not None else None,
            **kwargs
        )
        return instanced_obj
