"""simulation.py – Core of pySymGLSL.

This is deliberately *small* and *stateless*, merely orchestrating moderngl
objects laid out in simple dictionaries.  It supports:

*  Loading GLSL programs (`load_program`).  Vertex shader is optional – a
   trivial full-screen quad vertex shader is injected if none is supplied.
*  Creating *ping-pong* texture + FBO pairs with a single call
   (`create_texture`).
*  Baking a *render pass* definition (program, output, inputs, uniforms) into a
   lambda that can be executed with near-zero overhead (`bake_pass`).
*  Grouping passes returned by `bake_pass` into a *render graph* and executing
   them in sequence (`run_graph`).

"""
from __future__ import annotations
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple
import moderngl
import numpy as np

# -----------------------------------------------------------------------------
# Helpers – trivial vertex shader & fullscreen quad VAO
# -----------------------------------------------------------------------------

_DEFAULT_VS = """#version 330 core
in  vec2 in_position;
out vec2 v_texcoord;
void main(){
    v_texcoord = (in_position + 1.0) * 0.5;   // map [-1,1] square to [0,1]
    gl_Position = vec4(in_position, 0.0, 1.0);
}"""

_QUAD_VERTICES = np.array([[-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1.0], [-1.0, 1.0]], dtype='f4')

#_QUAD_VERTICES = np.array( [ [-1.0,-1.0,1.0,], [-1.0,1.0,1.0], [1.0,-1.0,1.0], [1.0,1.0,1.0] ], dtype="f4" )

import re

class GLSL_Simulation:
    """Main entry point.

    Parameters
    ----------
    ctx
        Existing *moderngl* context or *None*.  If *None*, a standalone context
        is created (headless).  Provide your own context if you need a window.
    sim_size
        Width × height of the simulation grid / textures.
    dtype
        Numpy dtype string for created textures.  Default 'f4' (32-bit float).
    """

    def __init__(self, sim_size: Tuple[int, int], ctx: moderngl.Context | None = None, *, dtype: str = "f4"):
        self.ctx = ctx or moderngl.create_standalone_context(require=330)
        self.sim_size = sim_size
        self.dtype = dtype

        # Resource registries --------------------------------------------------
        self.programs: Dict[str, moderngl.Program] = {}
        self.textures: Dict[str, moderngl.Texture] = {}
        self.samplers: Dict[str, moderngl.Sampler] = {}
        self.framebuffers: Dict[str, moderngl.Framebuffer] = {}

        self._current_pipeline = None
        self._current_base_dir = None
        self._current_graph_definition = None

        # One VAO that renders a fullscreen quad ------------------------------
        vbo = self.ctx.buffer(_QUAD_VERTICES.tobytes())
        vao_content = [(vbo, "2f", "in_position")]
        self.quad_vao = self.ctx.vertex_array(self.ctx.program(vertex_shader=_DEFAULT_VS, fragment_shader="void main(){}"), vao_content)
        # Replace the dummy program later when first render pass is executed.
        self._vao_prog_placeholder = self.quad_vao.program  # keep reference to dispose later
        # Store quad VBO and content for per-pass VAO creation
        self._quad_vbo = vbo
        self._quad_content = vao_content
        self.iFrame = 0

    def setup_texture(self, size: Tuple[int, int], channels=4, dtype="f4", repeat=(False, False), filter=moderngl.LINEAR, anisotropy=1.0):
        tex = self.ctx.texture(size, components=channels, dtype=dtype)
        sampler = self.ctx.sampler(
            repeat_x=repeat[0],
            repeat_y=repeat[1],
            filter=(filter, filter),
            anisotropy=anisotropy
        )
        return tex, sampler

    def build_pipeline(self, pipeline: list, base_dir):
        self._current_pipeline = pipeline
        self._current_base_dir = base_dir
        print("build_pipeline(): ", pipeline," base_dir: ", base_dir)
        # Auto-load any shader programs mentioned in Pipeline but not declared
        for prog_name, _, _, _ in pipeline:
            if prog_name not in self.programs:
                # treat prog_name as relative path
                self.load_program(prog_name, fragment_path=base_dir / prog_name)

        # 2. Determine buffer names and allocate single-framebuffer textures
        tex_names = set()
        for prog_name, out_name, inputs_map, dyn_uniforms in pipeline:
            tex_names.add(out_name)
            tex_names.update(inputs_map.values())

        # Allocate textures + FBOs
        for t in tex_names:
            if t not in self.textures:
                size_x, size_y = self.sim_size
                tex, sampler = self.setup_texture((size_x, size_y), dtype=self.dtype)
                fbo = self.ctx.framebuffer(color_attachments=[tex])
                self.textures[t] = tex
                self.samplers[t] = sampler
                self.framebuffers[t] = fbo
                
        # Initialize all textures to avoid undefined behavior in feedback loops
        self.initialize_textures(v=(0.0,0.0,0.0,0.0))

        # 3. Build pass descriptors ------------------------------------------------
        passes = []
        for prog_name, out_name, inputs_map, dyn_uniforms in pipeline:
            # Collect texture names in channel order: iChannel0, iChannel1, ...
            input_names = []
            if isinstance(inputs_map, dict):
                for ch in range(8):
                    key = f"iChannel{ch}"
                    if key in inputs_map:
                        input_names.append(inputs_map[key])
            # Ensure iResolution is available as dynamic uniform (convenience)
            if 'iResolution' not in dyn_uniforms: dyn_uniforms = dyn_uniforms + ['iResolution']
            if 'iFrame'      not in dyn_uniforms: dyn_uniforms = dyn_uniforms + ['iFrame']
            passes.append((prog_name, out_name, input_names, dyn_uniforms))

        baked_pipeline = self.bake_graph(passes)
        return baked_pipeline, tex_names


    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------
    def load_program(self, name: str, *, vertex_path: str | Path | None = None, fragment_path: str | Path) -> None:
        """Compile and register a GLSL program.

        If *vertex_path* is *None*, `_DEFAULT_VS` is used.
        Raises *RuntimeError* on compilation failure with error log attached.
        """
        vs_source = Path(vertex_path).read_text() if vertex_path else _DEFAULT_VS
        # load fragment source and strip single-line comments to avoid commented includes
        raw_fs = Path(fragment_path).read_text()
        fs_source = re.sub(r'//.*', '', raw_fs)
        try:
            program = self.ctx.program(vertex_shader=vs_source, fragment_shader=fs_source)
        except moderngl.Error as exc:
            # Provide nicer error with file path & log
            raise RuntimeError(f"Shader compilation failed for '{fragment_path}':\n{exc}") from exc
        self.programs[name] = program

    # ------------------------------------------------------------------
    # Render-pass baking & execution
    # ------------------------------------------------------------------
    def bake_pass(
        self,
        program_name: str,
        output_name: str,
        input_names: Sequence[str] | None = None,
        dynamic_uniforms: Sequence[str] | None = None,
    ) -> Callable[[dict], None]:
        """Return a callable that executes the given render pass.

        The callable has the signature ``func(dynamic_values: dict)`` where
        *dynamic_values* maps uniform names (given in *dynamic_uniforms*) to
        actual values.
        """
        input_names = input_names or []
        dynamic_uniforms = dynamic_uniforms or []

        program = self.programs[program_name]

        # Resolve uniform locations once to avoid dictionary lookups per frame
        uniform_setters: List[Callable[[dict], None]] = []
        for uni in dynamic_uniforms:
            if uni not in program:
                raise KeyError(f"Uniform '{uni}' not found in program '{program_name}'")
            loc = program[uni]
            # Pre-bind a setter that handles scalar and vector uniforms
            # Determine expected component count from loc.fmt if available
            fmt = getattr(loc, 'fmt', None)
            print(f"DEBUG: Uniform '{uni}', fmt: {fmt}")
            try:
                count = int(fmt[:-1]) if fmt and isinstance(fmt, str) and fmt[:-1].isdigit() else 1
            except:
                count = 1
            def setter(dv, loc=loc, name=uni, count=count):
                val = dv[name]
                if not isinstance(val, (tuple, list)):
                    # for single-component uniforms, pass scalar; otherwise wrap into tuple
                    if count > 1:
                        val = (val,) * count
                elif len(val) != count:
                    raise ValueError(f"Uniform '{name}' expects {count} components, got {len(val)}")
                loc.value = val
            uniform_setters.append(setter)

        # Bind sampler uniforms for texture units 0..N. Prefer iChannelN naming (Shadertoy style);
        # fall back to u_texture_N for legacy shaders.
        for i in range(len(input_names)):
            if f"iChannel{i}" in program:
                program[f"iChannel{i}"].value = i
            elif f"u_texture_{i}" in program:
                program[f"u_texture_{i}"].value = i

        # Set default resolution uniform if present
        #if "iResolution" in program:
        program["iResolution"].value = (self.sim_size[0], self.sim_size[1], 1.0)
        program["iFrame"]     .value = self.iFrame

        # Create a VAO for this render pass
        vao = self.ctx.vertex_array(program, self._quad_content)

        def _execute(dynamic_values: dict) -> None:
            # Dynamically retrieve output_fbo and input_textures
            output_fbo = self.framebuffers[output_name]
            # input_textures must be retrieved dynamically inside _execute
            input_textures_dynamic = [self.textures[n] for n in input_names]
            input_samplers_dynamic = [self.samplers[n] for n in input_names]

            # Use program for rendering
            # Program is bound via VAO

            # Dynamic uniform updates
            for setter in uniform_setters:
                setter(dynamic_values)

            # Bind input textures to units 0..N
            for i, (tex, sampler) in enumerate(zip(input_textures_dynamic, input_samplers_dynamic)):
                tex.use(location=i)
                sampler.use(location=i)

            # Render to output
            output_fbo.use()
            vao.render(moderngl.TRIANGLE_STRIP)

        return _execute

    # ------------------------------------------------------------------
    def bake_graph(self, graph_definition: Iterable[Tuple[str, str, List[str], List[str]]]):
        """Bake a sequence of passes into a list of callables."""
        self._current_graph_definition = graph_definition
        return [self.bake_pass(*g) for g in graph_definition]

    def initialize_textures(self, v=(0.0,0.0,0.0,0.0)): 
        """Initialize all textures with a default value to avoid undefined behavior in feedback loops."""
        # Save the current framebuffer binding
        default_framebuffer = self.ctx.fbo
        
        # Clear each texture with the specified color
        for name, fbo in self.framebuffers.items():
            fbo.use()
            self.ctx.clear(*v)
            print(f"Initialized texture '{name}' with data {v}")
        
        # Restore the original framebuffer
        default_framebuffer.use()
        
    def run_graph(self, baked_graph: Sequence[Callable[[dict], None]], dynamic_values: dict):
        """Execute all passes in *baked_graph* with shared *dynamic_values*."""
        for execute in baked_graph:
            execute(dynamic_values)
        self.iFrame += 1

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    def clear(self, color=(0.0, 0.0, 0.0, 0.0)):
        self.ctx.clear(*color)

    def release(self):
        """Free GL resources (optional)."""
        for prog in self.programs.values():
            prog.release()
        for tex in self.textures.values():
            tex.release()
        for sampler in self.samplers.values():
            sampler.release()
        for fbo in self.framebuffers.values():
            fbo.release()
        self.quad_vao.release()
        self._vao_prog_placeholder.release()
