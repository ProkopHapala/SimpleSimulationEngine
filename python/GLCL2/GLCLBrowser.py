import sys
import os
import importlib.util
import numpy as np
import pyopencl as cl
import copy

from PyQt5.QtWidgets import QApplication, QHBoxLayout, QWidget, QFileDialog, QVBoxLayout, QGroupBox
from PyQt5 import QtWidgets
from PyQt5.QtCore import QTimer, Qt

import OpenGL.GL as GL
import argparse

from ..pySymGLSL.BaseGUI import BaseGUI
from .GLCLGUI import GLCLWidget
from .OCLsystem import OCLSystem
from .OGLsystem import OGLSystem

class BakedKernelCall:
    """A fully pre-baked kernel execution object."""
    def __init__(self, ocl_system, kernel_obj, global_size, local_size, args_tuple):
        self.ocl_system = ocl_system
        self.kernel_obj = kernel_obj
        self.kernel_name = kernel_obj.function_name
        self.global_size = global_size
        self.local_size = local_size
        self.args_tuple = args_tuple
    
    def execute(self):
        """Execute the kernel with pre-baked arguments."""
        try:
            self.ocl_system.execute_kernel(
                self.kernel_obj, self.global_size, self.local_size, *self.args_tuple
            )
        except Exception as e:
            print(f"Error executing kernel {self.kernel_name}: {e}")
            import traceback
            traceback.print_exc()
            raise

class GLCLBrowser(BaseGUI):
    """Main browser class for loading and executing scientific simulation scripts."""
    
    def __init__(self, python_script_path=None, bDebugCL=False, bDebugGL=False, nDebugFrames=0, fps=60, start_paused=False, frame_limit=0, frame_delay_ms=None):
        super().__init__()
        self.setWindowTitle("GLCL Browser - Scientific Simulation Framework")
        self.setGeometry(100, 100, 1600, 900)

        self.ogl_system = OGLSystem()
        self.ocl_system = OCLSystem()
        
        self.bDebugCL = bDebugCL
        self.bDebugGL = bDebugGL
        self.nDebugFrames = nDebugFrames
        self.debug_frame_counter = 0
        self.fps = float(fps)
        self.start_paused = bool(start_paused)
        self.frame_limit = int(frame_limit)  # 0 -> unlimited
        # If frame_delay_ms is provided, it takes precedence over fps
        if frame_delay_ms is not None:
            self.frame_interval_ms = max(0, int(frame_delay_ms))
        else:
            self.frame_interval_ms = max(0, int(1000.0 / self.fps))

        self.current_script = None
        self.current_config = None
        self.baked_kernel_calls = []
        self.param_widgets = {}
        self.buffer_shapes = {}
        
        self.buffer_data = {}
        self.render_pipeline_info = []
        self.buffers_to_sync = []
        
        self.main_widget = QWidget()
        self.setCentralWidget(self.main_widget)
        
        self.glcl_widget = GLCLWidget(self, enable_opengl_debug=self.bDebugGL)
        self.glcl_widget.set_systems(self.ogl_system, self.ocl_system, self)

        self._build_ui()
        
        if python_script_path is not None:
            self.load_and_apply_script(python_script_path)
            
    def _build_ui(self):
        main_layout = QHBoxLayout()
        self.main_widget.setLayout(main_layout)

        control_panel = QWidget()
        control_layout = QVBoxLayout()
        control_panel.setLayout(control_layout)
        control_panel.setMaximumWidth(300)
        
        file_group = QGroupBox("Simulation Script")
        file_layout = QVBoxLayout()
        file_group.setLayout(file_layout)
        self.button("Load Script", self.load_simulation_script, layout=file_layout)
        self.button("Reload Current", self.reload_current_script, layout=file_layout)
        control_layout.addWidget(file_group)
        
        params_group = QGroupBox("Simulation Parameters")
        self.params_layout = QVBoxLayout()
        params_group.setLayout(self.params_layout)
        control_layout.addWidget(params_group)
        
        sim_group = QGroupBox("Simulation Control")
        sim_layout = QVBoxLayout()
        sim_group.setLayout(sim_layout)
        self.button("Start/Pause", self.toggle_simulation, layout=sim_layout)
        self.button("Reset", self.reset_simulation, layout=sim_layout)
        control_layout.addWidget(sim_group)
        control_layout.addStretch()

        main_layout.addWidget(self.glcl_widget, 1)
        main_layout.addWidget(control_panel)
        
        self.sim_timer = QTimer()
        self.sim_timer.timeout.connect(self.update_simulation)
        self.simulation_running = False

    def load_and_apply_script(self, script_path):
        try:
            self.simulation_running = False
            self.sim_timer.stop()

            spec = importlib.util.spec_from_file_location("simulation_script", script_path)
            module = importlib.util.module_from_spec(spec)
            sys.modules["simulation_script"] = module
            spec.loader.exec_module(module)
            
            config = getattr(module, 'config', None)
            init_func = getattr(module, 'init', None)

            if config is None:
                raise ValueError(f"Python script '{script_path}' must define a 'config' dictionary.")

            self.current_script = script_path
            self.apply_simulation_config(config, init_func, script_path)
            print(f"GLCLBrowser::load_and_apply_script() Successfully loaded: {script_path}")
            
        except Exception as e:
            print(f"GLCLBrowser::load_and_apply_script() Error loading script {script_path}: {e}")
            self.on_exception(e)

    def load_simulation_script(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Load Simulation Script", "", "Python Files (*.py)")
        if file_path:
            self.load_and_apply_script(file_path)

    def reload_current_script(self):
        if self.current_script:
            self.load_and_apply_script(self.current_script)

    def apply_simulation_config(self, config, init_func, script_path):
        print("--- Applying new simulation configuration ---")
        
        self.current_config = copy.deepcopy(config)
        script_dir = os.path.dirname(script_path) if script_path else "."

        self.ocl_system.clear()
        self.ogl_system.clear()
        for i in reversed(range(self.params_layout.count())):
            self.params_layout.itemAt(i).widget().setParent(None)
        self.param_widgets.clear()

        self.create_parameter_controls(self.current_config)
        
        self.setup_opencl_system(self.current_config, script_dir)
        
        self.init_simulation_data(self.current_config, init_func)
        
        self.bake_kernels(self.current_config)
        
        self.setup_opengl_system(self.current_config, script_dir)

        print("--- Configuration Summary ---")
        print(f"Buffers: {list(self.buffer_shapes.keys())}")
        print("Baked Kernels:")
        for kc in self.baked_kernel_calls:
            print(f"  - {kc.kernel_name}: global={kc.global_size}, local={kc.local_size}")
        print("-----------------------------")

        if not self.simulation_running and not self.start_paused:
            self.start_simulation()

    def init_simulation_data(self, config, init_func):
        print("Initializing simulation data and buffers...")
        buffers_config = config.get("buffers", {})
        parameters = config.get("parameters", {})
        
        init_data = init_func() if init_func else {}
        self.buffer_data = {}

        for name, (size_expr, stride, dtype) in buffers_config.items():
            size = self._resolve_expression(size_expr, parameters)
            shape = (size, stride)
            self.buffer_shapes[name] = shape
            
            if name in init_data:
                data = init_data[name].astype(np.float32)
            else:
                print(f"Warning: Buffer '{name}' not found in init() data. Defaulting to zeros.")
                data = np.zeros(shape, dtype=np.float32)

            self.ocl_system.create_buffer(name, data.nbytes)
            self.ocl_system.toGPU(name, data)
            self.buffer_data[name] = data
            print(f"  Created CL buffer '{name}' with shape {shape}")

    def setup_opencl_system(self, config, script_dir):
        print("Loading OpenCL source files...")
        for source_file in config.get("opencl_source", []):
            filepath = os.path.join(script_dir, source_file)
            program_name = os.path.splitext(os.path.basename(source_file))[0]
            self.ocl_system.load_program(program_name, filepath)

    def setup_opengl_system(self, config, script_dir):
        print("Preparing OpenGL shaders and render pipeline...")
        self.ogl_system.shader_configs = config.get("opengl_shaders", {})
        self.ogl_system.script_dir = script_dir
        self.render_pipeline_info = config.get("render_pipeline", [])
        
        self._precompute_buffer_sync_list(config)
        
        # Optional: FS pipeline + textures/FBOs
        # textures_2d: { name: (w_expr, h_expr) }
        # fbos: { fbo_name: color_texture_name }
        # fs_pipeline: [ (shader_name, out_fbo_name, [bind_texture_names]) ]
        textures_cfg_raw = config.get("textures_2d", {})
        fbos_cfg = config.get("fbos", {})
        fs_pipeline = config.get("fs_pipeline", [])

        # Resolve texture sizes (allow ints or parameter names)
        parameters = config.get("parameters", {})
        textures_cfg = {}
        for tname, size in textures_cfg_raw.items():
            try:
                w_expr, h_expr = size
                w = self._resolve_expression(w_expr, parameters)
                h = self._resolve_expression(h_expr, parameters)
                textures_cfg[tname] = (int(w), int(h))
            except Exception as e:
                print(f"ERROR: Failed to resolve size for texture '{tname}': {e}")
                raise

        # Store FS config in the widget
        self.glcl_widget.set_render_config(self.buffer_data, self.render_pipeline_info)
        self.glcl_widget.set_fs_config(textures_cfg, fbos_cfg, fs_pipeline)

        # After OGLSystem.clear() the programs are gone; force recompile & bake if context is valid
        self.glcl_widget.rebuild_gl_resources()
        
    def bake_kernels(self, config):
        print("Baking OpenCL kernels...")
        self.baked_kernel_calls.clear()
        kernels_config = config.get("kernels", {})
        pipeline = config.get("kernel_pipeline", [])
        parameters = config.get("parameters", {})

        for kernel_name in pipeline:
            if kernel_name not in kernels_config:
                raise ValueError(f"Kernel '{kernel_name}' from pipeline not found in kernels configuration.")
            
            local_size, global_size_expr, buffer_names, param_names = kernels_config[kernel_name]
            
            kernel_obj = self.ocl_system.get_kernel(kernel_name)
            
            global_size = (self._resolve_expression(global_size_expr, parameters),)
            
            args_list = []
            for buf_name in buffer_names:
                args_list.append(self.ocl_system.get_buffer(buf_name))
            
            for param_name in param_names:
                value, type_str, _ = parameters[param_name]
                if type_str == "int":
                    args_list.append(np.int32(value))
                elif type_str == "float":
                    args_list.append(np.float32(value))
                else:
                    args_list.append(value)
            
            baked_call = BakedKernelCall(self.ocl_system, kernel_obj, global_size, local_size, tuple(args_list))
            self.baked_kernel_calls.append(baked_call)

    def _resolve_expression(self, expr, parameters):
        if isinstance(expr, str):
            if expr in parameters:
                return parameters[expr][0]
            else:
                raise KeyError(f"Parameter '{expr}' not found in configuration.")
        return int(expr)

    def create_parameter_controls(self, config):
        print("Creating parameter controls...")
        parameters = config.get("parameters", {})
        for name, (value, type_str, step) in parameters.items():
            container = QWidget()
            layout = QHBoxLayout(container)
            layout.setContentsMargins(0, 0, 0, 0)
            label = QtWidgets.QLabel(name)
            
            if type_str == "int":
                widget = QtWidgets.QSpinBox()
                widget.setRange(-1_000_000, 1_000_000)
                widget.setSingleStep(int(step))
                widget.setValue(int(value))
                widget.valueChanged.connect(lambda v, p=name: self.update_parameter(p, v))
            elif type_str == "float":
                widget = QtWidgets.QDoubleSpinBox()
                widget.setRange(-1_000_000.0, 1_000_000.0)
                widget.setDecimals(6)
                widget.setSingleStep(float(step))
                widget.setValue(float(value))
                widget.valueChanged.connect(lambda v, p=name: self.update_parameter(p, v))

            layout.addWidget(label)
            layout.addWidget(widget)
            self.params_layout.addWidget(container)
            self.param_widgets[name] = widget

    def start_simulation(self):
        self.simulation_running = True
        self.debug_frame_counter = 0  # DEBUG reset on (re)start
        interval = self.frame_interval_ms
        self.sim_timer.start(interval)
        print(f"Simulation started (interval={interval} ms, fps~={(1000.0/interval) if interval>0 else 'max'}).")

    def pause_simulation(self):
        self.simulation_running = False
        self.sim_timer.stop()
        print("Simulation paused.")

    def toggle_simulation(self):
        # Use timer state as the source of truth to avoid desync
        if self.sim_timer.isActive():
            print("Toggle: pausing simulation (timer active)")
            self.pause_simulation()
        else:
            print("Toggle: starting simulation (timer inactive)")
            self.start_simulation()

    def reset_simulation(self):
        if self.current_script:
            print("Resetting simulation...")
            self.load_and_apply_script(self.current_script)

    def update_simulation(self):
        self.debug_frame_counter += 1
        try:
            if not self.bDebugGL:
                for kernel_call in self.baked_kernel_calls:
                    kernel_call.execute()
                
                for buf_name in self.buffers_to_sync:
                    gl_obj = self.glcl_widget.gl_objects.get(buf_name)
                    if gl_obj:
                        host_data = np.empty(self.buffer_shapes[buf_name], dtype=np.float32)
                        self.ocl_system.fromGPU(buf_name, host_data)
                        self.glcl_widget.update_buffer_data(buf_name, host_data)
            
            self.glcl_widget.update()
            
            if self.bDebugCL and self.debug_frame_counter <= self.nDebugFrames:
                for buf_name in self.buffer_shapes:
                    host = np.empty(self.buffer_shapes[buf_name], dtype=np.float32)
                    self.ocl_system.fromGPU(buf_name, host)
                    print(f"[Frame {self.debug_frame_counter}] '{buf_name}': {host[0]}")

            # Stop conditions
            # 1) Explicit frame limit (applies always if >0)
            if self.frame_limit > 0 and self.debug_frame_counter >= self.frame_limit:
                self.pause_simulation()
                print(f"Max frame limit ({self.frame_limit}) reached. Simulation stopped.")
                return
            # 2) Debug frame limit (only when in debug modes)
            if (self.bDebugCL or self.bDebugGL) and self.nDebugFrames > 0 and self.debug_frame_counter >= self.nDebugFrames:
                self.pause_simulation()
                print(f"Debug frame limit ({self.nDebugFrames}) reached. Simulation stopped.")
        except Exception as e:
            self.on_exception(e)

    def on_exception(self, e):
        print(f"FATAL SIMULATION ERROR: {e}")
        import traceback
        traceback.print_exc()
        self.simulation_running = False
        self.sim_timer.stop()
        # Use os._exit to prevent PyQt from catching the exception and hanging
        os._exit(1)

    def update_parameter(self, param_name, new_value):
        if self.current_config is None: return
        
        old_value, type_str, step = self.current_config["parameters"][param_name]
        
        # If particle count changes, a full reset is needed.
        if param_name == "particle_count":
            print(f"Particle count changed to {new_value}. Resetting simulation.")
            self.current_config["parameters"][param_name] = (new_value, type_str, step)
            self.reset_simulation()
        else:
            self.current_config["parameters"][param_name] = (new_value, type_str, step)
            self.bake_kernels(self.current_config)
            print(f"Parameter '{param_name}' updated to {new_value}. Kernels re-baked.")

    def _precompute_buffer_sync_list(self, config):
        print("Pre-computing buffer synchronization list...")
        render_buffers = {pass_info[2] for pass_info in self.render_pipeline_info if len(pass_info) > 2}
        
        kernel_buffers = set()
        for kernel_info in config.get("kernels", {}).values():
            kernel_buffers.update(kernel_info[2])
            
        self.buffers_to_sync = list(render_buffers.intersection(kernel_buffers))
        print(f"  Buffers to sync between CL and GL: {self.buffers_to_sync}")

if __name__ == '__main__':
    # -------------------------------------------------------------------------
    # Usage examples (development vs production):
    #
    # Production/interactive (continuous animation, ~60 FPS):
    #   python -m python.GLCL2.GLCLBrowser --script python.GLCL2/scripts/nbody.py --frame-delay 16
    #
    # Development: CL debugging (prints first element of each buffer for N frames and stops):
    #   python -m python.GLCL2.GLCLBrowser --script python.GLCL2/scripts/nbody.py --debug-cl --debug-frames 10 --frame-delay 100
    #
    # Development: GL-only (skip CL kernel execution; useful to test rendering/shaders):
    #   python -m python.GLCL2.GLCLBrowser --script python.GLCL2/scripts/nbody.py --debug-gl --start-paused
    #   # then press "Start/Pause" in the UI
    #
    # Production with explicit frame cap (e.g., batch run for 1000 frames):
    #   python -m python.GLCL2.GLCLBrowser --script python.GLCL2/scripts/nbody.py --max-frames 1000 --frame-delay 0
    # -------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description='GLCL Browser')
    parser.add_argument('--script', type=str, help='Path to simulation script')
    parser.add_argument('--debug-cl', action='store_true', help='Enable OpenCL debugging prints (applies debug-frames)')
    parser.add_argument('--debug-gl', action='store_true', help='Skip OpenCL execution for OpenGL debugging (applies debug-frames)')
    parser.add_argument('--debug-frames', type=int, default=0, help='Number of frames to run in debug mode (only when --debug-cl/--debug-gl active); 0 means unlimited')
    parser.add_argument('--fps', type=float, default=60.0, help='Deprecated: target frames per second; overridden by --frame-delay if provided')
    parser.add_argument('--frame-delay', type=float, default=None, help='Frame delay in milliseconds (0 for max speed). Takes precedence over --fps')
    parser.add_argument('--start-paused', action='store_true', help='Do not auto-start the simulation after loading the script')
    parser.add_argument('--max-frames', type=int, default=0, help='Explicit frame cap; 0 means unlimited')
    args = parser.parse_args()
    
    app = QApplication(sys.argv)
    
    default_script = os.path.join(os.path.dirname(__file__), "scripts", "nbody.py")
    script_path = args.script or default_script

    if not os.path.exists(script_path):
        print(f"Error: Script not found at {script_path}")
        sys.exit(1)

    print(f"Starting GLCL Browser with script: {script_path}")
    browser = GLCLBrowser(
        python_script_path=script_path,
        bDebugCL=args.debug_cl,
        bDebugGL=args.debug_gl,
        nDebugFrames=args.debug_frames,
        fps=args.fps,
        start_paused=args.start_paused,
        frame_limit=args.max_frames,
        frame_delay_ms=args.frame_delay
    )
    browser.show()
    sys.exit(app.exec_())