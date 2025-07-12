"""GLSL_GUI.py – minimal PyQt5 GUI around *GLSL_Simulation*.


run as:
python -m pySymGLSL.GLSL_GUI


The layout is purposely simple:

---------------------------------------------------------
|                GLRenderWidget (QOpenGLWidget)         |
---------------------------------------------------------
|  Uniforms list  |  Render-pipeline editor  |  Params  |
---------------------------------------------------------

Buttons:
• "Reload Shaders & Bake" – compiles shaders from disk and re-bakes the graph
  defined in the pipeline editor.
• "Update Uniforms" – regenerates parameter widgets from the *uniform list*.
• "Auto" checkbox – when checked, changing any parameter immediately triggers
  `run_graph` and an update of the GL view.

This is only a *starting point* – it demonstrates wiring GUI↔simulation; the
user can refine ergonomics later.
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List

from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtOpenGL import QGLFormat, QGLWidget
from PyQt5.QtWidgets import QFileDialog, QMessageBox

import moderngl
import numpy as np
import json
import re

from .GLSL_Simulation import GLSL_Simulation, _DEFAULT_VS

# --- helper to extract balanced block ---
def extract_json_block(src, start_pos, open_sym, close_sym):
    open_idx = src.find(open_sym, start_pos)
    if open_idx == -1: return None, None
    depth = 0
    for i in range(open_idx, len(src)):
        if src[i] == open_sym: depth += 1
        elif src[i] == close_sym: depth -= 1
        if depth == 0:
            return open_idx, i
    return None, None

# -----------------------------------------------------------------------------
# OpenGL rendering widget
# -----------------------------------------------------------------------------
class GLRenderWidget(QGLWidget):
    def __init__(self, parent=None, *, sim_size=(512, 512)):
        fmt = QGLFormat()
        fmt.setVersion(3, 3)
        fmt.setProfile(QGLFormat.CoreProfile)
        fmt.setSampleBuffers(False)
        QGLWidget.__init__(self, fmt, parent)

        self.sim_size = sim_size
        self.sim: GLSL_Simulation | None = None
        self.baked_graph = []
        self.dynamic_values: Dict[str, float] = {}
        # Disabled continuous redraw; UI will call updateGL on demand
        # self.timer = QtCore.QTimer()
        # self.timer.timeout.connect(self.updateGL)
        # self.timer.start(16)

    # ------------------------------------------------------------------
    def initializeGL(self):
        self.sim = GLSL_Simulation(self.sim_size, ctx=moderngl.create_context())
        self.ctx = self.sim.ctx  # convenience
        # fullscreen display program sampling uTex
        disp_fs = """#version 330 core
in vec2 v_texcoord;
uniform sampler2D uTex;
out vec4 fragColor;
void main(){ fragColor = texture(uTex, v_texcoord);}"""
        self.display_prog = self.ctx.program(vertex_shader=_DEFAULT_VS, fragment_shader=disp_fs)
        # create separate VAO for on-screen display (simulation's quad_vao is kept for internal passes)
        self.display_vao = self.ctx.vertex_array(self.display_prog, self.sim._quad_content)
        self.display_prog['uTex'].value = 0
        self.ctx.clear(0.0, 0.0, 0.0, 1.0)

    def paintGL(self):
        if self.baked_graph:
            # ensure default resolution uniform
            if 'iResolution' not in self.dynamic_values:
                w, h = self.sim_size
                self.dynamic_values['iResolution'] = (float(w), float(h), 1.0)
            self.sim.run_graph(self.baked_graph, self.dynamic_values)
            # blit latest output (first texture) to screen for now
            display_name = self.parent().cb_display.currentText() if hasattr(self.parent(), 'cb_display') else ''
            tex = self.sim.textures.get(display_name) if display_name else None
            if tex is None:
                tex = next(iter(self.sim.textures.values()))
            tex.use(location=0)
            self.ctx.screen.use()
            self.display_vao.render(moderngl.TRIANGLE_STRIP)
        else:
            self.ctx.clear(0.1, 0.1, 0.1)

    def resizeGL(self, w: int, h: int):
        if self.sim:
            self.ctx.viewport = (0, 0, w, h)

    # ------------------------------------------------------------------
    # External API for GUI
    # ------------------------------------------------------------------
    def rebuild(self, graph_lines: List[str]):
        """Recompile shaders from disk & bake graph from *graph_lines*.

        Each *graph_line* syntax:
          program_name fragment_path output_tex input0 input1 | uniform1 uniform2
        parts are whitespace-separated; a '|' separates input names from uniforms.
        """
        if not self.sim:
            return
        # Clear previous resources
        self.sim.release()
        self.sim = GLSL_Simulation(self.sim_size, ctx=self.ctx)

        passes = []
        for line in graph_lines:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # parse
            try:
                before_u, *after_u = line.split("|")
                tokens = before_u.strip().split()
                if len(tokens) < 3:
                    raise ValueError("Need at least program, frag_path, output")
                prog_name, frag_path, output = tokens[:3]
                inputs = tokens[3:]
                uniforms = after_u[0].strip().split() if after_u else []
            except Exception as exc:
                print("Error parsing line:", line, exc)
                continue

            frag_path = Path(frag_path).expanduser()
            self.sim.load_program(prog_name, fragment_path=frag_path)
            passes.append((prog_name, output, inputs, uniforms))

        self.baked_graph = self.sim.bake_graph(passes)
        self.updateGL()

    def update_uniforms(self, values: Dict[str, float]):
        self.dynamic_values.update(values)
        self.updateGL()


# -----------------------------------------------------------------------------
# Main window
# -----------------------------------------------------------------------------

class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        # this is relative to the script location
        #fname="pipelines/fluid_new.json"
        fname="pipelines/gauss.json"
        this_dir = Path(__file__).parent
        self.default_pipeline = this_dir / fname
        self.setWindowTitle("pySymGLSL GUI")
        self.resize(1200, 800)

        # Widgets -------------------------------------------------------------
        self.gl_view = GLRenderWidget(self)
        

        self.txt_uniforms = QtWidgets.QPlainTextEdit()
        self.txt_uniforms.setPlaceholderText("uniform_name min max default … one per line")
        btn_box1 = QtWidgets.QHBoxLayout()
        self.btn_update_uniforms = QtWidgets.QPushButton("Update Uniforms")
        self.btn_load_uniforms = QtWidgets.QPushButton("Load")
        self.btn_save_uniforms = QtWidgets.QPushButton("Save")
        btn_box1.addWidget(self.btn_update_uniforms)
        btn_box1.addWidget(self.btn_load_uniforms)
        btn_box1.addWidget(self.btn_save_uniforms)

        self.txt_pipeline = QtWidgets.QPlainTextEdit()
        self.txt_pipeline.setPlaceholderText("pipeline line syntax: prog frag_path output [inputs…] | [uniforms…]")
        self.txt_pipeline.setWordWrapMode(QtGui.QTextOption.NoWrap)
        btn_box2 = QtWidgets.QHBoxLayout()
        self.btn_rebake = QtWidgets.QPushButton("Reload & Bake")
        self.btn_load_pipeline = QtWidgets.QPushButton("Load")
        self.btn_save_pipeline = QtWidgets.QPushButton("Save")
        btn_box2.addWidget(self.btn_rebake)
        btn_box2.addWidget(self.btn_load_pipeline)
        btn_box2.addWidget(self.btn_save_pipeline)

        self.params_panel = QtWidgets.QWidget()
        self.params_layout = QtWidgets.QFormLayout(self.params_panel)
        self.chk_auto = QtWidgets.QCheckBox("auto")
        self.chk_auto.setChecked(True)

        self.cb_display = QtWidgets.QComboBox()  # choose texture to display

        # Layout --------------------------------------------------------------
        side_layout = QtWidgets.QVBoxLayout()
        side_layout.addWidget(QtWidgets.QLabel("Uniform list:"))
        side_layout.addWidget(self.txt_uniforms, 1)
        side_layout.addLayout(btn_box1)
        side_layout.addWidget(QtWidgets.QLabel("Render pipeline:"))
        side_layout.addWidget(self.txt_pipeline, 2)
        side_layout.addLayout(btn_box2)
        side_layout.addWidget(QtWidgets.QLabel("Parameters:"))
        side_layout.addWidget(QtWidgets.QLabel("Display texture:"))
        side_layout.addWidget(self.cb_display)
        side_layout.addWidget(QtWidgets.QLabel("Parameters:"))
        side_layout.addWidget(self.params_panel, 2)
        side_layout.addWidget(self.chk_auto)

        main_layout = QtWidgets.QHBoxLayout(self)
        main_layout.addWidget(self.gl_view, 4)
        main_layout.addLayout(side_layout, 1)

        # Signals -------------------------------------------------------------
        self.btn_rebake.clicked.connect(self.on_rebake)
        self.btn_update_uniforms.clicked.connect(self.on_update_uniforms)
        self.btn_load_pipeline.clicked.connect(self.on_load_json_pipeline)
        self.btn_save_pipeline.clicked.connect(self.on_save_pipeline)
        self.btn_load_uniforms.clicked.connect(self.on_load_uniforms)
        self.btn_save_uniforms.clicked.connect(self.on_save_uniforms)
        self.chk_auto.stateChanged.connect(self.on_auto_toggle)

        # internal
        # Schedule default pipeline load after GL context is ready
        QtCore.QTimer.singleShot(0, self._init_default)

        self.param_widgets: Dict[str, QtWidgets.QDoubleSpinBox] = {}

        

    # ------------------------------------------------------------------
    def _init_default(self):
        """Initialize GL simulation and load default pipeline after startup"""
        # Ensure simulation context is ready
        self.gl_view.initializeGL()
        # Load default pipeline
        self.load_pipeline(self.default_pipeline)

    def parse_uniform_lines(self):
        """Parse uniforms text box content as JSON and return parameters dict."""
        params_text = self.txt_uniforms.toPlainText().strip()
        #if not params_text:  return {}
        # Wrap in braces if not already a JSON object
        if not (params_text.startswith('{') and params_text.endswith('}')):
            params_text = '{' + params_text + '}'
        return json.loads(params_text)

    def on_update_uniforms(self):
        print("on_update_uniforms")
        for w in self.param_widgets.values():
            if isinstance(w, list):
                for spin in w: spin.deleteLater()
            else:
                w.deleteLater()
        self.param_widgets.clear()
        params_dict = self.parse_uniform_lines()
        self.populate_params_from_json(params_dict)
        self.update_sim_uniforms()

    def on_param_changed(self):
        print("on_param_changed()")
        if self.chk_auto.isChecked():
            self.update_sim_uniforms()

    def update_sim_uniforms(self):
        """Collect current parameter widget values (scalars or vectors) and send them to GL widget."""
        vals = {}
        for name, widget in self.param_widgets.items():
            if isinstance(widget, list):
                # vector uniform – gather components in original order
                vals[name] = [w.value() for w in widget]
            else:
                vals[name] = widget.value()
        self.gl_view.update_uniforms(vals)

    def on_rebake(self):
        """Re-compile shaders & bake graph using the contents of the two editors."""
        print("on_rebake()")
        params_raw   = self.txt_uniforms.toPlainText().strip()
        pipeline_raw = self.txt_pipeline.toPlainText().strip()
        # Build minimal JSON text from editors
        json_text = (
            '{ "parameters": {' + params_raw + '}, "Pipeline": [' + pipeline_raw + '] }'
        )
        try:
            data = json.loads(json_text)
        except Exception as exc:
            QMessageBox.critical(self, "JSON error", f"Failed to parse editors content as JSON:\n{exc}")
            return

        # Rebuild simulation using the new pipeline
        shader_dir = Path(__file__).parent / "shaders"
        try:
            baked, tex_names = self.gl_view.sim.build_pipeline(data["Pipeline"], shader_dir)
        except Exception as exc:
            QMessageBox.critical(self, "Bake error", f"Failed to build pipeline:\n{exc}")
            return

        self.gl_view.baked_graph = baked

        # Update display combo
        self.cb_display.clear()
        self.cb_display.addItems(tex_names)
        if tex_names:
            self.cb_display.setCurrentIndex(0)

        # (Re)create parameter widgets from current parameters
        self.populate_params_from_json(data.get("parameters", {}))
        self.gl_view.updateGL()


        # ------------------------------------------------------------------
    # Pipeline JSON helper
    # ------------------------------------------------------------------
    def load_pipeline_json(self, json_path: str | Path):
        """Load pipeline description from *json_path*.

        Returns a tuple (baked_graph, parameters:list[str], texture_names:list[str])
        """
        import re
        raw = Path(json_path).read_text()
        # strip // comments
        raw = re.sub(r'//.*$', '', raw, flags=re.MULTILINE)
        # remove trailing commas before } or ]
        raw = re.sub(r',\s*(?=[}\]])', '', raw)
        data = json.loads(raw)
        base_dir = Path(json_path).parent.parent / "shaders"
        baked_pipeline, tex_names= self.gl_view.sim.build_pipeline(data["Pipeline"], base_dir)        
        return baked_pipeline, data["parameters"], tex_names


    def load_pipeline(self, fname):
        # Read raw JSON file
        raw = Path(fname).read_text()
        
        # Parse to get data structure for simulation
        try:
            data = json.loads(raw)
            baked, params_dict, tex_names = self.load_pipeline_json(fname)
            self.gl_view.baked_graph = baked
            
            # Update display combo
            self.cb_display.clear()
            self.cb_display.addItems(tex_names)
            if tex_names:
                self.cb_display.setCurrentIndex(0)
            
            # Extract raw sections for display preserving original formatting
            try:
                # parameters { ... }
                key_pos = raw.find('"parameters"')
                p_start, p_end = extract_json_block(raw, key_pos, '{', '}')
                params_text = raw[p_start+1:p_end].strip('\n\r').strip() if p_start is not None else ''
                # Pipeline [ ... ]
                key_pos2 = raw.find('"Pipeline"')
                b_start, b_end = extract_json_block(raw, key_pos2, '[', ']')
                pipeline_text = raw[b_start+1:b_end].strip('\n\r ') if b_start is not None else ''                
                # Set text boxes
                self.txt_uniforms.setPlainText(params_text)
                self.txt_pipeline.setPlainText(pipeline_text)
            except Exception as e:
                print("Error extracting raw sections:", e)
                # Fallback to pretty-printed version
                self.txt_pipeline.setPlainText(json.dumps(data.get("Pipeline", []), indent=2))
                
            # Still need to populate param widgets from parsed data
            self.populate_params_from_json(params_dict)
            self.gl_view.updateGL()
        except Exception as e:
            print("Failed to load pipeline", e)

    # --- JSON pipeline loader ------------------------------------------
    def on_load_json_pipeline(self):
        print("on_load_json_pipeline")
        fname, _ = QFileDialog.getOpenFileName(self, 'Load JSON Pipeline', '', 'JSON (*.json)')
        if not fname:
            return
        self.load_pipeline(fname)

    def create_spin_box(self, value=0.0, step=0.1, max_width=80, vmin=-1e9, vmax=1e9, decimals=4):
        spin = QtWidgets.QDoubleSpinBox()
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setRange(vmin, vmax)
        spin.setValue(value)
        spin.setMaximumWidth(max_width)
        spin.valueChanged.connect(self.on_param_changed)
        return spin

    def populate_params_from_json(self, params_dict):
        """Create spin boxes from *params_dict* and show lines in txt_uniforms."""
        # Clear existing widgets and layout rows
        while self.params_layout.rowCount() > 0:
            self.params_layout.removeRow(0)
        self.param_widgets.clear()
        
        for name, spec in params_dict.items():
            print("name: ", name, "spec: ", spec)
            # Determine arity, defaults, step
            if isinstance(spec, str):
                typ = spec; cnt = 1 if typ == 'float' else int(typ[-1]); defaults = [0.0]*cnt; step = 0.1
            elif isinstance(spec, list) and len(spec) >= 3:
                typ, defaults, step = spec[0], spec[1], spec[2]; cnt = 1 if typ == 'float' else int(typ[-1])
            # Single value
            if cnt == 1:
                print("single value: ", name, defaults, step)
                spin = self.create_spin_box(defaults[0], step)
                self.params_layout.addRow(name, spin)
                self.param_widgets[name] = spin
            # Multiple values
            else:
                print("multiple values: ", name, defaults, step)
                container = QtWidgets.QWidget()
                hbox = QtWidgets.QHBoxLayout(container)
                hbox.setContentsMargins(0,0,0,0)
                hbox.setSpacing(2)
                spins = [self.create_spin_box(d, step, 70) for d in defaults]
                for spin in spins: hbox.addWidget(spin)
                self.params_layout.addRow(name, container)
                self.param_widgets[name] = spins
        self.update_sim_uniforms()

    def refresh_display_combo(self):
        if self.gl_view.sim:
            self.cb_display.clear()
            self.cb_display.addItems(self.gl_view.sim.textures.keys())

    def on_load_pipeline(self):
        print("on_load_pipeline")
        fname, _ = QFileDialog.getOpenFileName(self, 'Load Pipeline', '', 'Text (*.txt)')
        if fname:
            with open(fname) as f:
                self.txt_pipeline.setPlainText(f.read())

    def on_save_pipeline(self):
        print("on_save_pipeline")
        fname, _ = QFileDialog.getSaveFileName(self, 'Save Pipeline', '', 'JSON (*.json)')
        if fname:
            try:
                # Get raw sections from text boxes
                params_text   = self.txt_uniforms.toPlainText().strip()
                pipeline_text = self.txt_pipeline.toPlainText().strip()
                json_text = f'{{\n    "parameters": {{\n        {params_text}\n    }},\n    "Pipeline": [\n        {pipeline_text}\n    ]\n}}'
                print("saving json_text: ", json_text)
                json.loads(json_text)
                with open(fname, 'w') as f: f.write(json_text)
            except Exception as e:
                QMessageBox.critical(self, 'Error', f'Failed to save pipeline: {str(e)}')

    def on_load_uniforms(self):
        print("on_load_uniforms")
        fname, _ = QFileDialog.getOpenFileName(self, 'Load Uniforms', '', 'Text (*.txt)')
        if fname:
            with open(fname) as f:
                self.txt_uniforms.setPlainText(f.read())

    def on_save_uniforms(self):
        print("on_save_uniforms")
        fname, _ = QFileDialog.getSaveFileName(self, 'Save Uniforms', '', 'Text (*.txt)')
        if fname:
            with open(fname, 'w') as f:
                f.write(self.txt_uniforms.toPlainText())

    def on_auto_toggle(self, state):
        print("on_auto_toggle", state)
        if state and self.param_widgets:
            self.update_sim_uniforms()


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    mw = MainWindow()
    mw.show()
    sys.exit(app.exec_())
