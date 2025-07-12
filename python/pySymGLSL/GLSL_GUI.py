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
from PyQt5.QtWidgets import QFileDialog

import moderngl
import numpy as np
import json
import re

from .GLSL_Simulation import GLSL_Simulation, _DEFAULT_VS

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
        res = []
        for ln in self.txt_uniforms.toPlainText().splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            name = parts[0]
            default = float(parts[3]) if len(parts) > 3 else 0.0
            minv = float(parts[1]) if len(parts) > 1 else 0.0
            maxv = float(parts[2]) if len(parts) > 2 else 1.0
            res.append((name, minv, maxv, default))
        return res

    def on_update_uniforms(self):
        # Clear old
        for w in self.param_widgets.values():
            w.deleteLater()
        self.param_widgets.clear()

        for name, minv, maxv, default in self.parse_uniform_lines():
            spin = QtWidgets.QDoubleSpinBox()
            spin.setRange(minv, maxv)
            spin.setValue(default)
            spin.setDecimals(4)
            spin.valueChanged.connect(self.on_param_changed)
            self.params_layout.addRow(name, spin)
            self.param_widgets[name] = spin

        self.update_sim_uniforms()

    def on_param_changed(self):
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
        lines = self.txt_pipeline.toPlainText().splitlines()
        self.gl_view.rebuild(lines)
        self.refresh_display_combo()


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

        # # 1. Load shader programs (unique names)
        # for prog_name, (rel_path, _uniforms) in data.get("Shaders", {}).items():
        #     self.gl_view.sim.load_program(prog_name, fragment_path=base_dir / rel_path)
        
        baked_pipeline, tex_names= self.gl_view.sim.build_pipeline(data["Pipeline"], base_dir)        
        return baked_pipeline, data["parameters"], tex_names


    def load_pipeline(self, fname):
        # display Render pipeline JSON section
        try:
            import re, json as _json
            raw = Path(fname).read_text()
            raw = re.sub(r'//.*$', '', raw, flags=re.MULTILINE)
            raw = re.sub(r',\s*(?=[}\]])', '', raw)
            _data = _json.loads(raw)
            self.txt_pipeline.setPlainText(_json.dumps(_data.get("Pipeline", []), indent=2))
        except Exception as _exc:
            print("failed to show pipeline", _exc)
        # parse and build
        baked, params_dict, tex_names = self.load_pipeline_json(fname)
        self.gl_view.baked_graph = baked
        # populate param widgets
        self.populate_params_from_json(params_dict)
        # update display combo
        self.cb_display.clear()
        self.cb_display.addItems(tex_names)
        if tex_names:
            self.cb_display.setCurrentIndex(0)
        self.gl_view.updateGL()

    # --- JSON pipeline loader ------------------------------------------
    def on_load_json_pipeline(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Load JSON Pipeline', '', 'JSON (*.json)')
        if not fname:
            return
        self.load_pipeline(fname)

    def populate_params_from_json(self, params_dict):
        """Create spin boxes from *params_dict* and show lines in txt_uniforms."""
        for w in self.param_widgets.values():
            w.deleteLater()
        self.param_widgets.clear()
        lines = []  # for txt_uniforms

        for name, spec in params_dict.items():
            # Determine arity, defaults, step
            if isinstance(spec, str):
                typ = spec; cnt = 1 if typ == 'float' else int(typ[-1]); defaults = [0.0]*cnt; step = 0.1
            elif isinstance(spec, list) and len(spec) >= 3:
                typ, defaults, step = spec[0], spec[1], spec[2]; cnt = 1 if typ == 'float' else int(typ[-1])


            # build UI row -----------------------------------------------------
            if cnt == 1:
                spin = QtWidgets.QDoubleSpinBox(); spin.setDecimals(4); spin.setSingleStep(step); spin.setRange(-1e9, 1e9); spin.setValue(defaults[0]); spin.setMaximumWidth(80); spin.valueChanged.connect(self.on_param_changed)
                self.params_layout.addRow(name, spin)
                self.param_widgets[name] = spin
            else:
                container = QtWidgets.QWidget(); hbox = QtWidgets.QHBoxLayout(container); hbox.setContentsMargins(0,0,0,0); hbox.setSpacing(2)
                spins = []
                for i in range(cnt):
                    spin = QtWidgets.QDoubleSpinBox(); spin.setDecimals(4); spin.setSingleStep(step); spin.setRange(-1e9, 1e9); spin.setValue(defaults[i]); spin.setMaximumWidth(70); spin.valueChanged.connect(self.on_param_changed)
                    hbox.addWidget(spin); spins.append(spin)
                self.params_layout.addRow(name, container)
                self.param_widgets[name] = spins

            # Uniform list line
            line = f"{name} -inf inf " + " ".join(str(v) for v in defaults)
            lines.append(line)

        # show lists
        self.txt_uniforms.setPlainText("\n".join(lines))
        self.update_sim_uniforms()

    def refresh_display_combo(self):
        if self.gl_view.sim:
            self.cb_display.clear()
            self.cb_display.addItems(self.gl_view.sim.textures.keys())

    def on_load_pipeline(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Load Pipeline', '', 'Text (*.txt)')
        if fname:
            with open(fname) as f:
                self.txt_pipeline.setPlainText(f.read())

    def on_save_pipeline(self):
        fname, _ = QFileDialog.getSaveFileName(self, 'Save Pipeline', '', 'Text (*.txt)')
        if fname:
            with open(fname, 'w') as f:
                f.write(self.txt_pipeline.toPlainText())

    def on_load_uniforms(self):
        fname, _ = QFileDialog.getOpenFileName(self, 'Load Uniforms', '', 'Text (*.txt)')
        if fname:
            with open(fname) as f:
                self.txt_uniforms.setPlainText(f.read())

    def on_save_uniforms(self):
        fname, _ = QFileDialog.getSaveFileName(self, 'Save Uniforms', '', 'Text (*.txt)')
        if fname:
            with open(fname, 'w') as f:
                f.write(self.txt_uniforms.toPlainText())

    def on_auto_toggle(self, state):
        if state and self.param_widgets:
            self.update_sim_uniforms()


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    mw = MainWindow()
    mw.show()
    sys.exit(app.exec_())
