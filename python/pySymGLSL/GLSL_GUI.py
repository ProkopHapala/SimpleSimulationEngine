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
from .BaseGUI import BaseGUI, extract_json_block, strip_json_comments




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
        # Name of the texture currently shown in the widget (controlled by MainWindow)
        self.display_tex_name: str | None = None
        self.sim: GLSL_Simulation | None = None
        self.baked_graph = []
        self.dynamic_values: Dict[str, float] = {}
        # Timer for continuous redraw (disabled by default)
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.updateGL)


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

    # ------------------------------------------------------------------
    # Render -----------------------------------------------------------------
    def paintGL(self):
        if self.baked_graph:
            # ensure default resolution uniform
            # ensure default uniforms
            if 'iResolution' not in self.dynamic_values:
                w, h = self.sim_size
                self.dynamic_values['iResolution'] = (float(w), float(h), 1.0)
            # always provide frame counter
            self.dynamic_values['iFrame'] = self.sim.iFrame
            self.sim.run_graph(self.baked_graph, self.dynamic_values)
            # Choose texture to display. If *display_tex_name* was not set explicitly
            # (e.g. before GUI is fully initialized) fall back to the first texture.
            name = self.display_tex_name or ''
            tex = self.sim.textures.get(name) if name else None
            if tex is None and self.sim.textures:
                tex = next(iter(self.sim.textures.values()))
            tex.use(location=0)
            self.ctx.screen.use()
            self.ctx.viewport = (0, 0, self.width(), self.height())   # <- add this
            self.display_vao.render(moderngl.TRIANGLE_STRIP)
        else:
            self.ctx.clear(0.1, 0.1, 0.1)

    # ------------------------------------------------------------------
    # Play/Pause control ------------------------------------------------
    # ------------------------------------------------------------------
    # Helpers ------------------------------------------------------------
    def set_display(self, name: str | None):
        """Set the name of the texture that should be displayed on screen."""
        self.display_tex_name = name

    # Play/Pause control ------------------------------------------------
    def set_play(self, play: bool):
        if play:
            self.timer.start(16)   # ~60 FPS
        else:
            self.timer.stop()

    def resizeGL(self, w: int, h: int):
        if self.sim:
            self.ctx.viewport = (0, 0, w, h)
            # Update the baked graph with the new one generated during resize
            #self.baked_graph = self.sim.resize(w, h)

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

#class MainWindow(QtWidgets.QWidget):
class MainWindow(BaseGUI):
    
    def __init__(self):
        super().__init__("pySymGLSL GUI")
        
        # this is relative to the script location
        fname="pipelines/fluid.json"
        #fname="pipelines/gauss.json"
        this_dir = Path(__file__).parent
        self.default_pipeline = this_dir / fname
        self.resize(1200, 800)
        
        # Widgets -------------------------------------------------------------
        self.gl_view = GLRenderWidget(self)
                
        self.params_panel  = QtWidgets.QWidget()
        self.params_layout = QtWidgets.QFormLayout(self.params_panel)
        
        self.cb_display = QtWidgets.QComboBox()  # choose texture to display

        # Layout --------------------------------------------------------------
        side_layout = QtWidgets.QVBoxLayout()
        side_layout.addWidget(QtWidgets.QLabel("Uniform list:"))
        
        #side_layout.addWidget(self.txt_uniforms, 2)
        self.txt_uniforms  = self.textEdit(read_only=True, layout=side_layout)
                
        side_layout.addWidget(QtWidgets.QLabel("Render pipeline:"))
        self.txt_pipeline  = self.textEdit(read_only=True,       layout=side_layout)
        #side_layout.addWidget(self.txt_pipeline, 2)
        
        btn_box2 = QtWidgets.QHBoxLayout()
        #self.btn_update_uniforms  = self.button("Update Uniforms", self.on_update_uniforms, layout=btn_box1)
        #self.btn_rebake           = self.button("Reload & Bake",    self.on_rebake,             layout=btn_box2)
        self.btn_rebuild          = self.button("rebuild",         self.on_rebuild,             layout=btn_box2)
        self.btn_load_pipeline    = self.button("Load",            self.on_load_json_pipeline, layout=btn_box2)
        self.btn_save_pipeline    = self.button("Save",            self.on_save_pipeline,      layout=btn_box2)
        side_layout.addLayout(btn_box2)

        display_box = QtWidgets.QHBoxLayout()
        display_box.addWidget(QtWidgets.QLabel("Display texture:"))
        self.cb_display        = self.comboBox(layout=display_box, callback=self.on_display_changed)
        self.chk_play          = self.checkBox("Play", False, lambda s: self.gl_view.set_play(bool(s)), layout=display_box)
        side_layout.addLayout(display_box)
        
        side_layout.addWidget(QtWidgets.QLabel("Parameters:"))
        side_layout.addWidget(self.params_panel, 2)
        self.chk_auto          = self.checkBox("auto", True, self.on_auto_toggle, layout=side_layout)
        main_layout = QtWidgets.QHBoxLayout(self.main_widget)
        main_layout.addWidget(self.gl_view, 4)
        main_layout.addLayout(side_layout, 1)

        # internal
        # Schedule default pipeline load after GL context is ready
        QtCore.QTimer.singleShot(0, self.deferred_gl_init)

        self.param_widgets: Dict[str, QtWidgets.QDoubleSpinBox] = {}
        
    # ------------------------------------------------------------------
    def deferred_gl_init(self):
        """Initialize GL simulation and load default pipeline after startup"""
        self.gl_view.initializeGL()     # Ensure simulation context is ready
        self.load_pipeline(self.default_pipeline)  # Load default pipeline

    def parse_uniform_lines(self):
        """Parse uniforms text box content as JSON and return parameters dict."""
        params_text = self.txt_uniforms.toPlainText().strip()
        #if not params_text:  return {}
        # Wrap in braces if not already a JSON object
        if not (params_text.startswith('{') and params_text.endswith('}')):
            params_text = '{' + params_text + '}'
        return json.loads(params_text)

    def on_rebuild(self):
        print("on_rebuild()")
        self.on_update_uniforms()
        self.on_rebake()

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
        print("update_sim_uniforms()")
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
            #QMessageBox.critical(self, "JSON error", f"Failed to parse editors content as JSON:\n{exc}")
            print("on_rebake(): Failed to parse editors content as JSON: ", exc)
            print("json.loads() failed: ", json_text)
            return

        # Rebuild simulation using the new pipeline
        shader_dir = Path(__file__).parent / "shaders"
        try:
            baked, tex_names = self.gl_view.sim.build_pipeline(data["Pipeline"], shader_dir)
        except Exception as exc:
            #QMessageBox.critical(self, "Bake error", f"Failed to build pipeline:\n{exc}")
            print("on_rebake():  Failed to build pipeline: ", exc)
            print("baked: ", baked)
            print("tex_names: ", tex_names)
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

    def load_pipeline(self, fname):
        # Read raw JSON file
        print("load_pipeline(): ", fname)
        raw = Path(fname).read_text()
        raw = strip_json_comments(raw)
        #print("raw: ", raw)
        data = json.loads(raw)
        params_dict = data["parameters"]
        #baked, params_dict, tex_names = self.load_pipeline_json(fname)
        base_dir = Path(fname).parent.parent / "shaders"
        baked, tex_names= self.gl_view.sim.build_pipeline(data["Pipeline"], base_dir)
        self.gl_view.baked_graph = baked
        # Update display combo
        self.cb_display.clear()
        self.cb_display.addItems(tex_names)
        if tex_names:
            self.cb_display.setCurrentIndex(0)
        try:
            self.txt_uniforms.setPlainText(extract_json_block(raw, raw.find('"parameters"'), '{', '}', True))
            self.txt_pipeline.setPlainText(extract_json_block(raw, raw.find('"Pipeline"'), '[', ']', True))
        except Exception as e:
            print("Error in load_pipeline() extracting raw sections:", e)
            self.txt_pipeline.setPlainText(json.dumps(data.get("Pipeline", []), indent=2))
        self.populate_params_from_json(params_dict)
        self.gl_view.updateGL()
 
    # --- JSON pipeline loader ------------------------------------------
    def on_load_json_pipeline(self):
        print("on_load_json_pipeline()")
        fname, _ = QFileDialog.getOpenFileName(self, 'Load JSON Pipeline', '', 'JSON (*.json)')
        if not fname: return
        self.load_pipeline(fname)

    def refresh_display_combo(self):
        if self.gl_view.sim:
            self.cb_display.clear()
            self.cb_display.addItems(self.gl_view.sim.textures.keys())
    
    def on_display_changed(self):
        name = self.cb_display.currentText()
        print("on_display_changed() ->", name)
        # Inform GL widget which texture to show and request repaint
        self.gl_view.set_display(name)
        self.gl_view.updateGL()

    def on_load_pipeline(self):
        print("on_load_pipeline()")
        fname, _ = QFileDialog.getOpenFileName(self, 'Load Pipeline', '', 'Text (*.txt)')
        if fname:
            with open(fname) as f:
                self.txt_pipeline.setPlainText(f.read())

    def on_save_pipeline(self):
        print("on_save_pipeline()")
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
                #QMessageBox.critical(self, 'Error', f'Failed to save pipeline: {str(e)}')
                print("on_save_pipeline(): Failed to save pipeline: ", e)
                print("json.loads() failed: ", json_text)

    def on_load_uniforms(self):
        print("on_load_uniforms()")
        fname, _ = QFileDialog.getOpenFileName(self, 'Load Uniforms', '', 'Text (*.txt)')
        if fname:
            with open(fname) as f:
                self.txt_uniforms.setPlainText(f.read())

    def on_save_uniforms(self):
        print("on_save_uniforms()")
        fname, _ = QFileDialog.getSaveFileName(self, 'Save Uniforms', '', 'Text (*.txt)')
        if fname:
            with open(fname, 'w') as f:
                f.write(self.txt_uniforms.toPlainText())

    def on_auto_toggle(self, state):
        print("on_auto_toggle()", state)
        if state and self.param_widgets:
            self.update_sim_uniforms()


# -----------------------------------------------------------------------------

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    mw = MainWindow()
    mw.show()
    sys.exit(app.exec_())
