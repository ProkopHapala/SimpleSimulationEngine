"""
BaseGUI.py - Common GUI utilities and widget templates for PyQt5 applications
"""
from PyQt5 import QtWidgets
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtGui import QFont
from PyQt5 import QtCore, QtGui, QtWidgets



import json
import re

# --- helper to extract balanced block ---
def extract_json_block(src, start_pos, open_sym, close_sym, bToText=True ):
    #print("extract_json_block(): ", start_pos, open_sym, close_sym, bToText)
    open_idx = src.find(open_sym, start_pos)
    if open_idx == -1: return None
    depth = 0
    for i in range(open_idx, len(src)):
        if   src[i] == open_sym:  depth += 1
        elif src[i] == close_sym: depth -= 1
        if depth == 0:
            if bToText:
                b_start, b_end = open_idx, i
                text = src[b_start+1:b_end].strip('\n\r ') if b_start is not None else ''  
                text = "\n".join([ l.strip() for l in text.split('\n') ])
                print("extract_json_block() extracted text:\n", text)
                return text
            else:
                return (open_idx,i)
    print("extract_json_block() failed to extract block")
    return None

def strip_json_comments(json_str):
    """Remove // and /* */ style comments from JSON string"""
    pattern = r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"'
    return re.sub(
        pattern,
        lambda m: m.group(0) if m.group(0).startswith(('"', "'")) else '',
        json_str,
        flags=re.MULTILINE|re.DOTALL
    )

class BaseGUI(QtWidgets.QMainWindow):
    """Base class containing common GUI utilities and widget templates"""
    
    def __init__(self, title="Application GUI"):
        super().__init__()
        # set smaller default font
        app = QtWidgets.QApplication.instance()
        if app:
            app.setFont(QFont("Sans", 8))
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle(title)
        self.main_widget = QtWidgets.QWidget(self)
        self.setCentralWidget(self.main_widget)

    
    def button(self, text, callback=None, tooltip=None, layout=None):
       #print("button() ", text, "with callback", callback, "tooltip", tooltip, "layout", layout)
        btn = QtWidgets.QPushButton(text)
        if callback is not None: btn.clicked.connect(callback)
        if tooltip  is not None: btn.setToolTip(tooltip)
        if layout   is not None: layout.addWidget(btn)
        return btn

    def checkBox(self, text, checked=False, callback=None, layout=None):
        chk = QtWidgets.QCheckBox(text)
        chk.setChecked(checked)
        if callback is not None: chk.stateChanged.connect(callback)
        if layout   is not None: layout.addWidget(chk)
        return chk

    def comboBox(self, items=None, callback=None, layout=None):
        cb = QtWidgets.QComboBox()
        if items is not None:    cb.addItems(items)
        if callback is not None: cb.currentIndexChanged.connect(callback)
        if layout   is not None: layout.addWidget(cb)
        return cb

    def spinBox(self, value=0.0, step=0.1, max_width=80, vmin=-1e9, vmax=1e9, decimals=4, layout=None, label=None):
        spin = QtWidgets.QDoubleSpinBox()
        spin.setDecimals(decimals)
        spin.setSingleStep(step)
        spin.setRange(vmin, vmax)
        spin.setValue(value)
        spin.setMaximumWidth(max_width)
        spin.valueChanged.connect(self.on_param_changed)
        if layout is not None:
            if label is not None: 
                layout.addRow(label, spin)
            else:
                layout.addWidget(spin)
        return spin
    
    def textEdit(self, read_only=False, layout=None, wrap=False):
        txt = QtWidgets.QTextEdit()
        txt.setReadOnly(read_only)
        if not wrap: txt.setWordWrapMode(QtGui.QTextOption.NoWrap)
        if layout is not None: layout.addWidget(txt)
        return txt

    def spin_row(self,defaults, step, layout=None, label=None):
        container = QtWidgets.QWidget()
        hbox = QtWidgets.QHBoxLayout(container)
        hbox.setContentsMargins(0,0,0,0)
        hbox.setSpacing(2)
        spins = [self.spinBox(d, step, 70, layout=hbox) for d in defaults]
        #for spin in spins: hbox.addWidget(spin)
        self.params_layout.addRow(label, container)
        #if layout: layout.addWidget(container)
        return spins

    # def populate_params_from_json(self, params_dict):
    #     """Create spin boxes from *params_dict*"""
    #     # Should be implemented by child class
    #     raise NotImplementedError("populate_params_from_json must be implemented by child class")

    def populate_params_from_json(self, params_dict):
        """Create spin boxes from *params_dict* and show lines in txt_uniforms."""
        # Clear existing widgets and layout rows
        print("---------------\npopulate_params_from_json()")
        while self.params_layout.rowCount() > 0: self.params_layout.removeRow(0)
        self.param_widgets.clear()
        for name, (typ,defaults,step) in params_dict.items():
            print("name: ", name, "typ: ", typ, "defaults: ", defaults, "step: ", step)
            if len(defaults) == 1:
                print("single value: ", name, defaults, step)
                self.param_widgets[name] = self.spinBox(defaults, step, layout=self.params_layout, label=name)
            else:
                print("multiple values: ", name, defaults, step)
                self.param_widgets[name] = self.spin_row(defaults, step, layout=self.params_layout, label=name)
        #exit()
        self.update_sim_uniforms()
        print("populate_params_from_json() DONE")