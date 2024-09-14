from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtOpenGL import QGLWidget, QGLFormat
from PyQt5.QtCore import QTimer,QElapsedTimer

import moderngl
import numpy as np
import sys, os
from FormulaNode import *

'''
# Example of more complex scene
P1 = Plane(  [0.0,0.0,0.0] )
S2 = Sphere( [0.0,0.25, 0.0], 0.25 )
S3 = Sphere( [0.0,0.25    , 0.0], 0.2 )
S4 = Sphere( [0.0,0.25-0.2, 0.2], 0.2 )
S5 = Sphere( [-0.2,0.25+0.2, 0.0], 0.2 )  
scene = P1 + (S2 -( S3 +S4+S5 ))
'''




# Default Vertex Shader as a string
DEFAULT_VERTEX_SHADER = """
#version 330
in vec2 in_vert;
out vec2 v_text;
void main() {
    gl_Position = vec4(in_vert, 0.0, 1.0);
    v_text = in_vert;
}
"""

SHADERTOY_HEADER = """
#version 330
uniform vec3  iResolution;
uniform float iTime;
uniform float iTimeDelta;
uniform int   iFrame;
uniform vec4  iMouse;
uniform vec4  iDate;
uniform float iSampleRate;

layout(std140) uniform VecData {
    vec3 u_vec_data[100];  // Example: an array of 100 vec3s (adjust the size as needed)
};
"""
SHADERTOY_MAIN="""
void main(){
    vec2 fragCoord = gl_FragCoord.xy;
    vec4 fragColor; 
    mainImage( fragColor, fragCoord);
    gl_FragColor = fragColor;
}
"""

class GLPlotWidget(QGLWidget):
    # default window size
    width, height = 256, 256

    def __init__(self, parent):
        super(GLPlotWidget, self).__init__(parent)
        self.setMinimumSize(self.width, self.height)
        self.frameCount = 0

        self.elapsed_timer = QElapsedTimer()
        self.elapsed_timer.start()

        self.timer = QTimer()
        self.timer.setInterval(10)
        #self.timer.timeout.connect(self.updateGL)
        self.timer.timeout.connect(self.update_frame)
        self.ctx = None
        self.shader_program = None
        self.vbo = None
        self.vao = None

    def update_frame(self):
        self.vec_data[0,1] = 0.25
        self.vec_data[0,2] = 0.25
        self.vec_data[0,0] = 2.0*np.sin( self.elapsed_timer.elapsed()*0.001 )
        #print( "vec_data[0] ", self.vec_data[0])
        self.update_vec_data()

        self.updateGL()
        #self.update()

    def updateShader(self):
        if self.shader_program:
            self.shader_program.release()
        self.shader_program = self.ctx.program( vertex_shader=DEFAULT_VERTEX_SHADER,  fragment_shader=self.fragment_code )
        self.shader_program['iResolution'].value = (self.width, self.height, 1.0)
        self.update()

    def update_vec_data(self):
        self.vec_buffer.write(self.vec_data.tobytes())
        self.update()

    def initializeGL(self):
        self.ctx = moderngl.create_context()
        points = np.array([[-1.0, -1.0], [1.0, -1.0], [-1.0, 1.0], [1.0, 1.0], [1.0, -1.0], [-1.0, 1.0]], dtype='f4')
        self.vbo = self.ctx.buffer(points)
        #self.vao = self.ctx.simple_vertex_array(self.shader_program, self.vbo, 'in_vert')
        self.updateShader()
        self.vao = self.ctx.vertex_array(self.shader_program, [(self.vbo, '2f', 'in_vert')])
        self.frameCount = 0

        # Example array of vec3 (can also be vec4, just adjust the shape and dtype)
        #self.vec_data = np.zeros( (100,4), dtype='f4' )
        self.vec_data = np.zeros( (100,4), dtype=np.float32 )

        # Create a buffer and upload the vec3 array data to the GPU
        self.vec_buffer = self.ctx.buffer(self.vec_data.tobytes())
        
        # Ensure that the buffer is bound to the correct uniform block
        uniform_block = self.shader_program['VecData']  #    ;print( "uniform_block", uniform_block)
        uniform_block.binding = 0
        #self.shader_program.uniform_block_binding(block_index, 0)  # Bind the block to binding point 0
        self.vec_buffer.bind_to_uniform_block(0) 

    def link_program(self,vertex_shader, fragment_shader):
        """Link vertex and fragment shaders into a program"""
        return self.ctx.program(vertex_shader=vertex_shader, fragment_shader=fragment_shader)

    def paintGL(self):
        self.ctx.clear(0.0, 0.0, 0.0)
        try:
            self.shader_program['iTime'].value = 0.01 * self.frameCount
        except:
            pass

        mousePos = self.mapFromGlobal(QtGui.QCursor.pos())
        #self.shader_program['iMouse'].value = (0.1 * mousePos.x(), 0.1 * mousePos.y(), 0.0, 0.0)
        # Check the state of the mouse buttons
        left_button_pressed = QtGui.QGuiApplication.mouseButtons()  & QtCore.Qt.LeftButton
        right_button_pressed = QtGui.QGuiApplication.mouseButtons() & QtCore.Qt.RightButton
    
        iMouse_z = 1.0 if left_button_pressed else 0.0
        iMouse_w = 1.0 if right_button_pressed else 0.0
        H = self.height
        
        try:
            # Pass the mouse position and button states to the shader
            self.shader_program['iMouse'].value = (
                mousePos.x(),  # x component: scaled mouse x position
                H-mousePos.y(),  # y component: scaled mouse y position
                iMouse_z,            # z component: LMB state (1.0 if pressed, 0.0 if not)
                iMouse_w             # w component: RMB state (1.0 if pressed, 0.0 if not)
            )
            #print( "iMouse ", self.shader_program['iMouse'].value )
        except Exception as e:
            pass
            #print(f"Error setting iMouse: {e}")

        self.vao.render(moderngl.TRIANGLES)
        self.frameCount += 1

    def resizeGL(self, width, height):
        self.width, self.height = width, height
        self.ctx.viewport = (0, 0, width, height)
        self.shader_program['iResolution'].value = (self.width, self.height, 1.0)

class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        super(MainWindow, self).__init__()
        
        with open('primitives.glslf', 'r')         as f: self.src_primitives = f.read()
        with open('sdRayMarchRenderer.glslf', 'r') as f: self.src_renderer   = f.read()

        self.font =w= QtGui.QFont(); w.setFamily('Lucida'); w.setFixedPitch(True); w.setPointSize(10)

        l0 = QtWidgets.QHBoxLayout();
        self.setLayout(l0)

        qgl_format = QGLFormat()
        qgl_format.setSwapInterval(1)
        self.glview =w= GLPlotWidget(qgl_format); l0.addWidget(w); w.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        
        l1 = QtWidgets.QVBoxLayout(); l0.addLayout(l1)
        
        
        self.txtCode    =w= QtWidgets.QTextEdit(self);             l1.addWidget(w); w.setFont(self.font);                    w.setLineWrapMode(QtWidgets.QTextEdit.NoWrap)      
        l2 = QtWidgets.QHBoxLayout(); l1.addLayout(l2)

        self.btShUpdate =w= QtWidgets.QPushButton('update', self); l2.addWidget(w); w.setToolTip('update shader from code'); w.clicked.connect(self.updateShader);
        self.btShSave   =w= QtWidgets.QPushButton('save', self);   l2.addWidget(w); w.setToolTip('save shader from code');   w.clicked.connect(self.saveShaderDlg)
        self.btShLoad   =w= QtWidgets.QPushButton('load', self);   l2.addWidget(w); w.setToolTip('load shader from code');   w.clicked.connect(self.loadShaderDlg);
        self.btPlay     =w= QtWidgets.QPushButton('play', self);   l2.addWidget(w); w.setToolTip('on/off update on timer');  w.clicked.connect(self.play);

        # Create the dropdown list (QComboBox)
        self.cbMode =w= QtWidgets.QComboBox(self)
        w.addItems(["Plain", "SDF 2D", "SDF 3D", "RayTrace 3D"])
        #self.combo_box.setCurrentIndex(2)
        w.setCurrentText("SDF 3D")
        w.currentIndexChanged.connect(self.selectMode)
        l2.addWidget(w);

        self.cbDir  =w= QtWidgets.QComboBox(self); w.currentIndexChanged.connect(self.update_cbFile)
        self.cbFile =w= QtWidgets.QComboBox(self);

        # Set up the layout and main window
        #l2.addWidget(QLabel("Select Subfolder:"))
        l2.addWidget(self.cbDir)
        #l2.addWidget(QLabel("Select .glslf File:"))
        l2.addWidget(self.cbFile)

        self.setGeometry(100, 100, 1600, 1000) 
        self.setWindowTitle('Shader Toy')

        #self.initScene_SDF3D_Python()
        self.initScene_Plain()

        self.load_subfolders()
        self.update_cbFile()

        self.cbFile.currentIndexChanged.connect(self.load_selected_shader)

    def load_subfolders(self):
        # Load the subfolder names from the ./shaders directory
        shaders_dir = "./shaders"
        if os.path.exists(shaders_dir) and os.path.isdir(shaders_dir):
            subfolders = [f.name for f in os.scandir(shaders_dir) if f.is_dir()]
            self.cbDir.addItems(subfolders)
        else:
            print(f"Directory {shaders_dir} not found.")
    
    def update_cbFile(self):
        self.cbFile.clear()
        selected_folder = self.cbDir.currentText()
        #print(f"update_file_combo() Selected folder: {selected_folder}")
        if not selected_folder:
            return
        folder_path = os.path.join("./shaders", selected_folder)
        #print(f"update_file_combo() folder_path: {folder_path}")
        if os.path.exists(folder_path) and os.path.isdir(folder_path):
            glsl_files = [f for f in os.listdir(folder_path) if f.endswith('.glslf')]
            #print(f"update_file_combo() glsl_files: {glsl_files}")
            self.cbFile.addItems(glsl_files)

    def initScene_SDF3D_Python(self):
        #self.txtCode.setPlainText("d = opU( sdPlane(pos), sdSphere(pos-vec3( 0.0,0.25, 0.0), 0.25 ) );")
        #self.txtCode.setPlainText("d = opU( sdPlane(pos), sdSphere(pos-u_vec_data[0], 0.25 ) );")
        self.cbMode.setCurrentText("SDF 3D")
        scene_code= '''P1 = Plane(  [0.0,0.0,0.0] )\nS2 = Sphere( [0.0,0.25, 0.0], 0.25 )\nscene = P1 + S2'''
        self.txtCode.setPlainText(scene_code)
        # P1 = Plane(  [0.0,0.0,0.0] )
        # #P1 = Formula( f"sdPlane( pos)" )
        # S2 = Sphere( [0.0,0.25, 0.0], 0.25 )
        # scene = P1 + S2

        #S1 = Sphere( [1.3,0.5,.2], 1.5 )
        #S2 = Sphere( [0.3,0.5,2.2], 0.5 )
        #B1 = RoundBox( [-0.3,0.5,-2.2], [0.5,0.5,0.5], 1.5 )
        #scene = (B1 - S1) + S2
        #scene_str = "d="+scene.eval_CGS() + ";";  print( "scene_str: \n", scene_str )
        #self.txtCode.setPlainText( scene_str )
        self.setWindowTitle('Shader Toy')
        #self.setGeometry(100, 100, 1920, 1080)  # Set position (x, y) and size (width, height)
        self.setGeometry(100, 100, 1600, 1000) 
        self.makeShaderCode_SDF3D( self.txtCode.toPlainText() )

    
    def load_shader_Plain( self, filename ):
        scene_code = open(filename, 'r').read()
        self.txtCode.setPlainText(scene_code)
        self.makeShaderCode_Plain( self.txtCode.toPlainText() )
    
    def initScene_Plain(self):
        self.cbMode.setCurrentText("Plain")
        self.load_shader_Plain( 'shaders/SDF_2D/Moon.glslf' )

    def load_selected_shader(self):
        print(f"load_selected_shader()")
        filename =  os.path.join("./shaders", self.cbDir.currentText(), self.cbFile.currentText())
        self.load_shader_Plain( filename )
        self.glview.updateShader()

    def selectMode(self):
        selected_mode = self.cbMode.currentText()
        print(f"Selected Mode: {selected_mode}")

    def play(self):
        self.btPlay.setText("Stop")
        self.btPlay.clicked.connect(self.stop)
        self.glview.timer.start()
    
    def stop(self):
        self.btPlay.setText("Play")
        self.btPlay.clicked.connect(self.play)
        self.glview.timer.stop()

    def saveShaderDlg(self):
        fname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save Shader ... ', os.getcwd(), filter='*.glslf')
        if fname:
            print("saving shader to file:", fname)
            with open(fname, "w") as fo: fo.write(self.txtCode.toPlainText())
    
    def loadShaderDlg(self):
        fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Load Shader ... ', os.getcwd(), filter='*.glslf')
        if fname:
            print("loading shader from file:", fname)
            with open(fname, 'r') as f: self.txtCode.setPlainText(f.read())
            self.updateShader()

    def makeShaderCode_SDF3D(self, scene_code, bPython=True ):
        if bPython:
            local_context = {}
            try:
                exec(scene_code, globals(), local_context)
                if 'scene' in local_context:
                    scene = local_context['scene']
                    scene_code = "d="+scene.eval_CGS()+";"
                    print(f"makeShaderCode_SDF3D() Scene variable captured: {scene_code}")
                else:
                    print("makeShaderCode_SDF3D() Scene variable not found.")

            except Exception as e:
                print(f"makeShaderCode_SDF3D() Error executing code: {e}")

        #scene_code = self.txtCode.toPlainText()
        src_scene = "vec2 map( in vec3 pos ){\n float d;\n" + scene_code + "\n return vec2(d,1.0);\n}"
        src = SHADERTOY_HEADER + self.src_primitives + src_scene + self.src_renderer
        self.glview.fragment_code = src
        with open("fragment_code.glslf", "w") as fo: fo.write(self.glview.fragment_code)


    def makeShaderCode_Plain(self, scene_code ):
        src = SHADERTOY_HEADER + scene_code + SHADERTOY_MAIN
        self.glview.fragment_code = src
        with open("fragment_code.glslf", "w") as fo: fo.write(self.glview.fragment_code)

    def updateShader(self):
        

        selected_mode = self.cbMode.currentText()
        
        # Handle different cases based on the selected mode
        if selected_mode == "Plain":
            self.makeShaderCode_Plain()
            #print("updateShader(mode='RayTrace 3D'): NOT IMPLEMENTED")
        elif selected_mode == "SDF 2D":
            print("updateShader(mode='RayTrace 3D'): NOT IMPLEMENTED")
        elif selected_mode == "SDF 3D":
            self.makeShaderCode_SDF3D( self.txtCode.toPlainText() )
        elif selected_mode == "RayTrace 3D":
            print("updateShader(mode='RayTrace 3D'): NOT IMPLEMENTED")
        else:
            print("Unknown mode selected")

        

        self.glview.updateShader()

if __name__ == '__main__':
    app = QtWidgets.QApplication(['pyShaderToy'])
    window = MainWindow()
    window.show()
    app.exec_()
