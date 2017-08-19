#!/usr/bin/python

# from here : 
# http://cyrille.rossant.net/shaders-opengl/
# https://stackoverflow.com/questions/33201384/pyqt-opengl-drawing-simple-scenes
# from here : https://stackoverflow.com/questions/33201384/pyqt-opengl-drawing-simple-scenes

'''
TODO:
 - Mouse
 - textures iChannle
 - frameBuffer (multiplass rendering)

'''

# PyQt4 imports

from PyQt4          import QtGui, QtCore, QtOpenGL
from PyQt4.QtOpenGL import QGLWidget,QGLFormat
from PyQt4.QtCore   import QTimer
# PyOpenGL imports
import OpenGL.GL as gl
import OpenGL.arrays.vbo as glvbo

import GLUtils as glu

import sys,os
import numpy as np


# Window creation function.
def create_window(window_class):
    """Create a Qt window in Python, or interactively in IPython with Qt GUI
    event loop integration:
        # in ~/.ipython/ipython_config.py
        c.TerminalIPythonApp.gui = 'qt'
        c.TerminalIPythonApp.pylab = 'qt'
    See also:
        http://ipython.org/ipython-doc/dev/interactive/qtconsole.html#qt-and-the-qtconsole
    """
    app_created = False
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication(sys.argv)
        app_created = True
    app.references = set()
    window = window_class()
    app.references.add(window)
    window.show()
    if app_created:
        app.exec_()
    return window

class GLPlotWidget(QGLWidget):
    # default window size
    width, height = 256, 256

    def __init__(self, parent):
        QGLWidget.__init__(self, parent)
        self.setMinimumSize(self.width, self.height)
        self.frameCount = 0
        # Set up a timer to call updateGL() every 0 ms
        self.timer = QTimer()
        self.timer.setInterval(10)
        #self.timer.timeout.connect(self.updateWidget)
        self.timer.timeout.connect(self.updateGL)
        #timer.timeout.connect( self.widget.update)
        self.timer.start()

    def updateShader(self):
        try:gl.glDeleteShader(self.shader)
        except:pass
        self.vertex_code = glu.DEFAULT_VERTEX_SHADER
        vs = glu.compile_vertex_shader  (self.vertex_code)       # compile the vertex shader
        fs = glu.compile_fragment_shader(self.fragment_code)     # compile the fragment shader
        self.shader = glu.link_shader_program(vs, fs)            # compile the vertex shader
        gl.glUseProgram(self.shader)
        self.UNIFORM_LOCATIONS = { 
            'iResolution':  gl.glGetUniformLocation( self.shader, 'iResolution' ),
            'iTime':        gl.glGetUniformLocation( self.shader, 'iTime'       ), 
            'iTimeDelta':   gl.glGetUniformLocation( self.shader, 'iTimeDelta'  ), 
            'iFrame':       gl.glGetUniformLocation( self.shader, 'iFrame'      ), 
            'iMouse':       gl.glGetUniformLocation( self.shader, 'iMouse'      ), 
            'iDate':        gl.glGetUniformLocation( self.shader, 'iDate'       ), 
            'iSampleRate':  gl.glGetUniformLocation( self.shader, 'iSampleRate' ), 
        }
        #print self.UNIFORM_LOCATIONS
        gl.glUniform3f ( self.UNIFORM_LOCATIONS['iResolution'], 1.0*self.width, 1.0*self.height, 1.0 )

    def initializeGL(self):
        """Initialize OpenGL, VBOs, upload data on the GPU, etc."""
        self.points = np.array([[-1.0,-1.0],[1.0,-1.0],[-1.0,1.0],   [1.0,1.0],[1.0,-1.0],[-1.0,1.0]], dtype=np.float32)
        self.vbo = glvbo.VBO(self.points)    # create a Vertex Buffer Object with the specified data
        self.updateShader()
        gl.glClearColor(0, 0, 0, 0)          # background color
        self.frameCount = 0

    def paintGL(self):
        """Paint the scene."""
        #print "self.frameCount = ", self.frameCount
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)    # clear the buffer
        self.vbo.bind()                       # bind the VBO
        gl.glEnableVertexAttribArray(0)
        gl.glVertexAttribPointer(0, 2, gl.GL_FLOAT, gl.GL_FALSE, 0, None)   # these vertices contain 2 single precision coordinates
        gl.glUseProgram( self.shader)
        mousePos = self.mapFromGlobal( QtGui.QCursor.pos() )
        gl.glUniform1f ( self.UNIFORM_LOCATIONS['iTime'],       0.01*self.frameCount ) 
        #gl.glUniform3f ( self.UNIFORM_LOCATIONS['iResolution'], 1.0*self.width, 1.0*self.height, 1.0 )
        gl.glUniform4f ( self.UNIFORM_LOCATIONS['iMouse'],      0.1*mousePos.x(), 0.1*mousePos.y(), 0.0, 0.0 )
        gl.glDrawArrays(gl.GL_TRIANGLES, 0, len(self.points))   # draw "count" points from the VBO
        self.frameCount +=1

    def resizeGL(self, width, height):
        """Called upon window resizing: reinitialize the viewport."""
        self.width, self.height = width, height   # update the window size
        gl.glViewport(0, 0, width, height)        # paint within the whole window
        gl.glUniform3f ( self.UNIFORM_LOCATIONS['iResolution'], 1.0*self.width, 1.0*self.height, 1.0 )

class MainWindow(QtGui.QWidget):
    def __init__(self):
        super(MainWindow, self).__init__()
        qgl_format = QGLFormat()
        qgl_format.setSwapInterval(1)

        w = QtGui.QFont(); self.font=w; w.setFamily('Lucida'); w.setFixedPitch(True); w.setPointSize(10)

        l0 = QtGui.QHBoxLayout();  self.setLayout(l0)
        w = GLPlotWidget(qgl_format); self.glview = w; l0.addWidget(self.glview); w.setSizePolicy(QtGui.QSizePolicy.Expanding,QtGui.QSizePolicy.Expanding)
        
        l1 = QtGui.QVBoxLayout(); l0.addLayout(l1)
        w = QtGui.QTextEdit(self);            self.txtCode = w; l1.addWidget(w); w.setFont(self.font); w.setLineWrapMode(QtGui.QTextEdit.NoWrap)  

        l2 = QtGui.QHBoxLayout(); l1.addLayout(l2)
        w = QtGui.QPushButton('update',self); self.btShUpdate = w; l2.addWidget(w); w.setToolTip('update shader from code'); w.clicked.connect(self.updateShader); 
        w = QtGui.QPushButton('save',self);   self.btShSave   = w; l2.addWidget(w); w.setToolTip('save shader from code');   w.clicked.connect(self.saveShaderDlg); 
        w = QtGui.QPushButton('load',self);   self.btShLoad   = w; l2.addWidget(w); w.setToolTip('load shader from code');   w.clicked.connect(self.loadShaderDlg); 
        w = QtGui.QLineEdit(); self.bxShToy = w; l2.addWidget(w);  w.returnPressed.connect(self.downloadShaderEvent )
        #l3 = QtGui.QHBoxLayout(); l1.addLayout(l3)
        #w = QtGui.QLineEdit(); 
        #w = QtGui.QPushButton('fromShadeToy',self);   self.btSave  = w; l2.addWidget(w); w.setToolTip('load shader from code');   w.clicked.connect(self.loadShaderDlg);
        
        #self.loadShaderToy ("Torus_intersection.glslf")
        self.loadShader    ("Torus_intersection.glslf")
        #self.downloadShader("Xds3zN")

    def saveShaderDlg(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save Shader ... ', os.getcwd(), selectedFilter='*.glslf')
        if fname:
            print "saving shader to file : ", fname
            with open(fname, "w") as fo: fo.write(self.glview.fragment_code)
    
    def loadShaderDlg(self):
        fname= QtGui.QFileDialog.getOpenFileName(self, 'Load Shader ... ', os.getcwd(), selectedFilter='*.glslf')
        if fname:
            print "loading shader to file : ",fname
            self.loadShader(fname)
            self.updateShader()

    def updateShader(self):
        self.glview.fragment_code = unicode( self.txtCode.toPlainText()).encode('utf8');
        self.glview.updateShader()  

    def loadShader(self, fname ):
        w = self.glview
        with open(fname,"r") as f:w.fragment_code=f.read()
        self.txtCode.setPlainText( w.fragment_code )

    def downloadShader(self, id):
        w = self.glview 
        name, desc, codes = glu.downloadFromShaderToy(id)
        w.fragment_code=glu.fromShaderToy(codes[0].encode('utf8'))
        #with open("fragment_code.glslf", "w") as fo: fo.write(w.fragment_code)
        self.txtCode.setPlainText( w.fragment_code )
    
    def downloadShaderEvent(self):
        shId = str(self.bxShToy.text())
        fname = shId+".glslf"
        try:
            self.loadShade(fname)
        except:
            self.downloadShader( shId )
            with open(fname, "w") as fo: fo.write(self.glview.fragment_code)
        self.glview.updateShader( )

if __name__ == '__main__':
    app = QtGui.QApplication(['pyShaderToy'])
    window = MainWindow()
    window.show()
    app.exec_()
