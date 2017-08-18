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
from PyQt4 import QtGui, QtCore, QtOpenGL
from PyQt4.QtOpenGL import QGLWidget,QGLFormat
from PyQt4.QtCore import QTimer
# PyOpenGL imports
import OpenGL.GL as gl
import OpenGL.arrays.vbo as glvbo

import GLUtils as glu

import sys
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
    width, height = 800, 800

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


    def initializeGL(self):
        """Initialize OpenGL, VBOs, upload data on the GPU, etc."""
        gl.glClearColor(0, 0, 0, 0)          # background color

        self.points = np.array([[-1.0,-1.0],[1.0,-1.0],[-1.0,1.0],   [1.0,1.0],[1.0,-1.0],[-1.0,1.0]], dtype=np.float32)

        self.vbo = glvbo.VBO(self.points)    # create a Vertex Buffer Object with the specified data
        vs = glu.compile_vertex_shader  (self.vertex_code)       # compile the vertex shader
        fs = glu.compile_fragment_shader(self.fragment_code)     # compile the fragment shader
        self.shader = glu.link_shader_program(vs, fs)            # compile the vertex shader
        gl.glUseProgram(self.shader)
        '''
        uniform vec3      iResolution;           // viewport resolution (in pixels)
        uniform float     iTime;                 // shader playback time (in seconds)
        uniform float     iTimeDelta;            // render time (in seconds)
        uniform int       iFrame;                // shader playback frame
        uniform vec4      iMouse;                // mouse pixel coords. xy: current (if MLB down), zw: click
        uniform vec4      iDate;                 // (year, month, day, time in seconds)
        uniform float     iSampleRate;           // sound sample rate (i.e., 44100)
        //uniform float     iChannelTime[4];       // channel playback time (in seconds)
        //uniform vec3      iChannelResolution[4]; // channel resolution (in pixels)
        //uniform samplerXX iChannel0..3;          // input channel. XX = 2D/Cube
        '''
        self.UNIFORM_LOCATIONS = { 
            'iResolution':  gl.glGetUniformLocation( self.shader, 'iResolution' ),
            'iTime':        gl.glGetUniformLocation( self.shader, 'iTime'       ), 
            'iTimeDelta':   gl.glGetUniformLocation( self.shader, 'iTimeDelta'  ), 
            'iFrame':       gl.glGetUniformLocation( self.shader, 'iFrame'      ), 
            'iMouse':       gl.glGetUniformLocation( self.shader, 'iMouse'      ), 
            'iDate':        gl.glGetUniformLocation( self.shader, 'iDate'       ), 
            'iSampleRate':  gl.glGetUniformLocation( self.shader, 'iSampleRate' ), 
        }
        print self.UNIFORM_LOCATIONS
        self.frameCount = 0

        # Set up a timer to call updateGL() every 0 ms

    def paintGL(self):
        #print "self.frameCount = ", self.frameCount
        #try:
        """Paint the scene."""
        #print "QCursor.pos () : ", QtGui.QCursor.pos()
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)    # clear the buffer
        self.vbo.bind()                       # bind the VBO
        gl.glEnableVertexAttribArray(0)
        gl.glVertexAttribPointer(0, 2, gl.GL_FLOAT, gl.GL_FALSE, 0, None)   # these vertices contain 2 single precision coordinates
        gl.glUseProgram( self.shader)
        mousePos = self.mapFromGlobal( QtGui.QCursor.pos() )
        gl.glUniform1f ( self.UNIFORM_LOCATIONS['iTime'],       0.01*self.frameCount ) 
        gl.glUniform3f ( self.UNIFORM_LOCATIONS['iResolution'], 1.0*self.width, 1.0*self.height, 1.0 )
        gl.glUniform4f ( self.UNIFORM_LOCATIONS['iMouse'],      0.1*mousePos.x(), 0.1*mousePos.y(), 0.0, 0.0 )
        gl.glDrawArrays(gl.GL_TRIANGLES, 0, len(self.points))   # draw "count" points from the VBO
        self.frameCount +=1

    def resizeGL(self, width, height):
        """Called upon window resizing: reinitialize the viewport."""
        self.width, self.height = width, height   # update the window size
        gl.glViewport(0, 0, width, height)        # paint within the whole window

class MainWindow(QtGui.QWidget):
    def __init__(self):
        super(MainWindow, self).__init__()
        qgl_format = QGLFormat()
        qgl_format.setSwapInterval(1)

        self.widget = GLPlotWidget(qgl_format)
        self.widget.vertex_code   = glu.DEFAULT_VERTEX_SHADER
        #self.widget.fragment_code = glu.DEFAULT_FRAGMENT_SHADER

        '''
        with open("Torus_intersection.glslf","r") as f:
            #self.widget.fragment_code=f.read()
            self.widget.fragment_code=glu.fromShaderToy(f.read())
            #print self.widget.fragment_code
            with open("fragment_code.glslf", "w") as fo:
                fo.write(self.widget.fragment_code)
        '''

        name, desc, codes = glu.downloadFromShaderToy('Xds3zN')
        self.widget.fragment_code=glu.fromShaderToy(codes[0].encode('utf8'))
        with open("fragment_code.glslf", "w") as fo: fo.write(self.widget.fragment_code)

        self.button = QtGui.QPushButton('Test', self)

        mainLayout = QtGui.QHBoxLayout()
        mainLayout.addWidget(self.widget)
        mainLayout.addWidget(self.button)

        self.setLayout(mainLayout)

if __name__ == '__main__':
    app = QtGui.QApplication(['Yo'])
    window = MainWindow()
    window.show()
    app.exec_()
