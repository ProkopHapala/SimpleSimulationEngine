#!/usr/bin/python

# https://matplotlib.org/examples/user_interfaces/index.html
# https://matplotlib.org/examples/user_interfaces/embedding_in_qt5.html
# embedding_in_qt5.py --- Simple Qt5 application embedding matplotlib canvases

from __future__ import unicode_literals
import sys
import os
import numpy as np
import matplotlib; matplotlib.use('Qt5Agg')
from PyQt5 import QtCore, QtWidgets, QtGui
from numpy import arange, sin, pi
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import sys
#import pyRay       as ra
#import pyRay.scene as scn 
#import pyRay.ocl   as ocl 
#import pyRay.image as img 

import common      as ra
import scene as scn 
import ocl   as ocl 
import image as img 

import time

test_scenes = {

"user_func1" :
('''
U( obj1( pos, 0.5f, 0.6f ) );
S( obj1( pos+(float3)(0.0f,0.0f,+0.5f), 0.3f, 0.5 ) );
S( obj1( pos+(float3)(0.0f,0.0f,-0.5f), 0.3f, 0.7 ) );

''',
'''
float obj1(float3 pos, float R, float L){
float dist = 1000000.0f;
U( sdSphere  ( pos- (float3)(+L,0.0f,0.0), R )  );
U( sdSphere  ( pos- (float3)(-L,0.0f,0.0), R )  );
U( sdSphere  ( pos- (float3)(0.0f,+L,0.0), R )  );
U( sdSphere  ( pos- (float3)(0.0f,-L,0.0), R )  );
return dist;
}
'''),

"dirRep" :
('''
dist = sdSphere( opRepDir( pos, (float3)(0.7f,0.0f,0.0f), 5.0f ), 0.35f );
''',
""),

"objRep" :
('''
U( obj1( opRepDir( pos, (float3)(0.0f,0.0f,1.0f), 5.0f ), 0.5f, 0.6f ) );
pos.z+=0.5;
pos  = opRepDir( pos, (float3)(0.0f,0.0f,1.0f), 5.0f );
S( obj1( pos, 0.3f, 0.6f ) );
S( sdSphere  ( pos, 0.5 ) );
''',
'''
float obj1(float3 pos, float R, float L){
float dist = 1000000.0f;
U( sdSphere  ( pos- (float3)(+L,0.0f,0.0), R )  );
U( sdSphere  ( pos- (float3)(-L,0.0f,0.0), R )  );
U( sdSphere  ( pos- (float3)(0.0f,+L,0.0), R )  );
U( sdSphere  ( pos- (float3)(0.0f,-L,0.0), R )  );
return dist;
}
'''),

"granade" : 
(''' 
U( sdSphere  ( pos- (float3)(0.0f,0.5f,0.0), 0.5f )  );
U( sdCone    ( pos, (float3)(0.8f,0.6f,0.3f ) ) );
U( sdCylinder( pos, (float2)(0.1f,0.2f)       ) );
S( sdTorus   ( pos- (float3)(0.0f,0.1f,0.0), (float2)(0.25f,0.1f)     ) );
S( sdPlane   ( pos, (float3)(1.0,0.0,0.0), 0.0 ) ); 
S( sdCylinder( pos, (float2)(0.05f,0.9f)       ) );
S( sdSphere  ( pos- (float3)(0.0f,0.5f,0.0), 0.3f )  );
''',
""),

"RepRot" : 
(''' 
dist = opU( dist, sdSphere  ( opRepRot( pos, 6.0f )-(float3)(0.18f,0.0f,0.0f), 0.1f ) );
dist = opS( dist, sdTorus   ( pos, (float2)(0.25f,0.07f) ) );
dist = opU( dist, sdCone    ( pos, (float3)(0.8f,0.6f,0.3f) ) );
dist = opU( dist, sdCylinder( pos, (float2)(0.1f,0.2f) ) );
''',
""),

"primitives" :
(''' 
pos.y*=-1.0f;
pos.y+=0.5f;
U(sdPlaneY(pos));
dist = opU( dist, sdSphere(    pos-(float3)( 0.0,0.25, 0.0), 0.25 ) );
dist = opU( dist, sdBox(       pos-(float3)( 1.0,0.25, 0.0), (float3)(0.25) ) );
dist = opU( dist, udRoundBox(  pos-(float3)( 1.0,0.25, 1.0), (float3)(0.15), 0.1 ) );
dist = opU( dist, sdTorus(     pos-(float3)( 0.0,0.25, 1.0), (float2)(0.20,0.05) ) );
dist = opU( dist, sdCapsule(   pos,(float3)(-1.3,0.10,-0.1), (float3)(-0.8,0.50,0.2), 0.1  ) );
dist = opU( dist, sdTriPrism(  pos-(float3)(-1.0,0.25,-1.0), (float2)(0.25,0.05) ) );
dist = opU( dist, sdCylinder(  pos-(float3)( 1.0,0.30,-1.0), (float2)(0.1,0.2) ) );
dist = opU( dist, sdCone(      pos-(float3)( 0.0,0.50,-1.0), (float3)(0.8,0.6,0.3) ) );
dist = opU( dist, sdTorus82(   pos-(float3)( 0.0,0.25, 2.0), (float2)(0.20,0.05) ) );
dist = opU( dist, sdTorus88(   pos-(float3)(-1.0,0.25, 2.0), (float2)(0.20,0.05) ) );
dist = opU( dist, sdCylinder6( pos-(float3)( 1.0,0.30, 2.0), (float2)(0.1,0.2) ) );
dist = opU( dist, sdHexPrism(  pos-(float3)(-1.0,0.20, 1.0), (float2)(0.25,0.05) ) );
dist = opU( dist, sdPryamid4(  pos-(float3)(-1.0,0.15,-2.0), (float3)(0.8,0.6,0.25) ) );
dist = opU( dist, opS( udRoundBox(  pos-(float3)(-2.0,0.2, 1.0), (float3)(0.15),0.05), sdSphere(    pos-(float3)(-2.0,0.2, 1.0), 0.25)) );
//dist = opU( dist, opS( sdTorus82(  pos-(float3)(-2.0,0.2, 0.0), (float2)(0.20,0.1)),
//                           sdCylinder(  opRep( (float3)(atan(pos.x+2.0,pos.z)/6.2831, pos.y, 0.02+0.5*length(pos-(float3)(-2.0,0.2, 0.0))), (float3)(0.05,1.0,0.05)), (float2)(0.02,0.6))) );
dist = opU( dist, 0.5*sdSphere(    pos-(float3)(-2.0,0.25,-1.0), 0.2 ) + 0.03*sin(50.0*pos.x)*sin(50.0*pos.y)*sin(50.0*pos.z) );
dist = opU( dist, 0.5*sdTorus( opTwist(pos-(float3)(-2.0,0.25, 2.0)),(float2)(0.20,0.05)) );
dist = opU( dist, sdConeSection( pos-(float3)( 0.0,0.35,-2.0), 0.15, 0.2, 0.1 ) );
dist = opU( dist, sdEllipsoid( pos-(float3)( 1.0,0.35,-2.0), (float3)(0.15, 0.2, 0.05) ) );
dist = opU( dist, sdEllipsoid( pos-(float3)( 1.0,0.35,-2.0), (float3)(0.15, 0.2, 0.05) ) );
''',
""),

"bounding_sphere" : 
('''
dist = 0.2;
float bound = sdSphere( pos, 0.8f );
if( 0.2f>bound  ){
dist = opU( dist, sdSphere( pos-(float3)(0.7f,0.0f,0.0f), 0.4f ) ); 
dist = opU( dist, sdSphere( pos-(float3)(0.0f,0.7f,0.0f), 0.4f ) ); 
dist = opU( dist, sdSphere( pos-(float3)(0.0f,0.0f,0.7f), 0.4f ) ); 
dist = opU( dist, sdSphere( pos-(float3)(-0.7f,0.0f,0.0f), 0.4f ) ); 
dist = opU( dist, sdSphere( pos-(float3)(0.0f,-0.7f,0.0f), 0.4f ) ); 
dist = opU( dist, sdSphere( pos-(float3)(0.0f,0.0f,-0.7f), 0.4f ) ); 
dist = opI(dist, bound);
}''',
""),

}

class MyDynamicMplCanvas(FigureCanvas):
    """A canvas that updates itself every second with a new plot."""

    cbar = None 
    
    def __init__(self, parent=None, width=8, height=8, dpi=100 ):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        #self.compute_initial_figure()
        FigureCanvas.__init__(self, self.fig )
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,QtWidgets.QSizePolicy.Expanding,QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
            
    def plot(self, hitp, hitn, light_dir, rd, view_type ):
        self.axes.cla()
        if   view_type == 'iters':
            self.cax = self.axes.imshow( hitn[:,:,3] ,          interpolation='nearest' )
        elif view_type == 'AO':    
            AO = 1/(1+0.02*hitn[:,:,3]**2) 
            self.cax = self.axes.imshow( AO*0.5+0.5, interpolation='nearest', cmap='gray' )
        elif view_type == 'depth':
            mask         = hitp[:,:,3] > 20.0 
            hitp[mask,3] = np.NaN
            self.cax = self.axes.imshow( np.log( hitp[:,:,3] ), interpolation='nearest' )
        elif view_type == 'normal':
            self.cax = self.axes.imshow( hitn[:,:,0:3]*0.5+0.5, interpolation='nearest' )
        elif view_type == 'diffuse':
            diffuse  = img.getDiffuse ( hitn, light_dir )
            self.cax = self.axes.imshow( diffuse*0.5+0.5,       interpolation='nearest', cmap='gray' )
        elif view_type == 'phong':
            diffuse  = img.getDiffuse ( hitn, light_dir )
            specular = img.getSpecular( hitn, light_dir, rd, gloss=256.0, power=2 )
            self.cax = self.axes.imshow( (diffuse+specular)*0.5+0.5, interpolation='nearest', cmap='gray' )
        self.draw()

class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
       
        self.font = QtGui.QFont()
        self.font.setFamily('Lucida')
        self.font.setFixedPitch(True)
        self.font.setPointSize(10)
       
        # --- init QtMain
        QtWidgets.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("application main window")
        self.main_widget = QtWidgets.QWidget(self)
        l = QtWidgets.QVBoxLayout(self.main_widget)
        self.mplc1 = MyDynamicMplCanvas(self.main_widget, width=5, height=4, dpi=100)
        l.addWidget(self.mplc1)
        
        # --- bxPhi
        scaleLabel = QtWidgets.QLabel("Phi"); l.addWidget( scaleLabel )
        self.bxPhi = QtWidgets.QSpinBox()
        self.bxPhi.setRange(-360.0, 360.0)
        self.bxPhi.setSingleStep(15.0)
        self.bxPhi.setValue(0.0)
        self.bxPhi.valueChanged.connect(self.updateCamera)
        l.addWidget( self.bxPhi )
        
        # --- bxTheta
        scaleLabel = QtWidgets.QLabel("theta"); l.addWidget( scaleLabel )
        self.bxTheta = QtWidgets.QSpinBox()
        self.bxTheta.setRange(-90.0, 90.0)
        self.bxTheta.setSingleStep(15.0)
        self.bxTheta.setValue(0.0)
        self.bxTheta.valueChanged.connect(self.updateCamera)
        l.addWidget( self.bxTheta )
        
        # --- slView
        self.slView = QtWidgets.QComboBox(self)
        self.slView.addItem("iters")
        self.slView.addItem("AO")
        self.slView.addItem("depth")
        self.slView.addItem("normal")
        self.slView.addItem("diffuse")
        self.slView.addItem("phong")
        self.slView.setCurrentIndex(3)
        self.slView.currentIndexChanged.connect(self.updateView)
        l.addWidget( self.slView )
        
        # --- slScene
        self.slScene = QtWidgets.QComboBox(self)
        self.slScene.addItem("user_func1")
        self.slScene.addItem("dirRep")
        self.slScene.addItem("objRep") 
        self.slScene.addItem("RepRot")
        self.slScene.addItem("granade")
        self.slScene.addItem("primitives")
        self.slScene.addItem("bounding_sphere")
        self.slScene.setCurrentIndex(1)
        self.slScene.currentIndexChanged.connect(self.selectScene)
        l.addWidget( self.slScene )
        
        # TODO /home/prokop/Dropbox/MyDevSW/Python/_GUI/PyQt5/complex/syntaxhighlighter.py
        self.txtScene = QtWidgets.QPlainTextEdit(self)
        self.txtScene.setFont(self.font)
        #self.txtScene.returnPressed.connect(self.updateScene)
        l.addWidget( self.txtScene )
        
        self.txtFunctions = QtWidgets.QPlainTextEdit(self)
        self.txtFunctions.setFont(self.font)
        #self.txtScene.returnPressed.connect(self.updateScene)
        l.addWidget( self.txtFunctions )
        
        # --- btScene
        self.btScene = QtWidgets.QPushButton('render scene', self)
        self.btScene.setToolTip('compile scene program and run it')
        self.btScene.clicked.connect(self.updateScene)
        l.addWidget( self.btScene )
        
        # --- btSave
        self.btSave = QtWidgets.QPushButton('save fig', self)
        self.btSave.setToolTip('save current figure')
        self.btSave.clicked.connect(self.saveFig)
        l.addWidget( self.btSave )
        
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        
        self.light_dir  = np.array([0.5,-1.0,0.0]); 
        self.light_dir /=np.sqrt(np.dot(self.light_dir,self.light_dir));
        
        #ra.screen_resolution[0] = 16
        #ra.screen_resolution[1] = 16
        ra.screen_resolution[0] = 512
        ra.screen_resolution[1] = 512
        #ra.defines["HIT_PREC"]  = 0.01
        #ra.defines["MAX_STEPS"] = 16

        self.updateBuffers()
        #self.updateScene()
        self.selectScene()
        self.updateCamera()
    
    def updateView(self):
        view_kind = str(self.slView.currentText())
        #print view_kind
        t1 = time.clock() 
        self.mplc1.plot( self.hitp, self.hitn, self.light_dir, self.rd, view_kind )
        t2 = time.clock(); print "time(mplc1.plot) %f " %(t2-t1)
    
    def updateRender(self):
        t1 = time.clock()
        self.hitp,self.hitn = ocl.run_rayTraceBasic ( self.rd, self.ro, self.kargs )
        t2 = time.clock(); print "time(ocl.run_rayTraceBasic) %f " %(t2-t1)
        self.updateView()

    def updateCamera(self):
        phi   = self.bxPhi  .value()*(np.pi/180.0)
        theta = self.bxTheta.value()*(np.pi/180.0)
        fw = np.array([np.sin(phi)*np.cos(theta),np.sin(theta),np.cos(phi)*np.cos(theta)])
        #print theta,phi, fw
        self.cam          = ra.getCamMat( fw )
        t1 = time.clock()
        self.rd, self.ro  = ocl.getRays ( self.cam, fw0=-4.0, t0=8.0, tg=(0.25,0.25) )
        t2 = time.clock(); print "time(ocl.getRays) %f " %(t2-t1)
        self.updateRender()
   
    def updateScene(self):
        src_scene     = self.txtScene.toPlainText()
        src_functions = self.txtFunctions.toPlainText()
        print src_scene
        print src_functions
        prog_scr = scn.makeProgram( src_scene, src_user_func=src_functions )
        ocl.make_program(prog_scr)
        try:
            self.rd.shape # check if rd exists 
            self.updateRender()
        except:
            print "rd not initialized"
            pass
            
    def selectScene(self):
        scene_name = str(self.slScene.currentText())
        src_scene = test_scenes[scene_name]
        self.txtScene    .setPlainText(src_scene[0])
        self.txtFunctions.setPlainText(src_scene[1])
        self.updateScene()
            
    def updateBuffers(self):
        # TODO : we have to delete old buffers before creation new one
        self.kargs     = ocl.prep_rayTraceBasic( )
    
    def saveFig(self):
        self.mplc1.fig.savefig('render.png',bbox_inches='tight')
    
if __name__ == "__main__":

    #E,lvec, nDim, head = GU.loadXSF('ELJ_cl.xsf' );
    #import matplotlib.pyplot as plt
    #plt.imshow( E[50,:,:] )
    #plt.show()
    
    qApp = QtWidgets.QApplication(sys.argv)
    aw = ApplicationWindow()
    aw.show()
    sys.exit(qApp.exec_())
    
    
    

