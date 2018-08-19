import numpy as np
from   ctypes import c_int, c_double, c_double, c_bool, c_float, c_char_p, c_void_p
import ctypes
import os

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )

def recompile(path):
    print( path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make clean" )
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )

def init():
    LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/libFlight')
    recompile(LIB_PATH_CPP)
    return ctypes.CDLL( LIB_PATH_CPP+"/libFlight.so" )

def initViewLib():
    LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs_SDL/FlightView')
    recompile(LIB_PATH_CPP)
    return ctypes.CDLL( LIB_PATH_CPP+"/libFlightView.so" )

# =========== main

lib     = init()
libView = None

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

#char_ptr   = ctypes.POINTER(c_char  )
double_ptr = ctypes.POINTER(c_double)

# ========= C functions

#void loadFromFile( char* fname )
lib.loadFromFile.argtypes = [c_char_p]
lib.loadFromFile.restype  = None
def loadFromFile( fname ):
    lib.loadFromFile( fname )

#void setPose( double* pos, double* vel, double* rot )
lib.setPose.argtypes = [array1d,array1d,array2d]
lib.setPose.restype  = None
def setPose( pos, vel, rot ):
    lib.setPose( pos, vel, rot )

#void setTilt( int iwing, double tilt ){
lib.setTilt.argtypes = [c_int,c_double]
lib.setTilt.restype  = None
def setTilt( iwing, tilt ):
    lib.setTilt( iwing, tilt )

#void fly( int n, int nsub, double dt, double* pos_, double* vel_, double* rot_ )
lib.fly.argtypes = [c_int, c_int, c_double, array2d,array2d,array3d]
lib.fly.restype  = None
def fly( pos, vel, rot, nsub=10, dt=0.01 ):
    n = len(pos)
    lib.fly( n, nsub, dt, pos, vel, rot )

#void flyAndShootTargets( int nsub, double dt, double* controlBuff, double* stateBuff, int* targetShot ){
lib.flyAndShootTargets.argtypes = [c_int, c_double, array1d, array1d, array2d]
lib.flyAndShootTargets.restype  = None
def flyAndShootTargets( controlBuff, stateBuff, targetHits, nsub=10, dt=0.01 ):
    lib.flyAndShootTargets( nsub, dt, controlBuff, stateBuff, targetHits )

#void setTargets( int nsub, double dt, double* controlBuff, double* stateBuff, int* targetShot ){
lib.setTargets.argtypes = [c_int, array2d ]
lib.setTargets.restype  = None
def setTargets( targets ):
    lib.setTargets( len(targets), targets )


#void setTargets( int nsub, double dt, double* controlBuff, double* stateBuff, int* targetShot ){
lib.getWorldPointer.argtypes = []
lib.getWorldPointer.restype  = c_void_p
def getWorldPointer():
    return lib.getWorldPointer()

# ===== view

class FlightView():
    
    def __init__(self, work_dir, wh=(800,600) ):
        #void init( int w, int h, void* craft1_, void* bursts_){
        
        self.libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
        #self.libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )
        self.libGL  = ctypes.CDLL( "/usr/lib/nvidia-375/libGL.so",   ctypes.RTLD_GLOBAL )
        
        #self.libSDL = ctypes.CDLL( "libSDL2.so", ctypes.RTLD_GLOBAL )
        #self.libGL  = ctypes.CDLL( "libGL.so",   ctypes.RTLD_GLOBAL )

        self.lib = initViewLib()
        print "========= libView ", libView
        self.lib.init.argtypes = [ c_int, c_int, c_void_p, c_char_p ]
        self.lib.init.restype  = None
        self.lib.init( wh[0], wh[1], lib.getWorldPointer(), work_dir )

        #void fly( int n, int nsub, double dt, double* pos_, double* vel_, double* rot_ )
        self.lib.draw.argtypes = []
        self.lib.draw.restype  = None

    def draw(self):
        self.lib.draw()


