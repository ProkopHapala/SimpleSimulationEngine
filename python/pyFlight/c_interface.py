import numpy as np
from   ctypes import c_int, c_double, c_double, c_bool, c_float, c_char_p
import ctypes
import os

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/libFlight')

def recompile(path):
    print( path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )
    
# =========== main
recompile(LIB_PATH_CPP)

lib = ctypes.CDLL( LIB_PATH_CPP+"/libFlight.so" )


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


