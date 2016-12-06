import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float
import ctypes
import os

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/Vis3D')

def recompile(path):
    print( path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )
    
# =========== main
recompile(LIB_PATH_CPP)

libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )

lib = ctypes.CDLL( LIB_PATH_CPP+"/libVis3D.so" )

#lib.printHello()
#lib.initWindow()

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

default_icolor = int("0xFF101010", 0)

# ========= C functions

#int makeSpheres( int n, double * poss_, double * colors_, double * radius )
lib.spheres.argtypes   = [c_int, array2d, array2d, array1d ]
lib.spheres.restype    = c_int
def spheres( poss, colors, radius ):
	return lib.spheres( len(poss), poss, colors, radius )

#int polyline( int n, double * points_, int closed )
lib.polyline.argtypes   = [c_int, array2d, c_int, c_int ]
lib.polyline.restype    = c_int
def polyline( poss, closed=0, icolor=default_icolor ):
	return lib.polyline( len(poss), poss, closed, icolor )

#int lines( int nedges, int * edges, double * points_ )
lib.lines.argtypes   = [c_int, array2i, array2d, c_int ]
lib.lines.restype    = c_int
def lines( edges, points, icolor=default_icolor ):
	return lib.lines( len(edges), edges, points, icolor )    
	
#int triangles( int ntris, int * tris, double * points_ )
lib.triangles.argtypes   = [c_int, array2i, array2d, c_int ]
lib.triangles.restype    = c_int
def triangles( tris, points, icolor=default_icolor ):
	return lib.triangles( len(tris), tris, points, icolor )
	
# void setGobVar( double f )
lib.setGobVar.argtypes   = [c_double]
lib.setGobVar.restype    = None
def setGobVar( f ):
	return lib.setGobVar( fr )





