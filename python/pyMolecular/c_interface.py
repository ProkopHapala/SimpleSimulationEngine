import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float
import ctypes
import os

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/Molecular')

def recompile(path):
    print( path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )
    
# =========== main
recompile(LIB_PATH_CPP)

libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so",     ctypes.RTLD_GLOBAL )
libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",  ctypes.RTLD_GLOBAL )

lib = ctypes.CDLL( LIB_PATH_CPP+"/libMolecular.so" )

#lib.printHello()
lib.initWindow()



'''
print "KosmoSuiteCpp LIB_PATH     = ", LIB_PATH
print "KosmoSuiteCpp LIB_PATH_CPP = ", LIB_PATH_CPP 

lib    = ctypes.CDLL( LIB_PATH_CPP+"/lib"+name+ utils.ext )

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')


# ========= Nbody

# void set_Nbody( int n, double * mass, double * poss, double * vs ){
lib.nbody_setup.argtypes   = [c_int, array1d, array2d, array2d, array2d ]
lib.nbody_setup.restype    = None
def nbody_setup( mass, poss, vs, errs ):
	lib.nbody_setup( len(mass), mass, poss, vs, errs )
'''










