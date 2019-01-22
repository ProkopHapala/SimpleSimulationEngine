import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p
import ctypes
import os

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/Molecular')
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build-debug/libs/Molecular')

def recompile(path):
    print( "recompile path :", path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )
    
# =========== main
recompile(LIB_PATH_CPP)

lib = ctypes.CDLL( LIB_PATH_CPP+"/libReactiveFF.so" )

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

# ========= C functions

#void insertAtomType( int nbond, int ihyb, double rbond0, double aMorse, double bMorse, double c6, double R2vdW ){
lib.insertAtomType.argtypes = [ c_int, c_int, c_double, c_double, c_double, c_double, c_double ]
lib.insertAtomType.restype  = c_int
def insertAtomType( nbond, ihyb, rbond0, aMorse, bMorse, c6, R2vdW ):
    return lib.insertAtomType( nbond, ihyb, rbond0, aMorse, bMorse, c6, R2vdW )

#void ralloc(int natom){
lib.ralloc.argtypes = []
lib.ralloc.restype  = None
def ralloc(natom):
    lib.ralloc(natom)

#int*    getTypes(){ 
lib.getPoss.argtypes = []
lib.getPoss.restype  = ctypes.POINTER(c_int)
def getTypes(natom):
    return np.ctypeslib.as_array( lib.getTypes( ), shape=(natom))

#double* getPoss (){ 
lib.getPoss.argtypes = []
lib.getPoss.restype  = ctypes.POINTER(c_double)
def getPoss(natom):
    return np.ctypeslib.as_array( lib.getPoss( ), shape=(natom,3))

#double* getQrots(){ 
lib.getQrots.argtypes = []
lib.getQrots.restype  = ctypes.POINTER(c_double)
def getQrots(natom):
    return np.ctypeslib.as_array( lib.getQrots( ), shape=(natom,4))

#double* getHbonds(){ 
lib.getHbonds.argtypes = []
lib.getHbonds.restype  = ctypes.POINTER(c_double)
def getHbonds(natom):
    return np.ctypeslib.as_array( lib.getHbonds( ), shape=(natom,4,3) )

#double* getEbonds(){ 
lib.getEbonds.argtypes = []
lib.getEbonds.restype  = ctypes.POINTER(c_double)
def getEbonds(natom):
    return np.ctypeslib.as_array( lib.getEbonds( ), shape=(natom,4) )

#void setTypes( int natoms, int* types ){
lib.setTypes.argtypes = [c_int, array1i]
lib.setTypes.restype  = None
def setTypes(natom, itypes):
    lib.setTypes(natom, itypes)


#double relaxNsteps( int nsteps, double F2conf ){
lib.relaxNsteps.argtypes = [ c_int, c_double, c_double, c_double ]
lib.relaxNsteps.restype  = c_double
def relaxNsteps( nsteps=10, F2conf=0.0, dt=0.05, damp=0.9, ):
    return lib.relaxNsteps( nsteps, F2conf, dt, damp )

# ========= Python Functions





