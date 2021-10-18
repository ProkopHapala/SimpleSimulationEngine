
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p, c_char_p
import ctypes
import os
import sys

sys.path.append('../')
from pyMeta import cpp_utils 

c_double_p = ctypes.POINTER(c_double)
c_int_p    = ctypes.POINTER(c_int)

def _np_as(arr,atype):
    if arr is None:
        return None
    else: 
        return arr.ctypes.data_as(atype)

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

# ===== To generate Interfaces automatically from headers call:
header_strings = [
"void init_buffers(){",
"bool load_xyz( const char* fname ){",
"void init( int na, int ne ){",
"void eval(){",
"void info(){",
"double* getEnergyPointer(){",
"int*    getDimPointer   (){",
"double* getBuff(const char* name){",
"void setBuff(const char* name, double* buff){",
"int* getIBuff(const char* name){",
"void setIBuff(const char* name, int* buff){", 
"void setPauliModel(int i){",
"void setKPauli( double KPauli ){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )


#cpp_name='CombatModels'
#cpp_utils.make(cpp_name)
#LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/'+cpp_name )
#lib = ctypes.CDLL( LIB_PATH_CPP+("/lib%s.so" %cpp_name) )

cpp_utils.BUILD_PATH = os.path.normpath( cpp_utils.PACKAGE_PATH + '../../../cpp/Build/libs/Molecular' ) 
lib = cpp_utils.loadLib('eFF_lib')
array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')
# ========= C functions

#  void init_buffers(){
lib.init_buffers.argtypes  = [] 
lib.init_buffers.restype   =  None
def init_buffers():
    return lib.init_buffers() 

#  void load_xyz( const char* fname ){
lib.load_xyz.argtypes  = [c_char_p] 
lib.load_xyz.restype   =  c_bool
def load_xyz(fname):
    return lib.load_xyz(fname)
    #return lib.load_xyz(_np_as(fname,c_char_p)) 

#  void load_xyz( const char* fname ){
lib.load_fgo.argtypes  = [c_char_p] 
lib.load_fgo.restype   =  c_bool
def load_fgo(fname):
    return lib.load_fgo(fname);
    #return lib.load_xyz(_np_as(fname,c_char_p)) 

#  void init( int na, int ne ){
lib.init.argtypes  = [c_int, c_int] 
lib.init.restype   =  None
def init(na, ne):
    return lib.init(na, ne) 

#  void eval(){
lib.eval.argtypes  = [] 
lib.eval.restype   =  c_double
def eval():
    return lib.eval() 

#void evalFuncDerivs( int n, double* r, double* s, double* Es, double* Fs ){
lib.evalFuncDerivs.argtypes = [ c_int, array1d, array1d, array1d, array1d ]
lib.evalFuncDerivs.restype  = None
def evalFuncDerivs( r, s, Es=None, Fs=None ):
    r = r + s*0
    s = s + r*0
    n = len(r)
    if Es is None: Es=np.zeros(n)
    if Fs is None: Fs=np.zeros(n) 
    lib.evalFuncDerivs( n, r, s, Es, Fs )
    return Es,Fs

#  void info(){
lib.info.argtypes  = [] 
lib.info.restype   =  None
def info():
    return lib.info() 

#  double* getEnergyPointer(){
lib.getEnergyPointer.argtypes  = [] 
lib.getEnergyPointer.restype   = c_double_p
def getEnergyTerms( sh=(7,) ):
    # Ek=0, Eee EeePaul EeeExch Eae EaePaul Eaa
    ptr = lib.getEnergyPointer()
    return  np.ctypeslib.as_array( ptr, shape=sh )

#int*    getDimPointer   (){
lib.getDimPointer.argtypes  = [] 
lib.getDimPointer.restype   = c_int_p
def getDimPointer( sh=(3,) ):
    # ne=0 na=0 nDOFs=0
    ptr = lib.getDimPointer()
    return  np.ctypeslib.as_array( ptr, shape=sh )

#  double* getBuff(const char* name){
lib.getBuff.argtypes  = [c_char_p] 
lib.getBuff.restype   =  c_double_p
def getBuff( name, sh ):
    ptr = lib.getBuff(name)
    if not isinstance(sh, tuple): sh=(sh,)
    #sh_ = (natom,)
    #if sh is not None:
    #    sh_ = sh_ + sh
    print "DEBUG type( ptr ) ", type( ptr ), sh
    return np.ctypeslib.as_array( ptr, shape=sh)

#  void setBuff(const char* name, double* buff){
lib.setBuff.argtypes  = [c_char_p, c_double_p] 
lib.setBuff.restype   =  None
def setBuff(name, buff):
    return lib.setBuff(name, _np_as(buff,c_double_p)) 
    #return lib.setBuff(_np_as(name,c_char_p), _np_as(buff,c_double_p)) 

#  int* getIBuff(const char* name){
lib.getIBuff.argtypes  = [c_char_p] 
lib.getIBuff.restype   =  c_int_p
def getIBuff(name,sh):
    ptr = lib.getIBuff(name)
    if not isinstance(sh, tuple): sh=(sh,)
    return np.ctypeslib.as_array( ptr, shape=sh)
    #return lib.getIBuff(_np_as(name,c_char_p)) 

#  void setIBuff(const char* name, int* buff){
lib.setIBuff.argtypes  = [c_char_p, c_int_p] 
lib.setIBuff.restype   =  None
def setIBuff(name, buff):
    return lib.setIBuff(name, _np_as(buff,c_int_p)) 
    #return lib.setIBuff(_np_as(name,c_char_p), _np_as(buff,c_int_p)) 

#  void setPauliModel(int i){
lib.setPauliModel.argtypes  = [c_int] 
lib.setPauliModel.restype   =  None
def setPauliModel(i):
    return lib.setPauliModel(i) 

#  void setKPauli( double KPauli ){
lib.setKPauli.argtypes  = [c_double] 
lib.setKPauli.restype   =  None
def setKPauli(KPauli):
    return lib.setKPauli(KPauli) 

# ========= Python Functions

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    load_xyz("../../cpp/sketches_SDL/Molecular/data/e2_eFF.xyz")
    #load_xyz("../../cpp/sketches_SDL/Molecular/data/H2O_eFF.xyz")
    info()
    eval()

    plt.show()