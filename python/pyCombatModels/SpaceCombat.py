
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p
import ctypes
import os
import sys

#sys.path.append('../pyMeta')
import cpp_utils 
from cpp_utils import _np_as

cpp_utils.s_numpy_data_as_call = "_np_as(%s,%s)"

c_double_p = ctypes.POINTER(c_double)

# ===== To generate Interfaces automatically from headers call:
header_strings = [
    "void addTarget( char* str_target, char* str_shield ){",
    "void addGun( int n, char* str_gun, char* str_shot ){",
    "void evaluateCombat( int n, double* dists,  double* accels, double* out ){"
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_name='SpaceCombatLib'
cpp_utils.make(cpp_name)
lib    = ctypes.CDLL(  cpp_utils.BIN_PATH + "/lib" + cpp_name + cpp_utils.lib_ext )     # load dynamic librady object using ctypes 

# =============== C / Python interfaces 

#  void addTarget( const char* str_target, const char* str_shield ){
lib.addTarget.argtypes  = [c_char_p, c_char_p] 
lib.addTarget.restype   =  None
def addTarget(str_target, str_shield):
    return lib.addTarget(_np_as(str_target,c_char_p), _np_as(str_shield,c_char_p)) 

#  void addGun( int n, const char* str_gun, const char* str_shot ){
lib.addGun.argtypes  = [c_int, c_char_p, c_char_p] 
lib.addGun.restype   =  None
def addGun(n, str_gun, str_shot):
    return lib.addGun(n, _np_as(str_gun, c_char_p), _np_as(str_shot,c_char_p)) 

#  void evaluateCombat( int n, double* dists,  double* accels, double* out ){
lib.evaluateCombat.argtypes  = [c_int, c_double_p, c_double_p, c_double_p] 
lib.evaluateCombat.restype   =  None
def evaluateCombat( dists, accels, out ):
    n = len(dists)
    return lib.evaluateCombat(n, _np_as(dists,c_double_p), _np_as(accels,c_double_p), _np_as(out,c_double_p)) 

# =============== Test Run

if __name__ == "__main__":  
    import matplotlib.pyplot as plt

    nsamp  = 10
    ntg    = 1
    dists  = np.zeros((nsamp,))
    accels = np.zeros((nsamp,))
    out    = np.zeros((nsamp,ntg)) 

    addTarget(     "10.0  1.5e+8 2.5",        "2 10.0 0.5 150000" )
    addGun   ( 1,  "800 60000 1e+9 2e-4 10",  "0.15 0.12"         )
    evaluateCombat( dists, accels, out )
    
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
    plt.show()