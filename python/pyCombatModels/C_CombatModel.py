
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p
import ctypes
import os
import sys

sys.path.append('../pyMeta')
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
    "int insertAtomType( int nbond, int ihyb, double rbond0, double aMorse, double bMorse, double c6, double R2vdW, double Epz ){",
    "void reallocFF(int natom){",
    "int*    getTypes (){",
    "double* getDofs  (){",
    "double* getFDofs (){",
    "double* getEDofs (){",
    "void setupFF( int natom, int* types ){",
    "void setGridShape( int* n, double* cell ){",
    "void bindGrids( double* atomMap, double*  bondMap ){",
    "double setupOpt( double dt, double damp, double f_limit, double l_limit ){",
    "void setBox(double* pmin, double* pmax, double* k){",
    "double relaxNsteps( int nsteps, double Fconv, int ialg ){",
]
cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )

cpp_name='FARFF'
cpp_utils.make(cpp_name)
lib    = ctypes.CDLL(  cpp_utils.BIN_PATH + "/lib" + cpp_name + cpp_utils.lib_ext )     # load dynamic librady object using ctypes 

# ========= C functions

if __name__ == "__main__":
    import matplotlib.pyplot as plt