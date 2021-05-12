
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
"void load_xyz( const char* fname ){",
"void init( int na, int ne ){",
"void eval(){",
"void info(){",
"double* getBuff(const char* name){",
"void setBuff(const char* name, double* buff){",
"int* getIBuff(const char* name){",
"void setIBuff(const char* name, int* buff){", 
"void setPauliModel(int i){",
"void setKPauli( double KPauli ){",
]
cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )


#cpp_name='CombatModels'
#cpp_utils.make(cpp_name)
#LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/'+cpp_name )
#lib = ctypes.CDLL( LIB_PATH_CPP+("/lib%s.so" %cpp_name) )

lib = cpp_utils.loadLib('eFF_lib')

# ========= C functions

# ========= Python Functions

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    plt.show()