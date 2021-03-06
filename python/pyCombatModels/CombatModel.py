
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
"int newWeaponType( const char* str ){",
"int newUnitType( const char* str, int iprim, int isec ){",
"int int newUnit( int ityp, int n ){",
"void prepareCombat( int natt, int ndef, int* atts, int* defs ){",
"void solveCombat( int nrounds, int ndist, double* dists, double dt_max, double advance_dist, double maxCammo ){",
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

#libSDL = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libSDL2.so", ctypes.RTLD_GLOBAL )
#libGL  = ctypes.CDLL( "/usr/lib/x86_64-linux-gnu/libGL.so",   ctypes.RTLD_GLOBAL )


#cpp_name='CombatModels'
#cpp_utils.make(cpp_name)
#LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
#LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/'+cpp_name )
#lib = ctypes.CDLL( LIB_PATH_CPP+("/lib%s.so" %cpp_name) )

lib = cpp_utils.loadLib('CombatModels')

# ========= C functions

#  int newWeaponType( const char* str ){
lib.newWeaponType.argtypes  = [c_char_p] 
lib.newWeaponType.restype   =  c_int
def newWeaponType(str):
    #return lib.newWeaponType(_np_as(str,c_char_p)) 
    return lib.newWeaponType(str) 

#  int newUnitType( const char* str, int iprim, int isec ){
lib.newUnitType.argtypes  = [c_char_p, c_int, c_int] 
lib.newUnitType.restype   =  c_int
def newUnitType(str, iprim, isec):
    #return lib.newUnitType(_np_as(str,c_char_p), iprim, isec) 
    return lib.newUnitType(str, iprim, isec) 

#  int int newUnit( int ityp, int n ){
lib.newUnit.argtypes  = [c_int, c_int] 
lib.newUnit.restype   =  c_int
def newUnit(ityp, n):
    return lib.newUnit(ityp, n) 

#  void prepareCombat( int natt, int ndef, int* atts, int* defs ){
lib.prepareCombat.argtypes  = [c_int, c_int, c_int_p, c_int_p] 
lib.prepareCombat.restype   =  None
def prepareCombat(atts, defs):
    natt=len(atts); ndef=len(defs)
    atts=np.array(atts,dtype=np.int32)
    defs=np.array(defs,dtype=np.int32)
    return lib.prepareCombat(natt, ndef, _np_as(atts,c_int_p), _np_as(defs,c_int_p)) 

#  void solveCombat( int nrounds, int ndist, double* dists, double dt_max, double advance_dist, double maxCammo ){
lib.solveCombat.argtypes  = [c_int, c_int, c_double_p, c_double, c_double, c_double] 
lib.solveCombat.restype   =  None
def solveCombat(nrounds=1, dists=[20.], dt_max=0., advance_dist=0., maxCammo=1.):
    ndist=len(dists)
    dists=np.array(dists)
    return lib.solveCombat(nrounds, ndist, _np_as(dists,c_double_p), dt_max, advance_dist, maxCammo) 

# ========= MAIN

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    print( " ===== adding WeaponTypes" )
    kar98k = newWeaponType(" 400.   700. .1 .25 1.     2. 0. 0.  0. ") # rifle
    m67    = newWeaponType("  10.    30. .1 .25 5.    10. 1. 0.1 0.") # granade
    rpg17  = newWeaponType("  50.   250. .1 .25 10.   10. 3. .5 0. ")
    mg42   = newWeaponType(" 500.  1000. .1 .25 100.  60. 0. 0. 0. ") # machine gun
    gun37  = newWeaponType("1500.  5000. .2 .25 200.  10. 2. .1 0. ") # AT cannon light
    print( " ===== adding UnitTypes" )
    inf1   = newUnitType( "0   1.  0.1   1. 1000. 15.  2.", kar98k, m67   )
    infATG = newUnitType( "0   1.  0.1   1. 1000. 15.  2.", gun37, kar98k )
    infATR = newUnitType( "0   1.  0.1   1. 1000. 15.  2.", rpg17, kar98k )
    tank1  = newUnitType( "2  10. 11.0   4. 20000. 10. 6.", gun37, mg42   )
    print( " ===== adding Units" )
    unit_inf  = newUnit( inf1  ,1000  )
    unit_ATG  = newUnit( infATG,1000  )
    unit_ATR  = newUnit( infATR,1000  )
    unit_tank = newUnit( tank1 ,100   )
    print( " ===== to combat simulation" )
    print( " ===== assing units to combat " )
    prepareCombat( [unit_tank], [unit_inf,unit_ATR] )
    print( " ===== run combat " )
    exit(0)
    solveCombat()