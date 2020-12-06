
import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p
import ctypes
import os
import sys
#from path import path

work_dir    = os.path.dirname( os.path.realpath( __file__ ) )
python_path = os.path.normpath( work_dir  + '../../' );    print "python_path :  ", python_path
sys.path.append( python_path   );                          print "sys.path    :  ", sys.path

from   pyMeta import cpp_utils
from   pyMeta.cpp_utils import _np_as
#import cpp_utils 
#from cpp_utils import _np_as

c_double_p = ctypes.POINTER(c_double)

# ===== To generate Interfaces automatically from headers call:
header_strings = [
    "void clearTargets( )",
    "void clear( )",
    "void addTarget( char* str_target, char* str_shield ){",
    "void addGun( int n, char* str_gun, char* str_shot ){",
    "void evaluateCombat( int n, double* dists,  double* accels, double* out ){"
]
#cpp_utils.writeFuncInterfaces( header_strings );        exit()     #   uncomment this to re-generate C-python interfaces

cpp_name='SpaceCombatLib'
lib  = cpp_utils.loadLib( cpp_name )

# =============== C / Python interfaces 

#  void clearTargets( )
lib.clearTargets.argtypes  = [] 
lib.clearTargets.restype   =  None
def clearTargets():
    return lib.clearTargets() 

#  void clear( )
lib.clearGuns.argtypes  = [] 
lib.clearGuns.restype   =  None
def clearGuns():
    return lib.clearGuns() 

#  void addTarget( const char* str_target, const char* str_shield ){
lib.addTarget.argtypes  = [c_char_p, c_char_p] 
lib.addTarget.restype   =  None
def addTarget(str_target, str_shield):
    return lib.addTarget(_np_as(str_target,c_char_p), _np_as(str_shield,c_char_p)) 
    #return lib.addTarget(str_target.encode('utf-8'), str_shield.encode('utf-8') ) 

#  void addGun( int n, const char* str_gun, const char* str_shot ){
lib.addGun.argtypes  = [c_int, c_char_p, c_char_p, c_double] 
lib.addGun.restype   =  None
def addGun(n, str_gun, str_shot, burstTime=1.0):
    return lib.addGun(n, _np_as(str_gun, c_char_p), _np_as(str_shot,c_char_p), burstTime ) 

#  void evaluateCombat( int n, double* dists,  double* accels, double* out ){
lib.evaluateCombat.argtypes  = [c_int, c_double_p, c_double_p, c_double_p] 
lib.evaluateCombat.restype   =  None
def evaluateCombat( dists, accels, out ):
    n = len(dists)
    return lib.evaluateCombat(n, _np_as(dists,c_double_p), _np_as(accels,c_double_p), _np_as(out,c_double_p)) 

# =============== Test Run

if __name__ == "__main__":  
    import matplotlib.pyplot as plt

    nsamp  = 100
    ntg    = 1
    dists  = (10.0**np.linspace( 1., 5., nsamp ))*1e+3
    accels = np.zeros((nsamp,)) + 0.1
    out    = np.zeros((nsamp,ntg))
    out2   = np.zeros((nsamp,ntg))
    out12  = np.zeros((nsamp,ntg))

    addTarget(     "10.0  1.5e+8 2.5 5.0",        "2 10.0 0.5 150000" )
    #               length,  maxForce,   maxPower,  scatter,  fireRate       mass    caliber
    addGun   ( 50,  "20     50000       3e+9       5e-5     100",          "0.03    0.01", 1.0   )
    evaluateCombat( dists, accels, out )
    clearGuns()
    addGun   ( 1 ,  "1500    160000      10e+9       5e-5     10",           "0.15    0.12", 1.0   )
    evaluateCombat( dists, accels, out2,  )
    addGun   ( 50,  "20     50000       3e+9       5e-5     100",          "0.03    0.01", 1.0   )
    evaluateCombat( dists, accels, out12  )
    
    plt.plot(dists*1e-3, out *1e-6, label='railGun: 50x,100rps 20m   ' ); 
    plt.plot(dists*1e-3, out2*1e-6, label='railGun: 1x,10rps   1500m ' );
    plt.plot(dists*1e-3, out12*1e-6, label='both' );
    #plt.xlabel('distance [km]'); plt.ylabel('health [1]'); plt.xscale('log');
    plt.xlabel('distance [km]'); plt.ylabel('damage [MJ]'); plt.xscale('log'); plt.yscale('log'); plt.ylim(1e-2,1e+4)
    plt.legend()
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
    plt.show()