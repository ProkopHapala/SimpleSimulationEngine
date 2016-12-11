import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p
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

lib = ctypes.CDLL( LIB_PATH_CPP+"/libMolecular.so" )

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

default_icolor = int("0xFF101010", 0)

# ========= C functions

#void initWorld( char * workdir ){
lib.initWorld.argtypes   = [ c_char_p ]
lib.initWorld.restype    = None
def initWorld( workdir ):
	return lib.initWorld( workdir )

#double relax( int niter, double fmaxConverg ){
lib.relax.argtypes   = [ c_int, c_double ]
lib.relax.restype    = c_double
def relax( n, fconv=1e-5 ):
    return lib.relax( n, fconv )

#void exportAtoms( char * fname ){
lib.exportAtoms.argtypes   = [ c_char_p ]
lib.exportAtoms.restype    = None
def exportAtoms( fname ):
    return lib.exportAtoms( fname )
	

# ---- testing	

# void initComparator( int n, double * points ){
lib.initComparator.argtypes   = [ c_int, array2d ]
lib.initComparator.restype    = None
def initComparator( points ):
    n = len( points )
    return lib.initComparator( n, points )

# void compDistance( double * points ){
lib.compDistance.argtypes   = [ array2d ]
lib.compDistance.restype    = c_double
def compDistance( points ):
    return lib.compDistance( points )
	
# void initComparatorT( int n, double * points, int * types ){
lib.initComparatorT.argtypes   = [ c_int, array2d, array1i ]
lib.initComparatorT.restype    = None
def initComparatorT( points, types ):
    n = len( points )
    return lib.initComparatorT( n, points, types )

# double compDistanceT( int n, double * points, int * types ){
lib.compDistanceT.argtypes   = [ c_int, array2d, array1i ]
lib.compDistanceT.restype    = c_double
def compDistanceT( points, types ):
    n = len( points )
    return lib.compDistanceT( n, points, types )

#double getPlaneWaveDescriptor( double * center_, int np, double * points, int nk,  double * ks, double * coefs )
lib.getPlaneWaveDescriptor.argtypes   = [  array1d, c_int, array2d, c_int, array2d, array1d ]
lib.getPlaneWaveDescriptor.restype    = None
def getPlaneWaveDescriptor( points, ks, center=np.zeros(3) ):
    nk = len( ks )
    coefs = np.zeros( nk )
    lib.getPlaneWaveDescriptor( center, len(points), points, nk, ks, coefs )
    return coefs
	
#void testMultipole( int order, int np, double * ps, double * Qs,    int nsamples, double * psamples, double * Eref, double * Eaprox
lib.testMultipole.argtypes   = [ c_int, c_int, array2d, array1d,  c_int, array2d, array1d, array1d ]
lib.testMultipole.restype    = None
def testMultipole( ps, Qs, psamples, order=2 ):
    nsamples = len(psamples)
    Eref     = np.zeros(nsamples); 
    Eaprox   = np.zeros(nsamples);
    lib.testMultipole( order, len(ps), ps, Qs, nsamples, psamples, Eref, Eaprox )
    return Eaprox, Eref






