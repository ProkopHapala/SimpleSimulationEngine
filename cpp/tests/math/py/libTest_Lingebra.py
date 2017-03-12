#!/usr/bin/python

import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p
import ctypes
import os

LIB_PATH       = os.path.dirname ( os.path.realpath(__file__) )
print "LIB_PATH : ",LIB_PATH
LIB_PATH_CPP   = os.path.normpath(LIB_PATH+'../../../../'+'/Build/tests/math')   # should be make it more robust in future ?
#LIB_PATH_CPP  = LIB_PATH

def recompile(path):
    print( path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )

recompile(LIB_PATH_CPP)
lib = ctypes.CDLL( LIB_PATH_CPP+"/libTest_Lingebra.so" )    

array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')    
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')

#    void   eig_Jacobi_init( int n, double* A, double* V, int* mjs   ){
lib.eig_Jacobi_init.argtypes   = [ c_int, array2d,array2d,array1i,array1i ]
lib.eig_Jacobi_init.restype    =  c_double

#    double eig_Jacobi_step( int n, double* A, double* V, int* mjs, int* ijmax, double vmax ){
lib.eig_Jacobi_step.argtypes   = [ c_int, array2d,array2d,array1i,array1i, c_double ]
lib.eig_Jacobi_step.restype    = c_double

# void eig_Jacobi( int n, double* A, double* V, double* es, double tol, int nMaxIter ){
lib.eig_Jacobi.argtypes   = [ c_int, array2d,array2d,array1d, c_double, c_int ]
lib.eig_Jacobi.restype    = c_int
def eig_Jacobi( A, V=None, es=None, tol=1e-9, nMaxIter=1000 ):
    n=len(A)
    if V is None:
        V=np.zeros((n,n))
    if es is None:
        es=np.zeros(n)
    lib.eig_Jacobi(n,A,V,es,tol,nMaxIter)
    return es,V
