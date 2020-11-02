import numpy as np
from   ctypes import c_int, c_double, c_bool, c_float, c_char_p, c_bool, c_void_p
import ctypes
import os

LIB_PATH      = os.path.dirname( os.path.realpath(__file__) )
LIB_PATH_CPP  = os.path.normpath(LIB_PATH+'../../../'+'/cpp/Build/libs/Molecular')

def recompile(path):
    print( "recompile path :", path )
    dir_bak = os.getcwd()
    os.chdir( path)
    os.system("make" )
    os.chdir( dir_bak )
    print( os.getcwd() )
    
# =========== main
recompile(LIB_PATH_CPP)

lib = ctypes.CDLL( LIB_PATH_CPP+"/libCLCFGO_lib.so" )

array1ui = np.ctypeslib.ndpointer(dtype=np.uint32, ndim=1, flags='CONTIGUOUS')
array1i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=1, flags='CONTIGUOUS')
array2i  = np.ctypeslib.ndpointer(dtype=np.int32,  ndim=2, flags='CONTIGUOUS')
array1d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags='CONTIGUOUS')
array2d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags='CONTIGUOUS')
array3d  = np.ctypeslib.ndpointer(dtype=np.double, ndim=3, flags='CONTIGUOUS')

# ========= C functions

#void init( int natom_, int nOrb_, int perOrb_, int natypes_  ){
lib.init.argtypes = [ c_int, c_int, c_int, c_int ]
lib.init.restype  = None
def init( natom, nOrb, perOrb, natypes ):
    lib.init( natom, nOrb, perOrb, natypes )

#int* getIBuff(const char* name){ 
lib.getIBuff.argtypes = [c_char_p]
lib.getIBuff.restype  = ctypes.POINTER(c_int)
def getIBuff(name,natom):
    ptr = lib.getIBuff(name)
    return np.ctypeslib.as_array( ptr, shape=(natom,))


#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = ctypes.POINTER(c_double)
def getBuff(name,sh):
    ptr = lib.getBuff(name)
    #sh_ = (natom,)
    #if sh is not None:
    #    sh_ = sh_ + sh
    return np.ctypeslib.as_array( ptr, shape=sh)

#void testDerivs_Coulomb_model( int n, double x0, double dx ){
lib.testDerivs_Coulomb_model.argtypes = [ c_int, c_double, c_double ]
lib.testDerivs_Coulomb_model.restype  = c_double
def testDerivs_Coulomb_model( n=100, x0=0.0, dx=0.1 ):
    return lib.testDerivs_Coulomb_model( n, x0, dx )

# ========= Python Functions


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    init(2,2,2,1)
    ecoef = getBuff("ecoef",(2,2)  )
    esize = getBuff("esize",(2,2)  )
    epos  = getBuff("epos" ,(2,2,3))

    n = 30
    testDerivs_Coulomb_model( n=n, x0=0.0, dx=0.1 )

    l_xs      = getBuff( "l_xs",     (n,) )
    l_Pi      = getBuff( "l_Pi",     (n,) )
    l_Qi      = getBuff( "l_Qi",     (n,) )
    l_dQi_ana = getBuff( "l_dQi_ana",(n,) )
    l_dQi_num = getBuff( "l_dQi_num",(n,) )
    l_E       = getBuff( "l_E",      (n,) )
    l_Fana    = getBuff( "l_Fana",   (n,) )
    l_Fnum    = getBuff( "l_Fnum",   (n,) )

    plt.figure(figsize=(14,8))
    plt.subplot(1,2,1)
    plt.plot(l_xs,l_Pi,label="Pi")
    plt.plot(l_xs,l_Qi,label="Qi")
    plt.plot(l_xs,l_dQi_ana,label="dQi_ana")
    plt.plot(l_xs,l_dQi_num,label="dQi_num", ls=':',lw=3)
    plt.plot(l_xs,l_E,label="E")
    plt.plot(l_xs,l_Fana,label="Fana")
    plt.plot(l_xs,l_Fnum,label="Fnum",ls=':',lw=3)

    plt.legend()
    plt.ylim(-30,40)
    plt.grid()


    # ==== Derivative of Coulomb term with considering the Charges
    import CLCFGO_coulomb_derivs as ref
    
    ca = 1.0
    cb = 1.0
    cc = 1.0
    cd = 1.0

    sa = 1.0
    sb = 1.0
    sc = 1.0
    sd = 1.0

    dx =  0.1
    xa =  np.arange( 0.01, 3.0, dx )
    xb =  0.0
    xc = -1.5
    xd =  0.0

    xs_ = (xa[1:]+xa[:-1])*0.5

    # overlaps
    Sab, si, xab, dQab, dA, dB = ref.product3D_s_deriv( sa,xa, sb,xb )
    Scd, sj, xcd, dQcd, dC, dD = ref.product3D_s_deriv( sc,xc, sd,xd )
    # coulomb
    s2        = si*si + sj*sj
    s         = np.sqrt(s2)
    r         = xab-xcd
    e, fx, fs = ref.Coulomb( r, s )
    dXxi      = dA[2] + xa*0 

    plt.subplot(1,2,2)
    plt.plot( xa, Sab , label='Sab' )
    plt.plot( xa, r   , label='r'   )
    #plt.plot( xa, dQab, label='dSab_ana' )
    #plt.plot( xs_, (Sab[1:]-Sab[:-1])/dx,':', label='dSab_num' )

    qij  = 4*Scd*Sab
    #qij  = Sab
    dQij = 4*Scd*dQab

    # Q: Why we dont need derivatives of charge ????
    #Fx = -fx*0.5*dA[1] # This works for zero initial distance between blobs
    Fx  = fx*r*dXxi
    Fpi = fx*r*qij # see 
    fxi = Fpi*dXxi
    
    print "Scd, 4*Scd ", Scd, 4*Scd
    print "For some reason each charge is scaled by 2.0"

    E = e*qij
    F = fxi + e*dQij   # total derivtive F = dE/dx = d(e*qi)/dx

    #plt.figure(figsize=(12,10))
    plt.plot( xa, E,  label='E' )
    plt.plot( xa, F,  label='dEdx_ana' )
    plt.plot( xs_, (E[1:]-E[:-1])/dx,':', label='dEdx_num', lw=3 )
    plt.plot( xa, fxi, label='fxi' )

    plt.legend()
    plt.ylim(-30,40)
    plt.grid()


    plt.show()