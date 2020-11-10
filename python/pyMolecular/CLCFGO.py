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

c_int_p    = ctypes.POINTER(c_int)
c_double_p = ctypes.POINTER(c_double)


# ========= C functions

#void init( int natom_, int nOrb_, int perOrb_, int natypes_  ){
lib.init.argtypes = [ c_int, c_int, c_int, c_int ]
lib.init.restype  = None
def init( natom, nOrb, perOrb, natypes ):
    lib.init( natom, nOrb, perOrb, natypes )

#int* getIBuff(const char* name){ 
lib.getIBuff.argtypes = [c_char_p]
lib.getIBuff.restype  = c_int_p
def getIBuff(name,natom):
    ptr = lib.getIBuff(name)
    return np.ctypeslib.as_array( ptr, shape=(natom,))


#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p 
def getBuff(name,sh):
    ptr = lib.getBuff(name)
    #sh_ = (natom,)
    #if sh is not None:
    #    sh_ = sh_ + sh
    return np.ctypeslib.as_array( ptr, shape=sh)

#void testDerivs_Coulomb_model( int n, double x0, double dx ){
lib.testDerivsP_Coulomb_model.argtypes = [ c_int, c_double, c_double ]
lib.testDerivsP_Coulomb_model.restype  = c_double
def testDerivsP_Coulomb_model( n=100, x0=0.0, dx=0.1 ):
    return lib.testDerivsP_Coulomb_model( n, x0, dx )

#void testDerivs_Coulomb_model( int n, double x0, double dx ){
lib.testDerivsS_Coulomb_model.argtypes = [ c_int, c_double, c_double ]
lib.testDerivsS_Coulomb_model.restype  = c_double
def testDerivsS_Coulomb_model( n=100, x0=0.0, dx=0.1 ):
    return lib.testDerivsS_Coulomb_model( n, x0, dx )

#void testDerivsP_Total( int n, double x0, double dx ){
lib.testDerivsP_Total.argtypes = [ c_int, c_double, c_double ]
lib.testDerivsP_Total.restype  = c_double
def testDerivsP_Total( n=100, x0=0.0, dx=0.1 ):
    return lib.testDerivsP_Total( n, x0, dx )

#void testDerivsS_Total( int n, double x0, double dx ){
lib.testDerivsS_Total.argtypes = [ c_int, c_double, c_double ]
lib.testDerivsS_Total.restype  = c_double
def testDerivsS_Total( n=100, x0=0.0, dx=0.1 ):
    return lib.testDerivsS_Total( n, x0, dx )




#void testDerivsTotal( int n, double* xs, double* Es, double* Fs, int what ){
lib.testDerivsTotal.argtypes = [ c_int, array1d, array1d, array1d, c_int ]
lib.testDerivsTotal.restype  = c_double
def testDerivsTotal( xs, Es=None, Fs=None, what=0 ):
    n = len(xs)
    if Es is None: Es = np.zeros(n)
    if Fs is None: Fs = np.zeros(n)
    lib.testDerivsTotal( n, xs, Es, Fs, what )
    return Es,Fs

#void testDerivsCoulombModel( CLCFGO& solver, int n, double* xs, double* Es, double* Fs, int what  ){
lib.testDerivsCoulombModel.argtypes = [ c_int, array1d, array1d, array1d, c_int ]
lib.testDerivsCoulombModel.restype  = c_double
def testDerivsCoulombModel( xs, Es=None, Fs=None, what=0 ):
    n = len(xs)
    if Es is None: Es = np.zeros(n)
    if Fs is None: Fs = np.zeros(n)
    lib.testDerivsCoulombModel( n, xs, Es, Fs, what )
    return Es,Fs

# ========= Python Functions


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import CLCFGO_coulomb_derivs as ref

    init(2,2,2,1)
    ecoef = getBuff("ecoef",(2,2)  )
    esize = getBuff("esize",(2,2)  )
    epos  = getBuff("epos" ,(2,2,3))

    ecoefs = [[1.0,1.0],[1.0,1.0] ]
    esizes = [[1.0,1.0],[1.0,1.0] ]
    eXpos  = [[0.,+0.5],[-3.5,-1.5]]
    eYpos  = [[0.,+0.0],[ 0.0, 0.0]]
    eZpos  = [[0.,+0.0],[ 0.0, 0.0]]

    #ecoefs = [[+0.93,+0.68],[+0.65,+1.3]]
    #esizes = [[+1.30,+0.90],[+1.60,+0.7]]
    #eXpos  = [[+0.00,+0.50],[-3.50,-2.0]]
    #eYpos  = [[+0.00,+0.00],[+0.00,+0.0]]
    #eZpos  = [[+0.50,-0.30],[-0.40,+0.8]]

    '''
    dx = 0.025
    xs =  np.arange( 0.0, 3.0, dx )
    xs_ = 0.5*(xs[1:]+xs[:-1])

    # ============== Gaussian Overlap Product derivatives
    #esizes[0][0] = xs
    print "esizes", esizes

    C,s,p, dCr, dA,dB = ref.product3D_s_deriv( xs,eXpos[0][0], eXpos[0][1],eXpos[0][1] )
    (dSsi,dXsi,dXxi,dCsi) = dA

    plt.figure()
    plt.plot( xs, C,    label = "C" )
    plt.plot( xs, dCsi, label = "dC/dsa" )
    plt.plot( xs_, (C[1:]-C[:-1])/dx, label = "dC/dsa_num", lw=3,ls=":" )
    
    plt.plot( xs, p,    label = "p")
    plt.plot( xs, dXsi, label = "dp/dsa")
    plt.plot( xs_, (p[1:]-p[:-1])/dx, label = "dp/dsa_num", lw=3,ls=":" )
    
    plt.plot( xs, s,    label = "s")
    plt.plot( xs, dSsi, label = "ds/dsa")
    plt.plot( xs_, (s[1:]-s[:-1])/dx, label = "ds/dsa_num", lw=3,ls=":" )

    plt.legend()
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
    plt.title('Gaussian Overlap Derivative')

    # ============== Gaussian Electrostatics derivatives
    print "esizes", esizes

    plt.figure()
    E,fr,fs = ref.Coulomb( xs, 1.0 )
    plt.plot( xs, E,  label = "E(r)")
    plt.plot( xs, fr*xs, label = "dE/dr")
    plt.plot( xs_, (E[1:]-E[:-1])/dx, label = "dE/dr_num", lw=3,ls=":" )

    E,fr,fs = ref.Coulomb( 1.0, xs  )
    plt.plot( xs, E,  label = "E(s)")
    plt.plot( xs, fs*xs, label = "dE/ds")
    plt.plot( xs_, (E[1:]-E[:-1])/dx, label = "dE/ds_num", lw=3,ls=":" )
    
    plt.legend()
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')
    plt.title('Gaussian Coulomb Derivative')

    '''

    #plt.show()
    #exit()
    # =====================

    #ylims=[-5,5]
    #ylims=[-5,20]
    ylims=[-30,80]
    plt.figure(figsize=(14,8))
    
    # =========================================
    # ============== Derivs in Python =========
    # =========================================
    print " ========== Derivs in Python "

    x0 =  0.2
    dx =  0.025
    xs =  np.arange( x0, 3.0, dx )
    xs_ = (xs[1:]+xs[:-1])*0.5
    #eXpos[0][0] = xa 

    #(E, Fp,Fs) = ref.evalEFtot( ecoefs, esizes, eXpos, xa=xs ); F=Fp
    (E, Fp,Fs) = ref.evalEFtot( ecoefs, esizes, eXpos, sa=xs ); F=Fs
    #(E, F) = ref.evalEF_S_off ( xs, ecoefs, esizes, eXpos )

    plt.subplot(1,2,2)
    #plt.plot( xs, r   , label='r'   )
    #plt.plot( xs, r , label='r'  )
    #plt.plot( xs, Sab , label='Sab' )
    #plt.plot( xs, Qab , label='Qab'  )
    #plt.plot( xs, dQab, label='dSab_ana' )
    #plt.plot( xs_, (Sab[1:]-Sab[:-1])/dx,':', label='dSab_num' )
    #plt.figure(figsize=(12,10))
    plt.plot( xs, E,  label='E' )
    plt.plot( xs, F,  label='F_ana' )
    plt.plot( xs_,(E[1:]-E[:-1])/dx,':', label='F_num', lw=3 )
    #plt.plot( xs, fxi, label='fxi' )

    plt.title('Python')
    plt.legend()
    #plt.ylim(-30,40)
    #plt.ylim(-5,30)
    plt.ylim( ylims[0], ylims[1] )
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    # =========================================
    # ============== Derivs in C++ ============
    # =========================================


    plt.subplot(1,2,1)
    plt.title('C++')
    #plt.plot(l_xs,l_r,label="r")

    ecoef[:,:]   = np.array(ecoefs)[:,:]
    esize[:,:]   = np.array(esizes)[:,:]
    epos [:,:,0] = np.array(eXpos)[:,:]
    epos [:,:,1] = np.array(eYpos)[:,:]
    epos [:,:,2] = np.array(eZpos)[:,:]

    n = len(xs)
    #testDerivs_Coulomb_model( n=n, x0=0.0, dx=0.1 )    
    print "===>> RUN  C++ test : testDerivs_Total "
    #testDerivsP_Coulomb_model( n=n, x0=x0, dx=dx )
    #testDerivsS_Coulomb_model( n=n, x0=x0, dx=dx )
    #testDerivsP_Total        ( n=n, x0=x0, dx=dx )
    #testDerivsS_Total        ( n=n, x0=x0, dx=dx )

    Es,Fs = testDerivsTotal( xs, what=0 ) # position deriv
    #Es,Fs = testDerivsTotal( xs, what=1 ) # size     deriv
    print "===<< DONE C++ test : testDerivs_Total "

    '''
    l_xs     = getBuff( "l_xs",    (n,) )
    #l_r      = getBuff( "l_r",     (n,) )
    l_Q      = getBuff( "l_Q",     (n,) )
    #l_dQ_ana = getBuff( "l_dQ_ana",(n,) )
    #l_dQ_num = getBuff( "l_dQ_num",(n,) )
    l_E      = getBuff( "l_E",     (n,) )
    l_Fana   = getBuff( "l_Fana",  (n,) )
    l_Fnum   = getBuff( "l_Fnum",  (n,) )

    plt.plot(l_xs,l_E,label="E" )
    plt.plot(l_xs,l_Fana,label="Fana")
    plt.plot(l_xs,l_Fnum,label="Fnum",ls=':',lw=3)
    #plt.plot(l_xs,l_Q,label="Q")
    #plt.plot(l_xs,l_dQ_ana,label="dQ_ana")
    #plt.plot(l_xs,l_dQ_num,label="dQ_num", ls=':',lw=3)
    '''

    plt.plot(xs ,Es,label="E" )
    plt.plot(xs ,Fs,label="Fana")
    plt.plot(xs_,(Es[1:]-Es[:-1])/dx,label="Fnum",ls=':',lw=3)


    plt.legend()
    #plt.ylim(-30,40)
    #plt.ylim(-5,30)
    #plt.xlim(0,l_xs[-3])
    plt.ylim( ylims[0], ylims[1] ) 
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    #print "Fc++ %g Fpy %g " %(l_Fana[0],F[0])

    plt.show()