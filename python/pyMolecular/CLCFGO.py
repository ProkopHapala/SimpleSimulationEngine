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

#void loadFromFile( char const* filename, bool bCheck ){
lib.loadFromFile.argtypes = [ c_char_p, c_bool ]
lib.loadFromFile.restype  = c_bool
def loadFromFile( fname, bCheck=True ):
    return lib.loadFromFile( fname, bCheck )

#void init( int natom_, int nOrb_, int perOrb_, int natypes_  ){
lib.init.argtypes = [ c_int, c_int, c_int, c_int ]
lib.init.restype  = None
def init( natom, nOrb, perOrb, natypes ):
    lib.init( natom, nOrb, perOrb, natypes )

# void eval(){
lib.eval.argtypes = [ ]
lib.eval.restype  = c_double
def eval( ):
    return lib.eval( )

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
def getDimPointer( sh=(8,) ):
    # natypes   natom  nOrb  nBas  perOrb  perOrb2   nqOrb  nQtot 
    ptr = lib.getDimPointer()
    return  np.ctypeslib.as_array( ptr, shape=sh )

# void printAtomsAndElectrons(){
lib.printAtomsAndElectrons.argtypes = [ ]
lib.printAtomsAndElectrons.restype  = None
def printAtomsAndElectrons( ):
    lib.printAtomsAndElectrons( )

# void printSetup(){
lib.printSetup .argtypes = [ ]
lib.printSetup .restype  = None
def printSetup ( ):
    lib.printSetup ( )

#int* getIBuff(const char* name){ 
lib.getIBuff.argtypes = [c_char_p]
lib.getIBuff.restype  = c_int_p
def getIBuff(name,sh):
    if not isinstance(sh, tuple): sh=(sh,)
    ptr = lib.getIBuff(name)
    return np.ctypeslib.as_array( ptr, shape=sh)


#double* getBuff(const char* name){ 
lib.getBuff.argtypes = [c_char_p]
lib.getBuff.restype  = c_double_p 
def getBuff(name,sh):
    ptr = lib.getBuff(name)
    if not isinstance(sh, tuple): sh=(sh,)
    #sh_ = (natom,)
    #if sh is not None:
    #    sh_ = sh_ + sh
    #print "DEBUG type( ptr ) ", type( ptr ), sh
    return np.ctypeslib.as_array( ptr, shape=sh)

#void atomsPotAtPoints(         int n, double* ps, double* out, double s, double Q )
#void orbAtPoints     ( int io, int n, double* ps, double* out )
#void rhoAtPoints     ( int io, int n, double* ps, double* out )
#void hartreeAtPoints ( int io, int n, double* ps, double* out )

#void atomsPotAtPoints(         int n, double* ps, double* out, double s, double Q )
lib.atomsPotAtPoints.argtypes = [ c_int, array2d, array1d, c_double, c_double ]
lib.atomsPotAtPoints.restype  = c_double
def atomsPotAtPoints( ps, out=None, s=0.0, Q=1.0 ):
    n = len(ps)
    if out is None: out=np.zeros(n)
    lib.atomsPotAtPoints( n, ps, out, s, Q )
    return out

#void orbAtPoints     ( int io, int n, double* ps, double* out )
lib.orbAtPoints.argtypes = [ c_int, c_int, array2d, array1d ]
lib.orbAtPoints.restype  = c_double
def orbAtPoints( ps, io=0, out=None ):
    n = len(ps)
    if out is None: out=np.zeros(n)
    lib.orbAtPoints( io, n, ps, out )
    return out

#void rhoAtPoints     ( int io, int n, double* ps, double* out )
lib.rhoAtPoints.argtypes = [ c_int, c_int, array2d, array1d ]
lib.rhoAtPoints.restype  = c_double
def rhoAtPoints( ps, io=0, out=None ):
    n = len(ps)
    if out is None: out=np.zeros(n)
    lib.rhoAtPoints( io, n, ps, out )
    return out

#void hartreeAtPoints     ( int io, int n, double* ps, double* out )
lib.hartreeAtPoints.argtypes = [ c_int, c_int, array2d, array1d ]
lib.hartreeAtPoints.restype  = c_double
def hartreeAtPoints( ps, io=0, out=None ):
    n = len(ps)
    if out is None: out=np.zeros(n)
    lib.hartreeAtPoints( io, n, ps, out )
    return out

#double test_GaussIntegral_ST( int iMODE, int n, double sj, double* sis, double* rs, double* E, double* fr, double* fs ){
lib.test_GaussIntegral_ST.argtypes = [ c_int, c_int, c_double, array1d, array1d, array1d, array1d, array1d ]
lib.test_GaussIntegral_ST.restype  = c_double
def test_GaussIntegral_ST( iMODE=0, sj=1.0, sis=None, rs=None, r0=0.0, si=1.0 ):
    if rs  is None:
        n=len(sis)
        rs =np.zeros(n); rs[:]=r0
    if sis is None:
        n=len(rs) 
        sis=np.zeros(n); sis[:]=si
    E=np.zeros(n); fr=np.zeros(n); fs=np.zeros(n)
    lib.test_GaussIntegral_ST( iMODE, n, sj, sis, rs, E, fr, fs )
    return E, fr, fs

#test_Poisson( double Rmax, double gStep, double * line_rho=0, double* line_rho_=0, bool bPrint=0, bool bSave=0, useWf=true ){
lib.test_Poisson.argtypes = [ c_int, c_double, c_double, array1d, array1d, c_bool, c_bool, c_bool ]
lib.test_Poisson.restype  = c_double
def test_Poisson( io=0, Rmax=5.0, dx=0.1, bPrint=False, bSave=False, useWf=True, line_rho=None, line_rho_=None ):
    n = int(2*Rmax/dx)
    if line_rho  is None: line_rho =np.zeros(n)
    if line_rho_ is None: line_rho_=np.zeros(n)
    err2 = lib.test_Poisson( io, Rmax, dx, line_rho, line_rho_, bPrint, bSave, useWf )
    return err2, line_rho, line_rho_

#double test_CrossKinetic( int io, int jo, int nint, double dx, double Rmax, double gStep, double * line_Ek=0, double* line_Ek_g=0,  double * line_f1=0, double* line_f2=0, int bPrint=0, bool bSave=0 ){
lib.test_OrbInteraction.argtypes = [ c_int, c_int, c_int, c_int, c_double, c_double, c_double, array1d, array1d, array1d, array1d, c_int, c_bool ]
lib.test_OrbInteraction.restype  = c_double
def test_OrbInteraction( iMODE=1, io=0, jo=0, nint=40, dx=0.2, Rmax=5.0, gStep=0.1, bPrint=0, bSave=False, line_Ek=None, line_Ek_g=None, line_f1=None, line_f2=None ):
    '''
    iMODE :  1) Overlap S12 2) Kinetic T12 3) Colomb K12 
    '''
    ng = int( ( 2*Rmax + nint*dx )/gStep )
    #print(" test_OrbInteraction ng ", ng)
    if line_Ek   is None: line_Ek   =np.zeros(nint)
    if line_Ek_g is None: line_Ek_g =np.zeros(nint)
    if line_f1   is None: line_f1   =np.zeros(ng)
    if line_f2   is None: line_f2   =np.zeros(ng)
    err2 = lib.test_OrbInteraction( iMODE, io, jo, nint, dx, Rmax, gStep, line_Ek, line_Ek_g, line_f1, line_f2, bPrint, bSave )
    #print "line_Ek_g ", line_Ek_g
    #print "line_Ek   ", line_Ek
    return err2, line_Ek, line_Ek_g, line_f1, line_f2

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

#void setSwitches(bool bNormalize, bool bEvalKinetic, bool bEvalCoulomb, bool  bEvalExchange, bool  bEvalPauli, int iPauliModel,  bool bEvalAA, bool  bEvalAE, bool  bEvalAECoulomb, bool  bEvalAEPauli ){
lib.setSwitches.argtypes = [ c_bool, c_bool, c_bool, c_bool, c_bool, c_int, c_bool, c_bool, c_bool, c_bool ]
lib.setSwitches.restype  = None
def setSwitches( normalize=True, kinetic=True, coulomb=True, exchange=True, pauli=True, pauliModel=0, AA=True, AE=True, AECoulomb=True, AEPauli=True ):
    lib.setSwitches( normalize, kinetic, coulomb, exchange, pauli, pauliModel, AA, AE, AECoulomb, AEPauli )

#void setSwitches_(int bNormalize, int bEvalKinetic, int bEvalCoulomb, int  bEvalExchange, int  bEvalPauli, int bEvalAA, int bEvalAE, int bEvalAECoulomb, int bEvalAEPauli ){
lib.setSwitches_.argtypes = [ c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int ]
lib.setSwitches_.restype  = None
def setSwitches_( normalize=0, kinetic=0, coulomb=0, exchange=0, pauli=0, AA=0, AE=0, AECoulomb=0, AEPauli=0 ):
    lib.setSwitches_( normalize, kinetic, coulomb, exchange, pauli, AA, AE, AECoulomb, AEPauli )

#void setPauli( int iPauli ){ 
lib.setPauli.argtypes = [ c_int ]
lib.setPauli.restype  = None
def setPauli( iPauli ):
    lib.setPauli( iPauli )

# ========= Python Functions

def test_Gaussian_Overlap_Product_derivatives():
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

def test_Gaussian_Electrostatics_derivatives():
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

# =========================================
# ============== Derivs in Python ============
# =========================================
def plot_Derivs_Python():
    print " ========== Derivs in Python "
    #eXpos[0][0] = xa 
    (E, Fp,Fs) = ref.evalEFtot( ecoefs, esizes, eXpos, xa=xs ); F=Fp
    #(E, Fp,Fs) = ref.evalEFtot( ecoefs, esizes, eXpos, sa=xs ); F=Fs
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
def plot_Derivs_Cpp():

    init(natom,norb,perORb,1)  #  natom, nOrb, perOrb, natypes
    ecoef = getBuff("ecoef",(norb,perORb)  )
    esize = getBuff("esize",(norb,perORb)  )
    epos  = getBuff("epos" ,(norb,perORb,3))

    aQ     = getBuff("aQs",   (natom,)  )
    aQsize = getBuff("aQsize",(natom,)  )
    aPcoef = getBuff("aPcoef",(natom,)  )
    aPsize = getBuff("aPsize",(natom,)  )
    apos   = getBuff("apos"  ,(natom,3) )

    plt.subplot(1,2,1)
    plt.title('C++')
    #plt.plot(l_xs,l_r,label="r")

    ecoef[:,:]   = np.array(ecoefs)[:,:]
    esize[:,:]   = np.array(esizes)[:,:]
    epos [:,:,0] = np.array(eXpos)[:,:]
    epos [:,:,1] = np.array(eYpos)[:,:]
    epos [:,:,2] = np.array(eZpos)[:,:]

    aQ    [:]   = np.array(aQs    )
    aQsize[:]   = np.array(aQsizes)
    aPcoef[:]   = np.array(aPcoefs)
    aPsize[:]   = np.array(aPsizes)
    apos [:,:]= np.array(aposs)[:,:]
    #aposs_ = np.array(aposs)[:,:]
    #apos  [:,1] = aposs_[:,1]
    #apos  [:,2] = aposs_[:,2]

    n = len(xs)
    #testDerivs_Coulomb_model( n=n, x0=0.0, dx=0.1 )    
    print "===>> RUN  C++ test : testDerivs_Total "

    #setSwitches( normalize=False, kinetic=False, coulomb=False, exchange=False, pauli=False, pauliModel=1, AA=False, AE=True, AECoulomb=False, AEPauli=True );
    #setSwitches( normalize=False, kinetic=False, coulomb=False, exchange=False, pauli=False, pauliModel=1, AA=False, AE=True, AECoulomb=True, AEPauli=False );
    setSwitches( normalize=False, kinetic=False );

    Es,Fs = testDerivsTotal( xs, what=0 ) # position deriv
    #Es,Fs = testDerivsTotal( xs, what=1 ) # size     deriv
    print "===<< DONE C++ test : testDerivs_Total "

    plt.plot(xs ,Es,label="E" )
    plt.plot(xs ,-Fs,label="Fana")
    plt.plot(xs_,(Es[1:]-Es[:-1])/dx,label="Fnum",ls=':',lw=3)


if __name__ == "__main__":

    loadFromFile( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o.fgo" )
    #loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H2O_1g_8o.fgo" )
    printSetup()
    printAtomsAndElectrons()
    eval()
    exit()

    import matplotlib.pyplot as plt
    import CLCFGO_coulomb_derivs as ref

    natom   = 1
    aposs   = [[0.0,0.0,0.0],]
    aQs     = [4.0,]  # atomic nuclei charge (after screening core electrons)
    aPcoefs = [500.0,]  # atomic core pauli repulsion coeficient (strenght)
    aQsizes = [0.5,]  # atomic nuclei/pseudopotential size (radius of core electrons )
    aPsizes = [0.1,]  # atomic nuclei/pseudopotential size (radius of core electrons )

    norb   = 1
    perORb = 1
    ecoefs = [[1.0,], ]
    esizes = [[0.1,],]
    #eXpos  = [[0.,],]
    eXpos  = [[0.,],]
    eYpos  = [[0.,],]
    eZpos  = [[0.,],]

    '''
    norb   = 2
    perORb = 2
    ecoefs = [[1.0,1.0],[1.0,1.0] ]
    esizes = [[1.0,1.0],[1.0,1.0] ]
    #eXpos  = [[0.,+0.5],[-3.5,-1.5]]
    eXpos  = [[0.,+0.0],[ -0.5, 0.5]]
    eYpos  = [[0.,+0.0],[ 0.0, 0.0]]
    eZpos  = [[0.,+0.0],[ 0.0, 0.0]]
    '''

    '''
    norb   = 2
    perORb = 2
    #ecoefs = [[+0.93,+0.68],[+0.65,+1.3]]
    #esizes = [[+1.30,+0.90],[+1.60,+0.7]]
    #eXpos  = [[+0.00,+0.50],[-3.50,-2.0]]
    #eYpos  = [[+0.00,+0.00],[+0.00,+0.0]]
    #eZpos  = [[+0.50,-0.30],[-0.40,+0.8]]
    '''

    x0 =  -1.0
    dx =  0.05
    xs =  np.arange( x0, 4.0, dx )
    xs_ = (xs[1:]+xs[:-1])*0.5

    #test_Gaussian_Overlap_Product_derivatives()
    test_Gaussian_Electrostatics_derivatives()

    #plt.show()
    #exit()
    # =====================

    #ylims=[-5,5]
    #ylims=[-5,20]
    #ylims=[-30,80]
    ylims=[-200,400]
    plt.figure(figsize=(14,8))

    plot_Derivs_Python()
    plot_Derivs_Cpp()

    plt.legend()
    #plt.ylim(-30,40)
    #plt.ylim(-5,30)
    #plt.xlim(0,l_xs[-3])
    #plt.ylim( ylims[0], ylims[1] ) 
    plt.grid()
    plt.minorticks_on()
    plt.grid(which='minor', linestyle=':', linewidth='0.5', color='gray')

    #print "Fc++ %g Fpy %g " %(l_Fana[0],F[0])

    plt.show()