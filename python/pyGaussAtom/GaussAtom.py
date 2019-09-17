
import numpy as np
import scipy.linalg as linalg

np.set_printoptions(linewidth=200)

#invSqrt2pi = np.sqrt(np.pi*2)

#invSqrtPiHalf = 1/np.sqrt(np.pi*0.5)

gauss_renorm0 = 1/( 4*np.pi * np.sqrt(np.pi*0.5) )

# =========================== Functions

def applyH( f, ddfR, ddfT, V, k_h2m=0.1, bDebug=False ):
    Tf = -( ddfR + 2*ddfT )*k_h2m
    Vf = f*V
    if bDebug:
        ff  = np.trapz(  f*f*S, r )
        fTf = np.trapz( Tf*f*S, r )
        fVf = np.trapz( Vf*f*S, r )
        print "<f|f>", ff ,"<f|T|f> : ", fTf, " <f|V|f> ", fVf, " E tot ", fTf + fVf
    return Vf + Tf

def Gauss( r, r2=None, s=1.0, pre=gauss_renorm0, bNumRenorm=True ):
    s2    = s*s
    b     = 1./(-2*s**2)
    #invS  = 1./s
    #invS2 = invS*invS
    #invS3 = invS*invS2
    if r2 is None:
        r2 = r**2
    g    = np.exp(b*r2)   #   *invS3*pre
    if bNumRenorm:
        rho = g*g
        S   = 4*np.pi*r2
        norm = np.sqrt( np.trapz( rho*S, r ) )
        print "norm ", norm 
        g/= norm

    dg   =  g*(     r       )*2*b
    ddgR =  g*( 2*b*r*r + 1 )*2*b
    ddgT =  g*(         + 1 )*2*b
    return g,dg,ddgR,ddgT

def makeBasis( r, sigmas=[0.2,0.5,0.9] ):
    r2=r**2
    basis = []
    for s in sigmas:
        basis.append( Gauss(r, r2, s=s ) )
    return basis

def Hbasis_1D(basis, V, k_h2m=0.1 ):
    Hchis = []
    for bas in basis:
        f    = bas[0]
        ddfR = bas[2]
        Hchis.append( applyH( f, ddfR, 0, V, k_h2m ) )
    return Hchis

def Hbasis_3D(basis, V, k_h2m=0.1, bDebug=False ):
    Hchis = []
    for bas in basis:
        f    = bas[0]
        ddfR = bas[2]
        ddfT = bas[3]
        Hchis.append( applyH( f, ddfR, ddfT, V, k_h2m, bDebug=bDebug ) )
    return Hchis

def numDeriv( r, f ):
    return (f[2:]-f[:-2])/(r[2]-r[0])

# =========================== Main

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    xmax = 5.0

    N    = 1000+1
    Rcut = 4.0
    Rmax = 10.0
    r    = np.linspace(0,Rmax,N)
    r2   = r**2
    S    = 4*np.pi*r2

    w    =  0.3
    w2 = w**2

    COULOMB_CONST = 14.399644

    V = -np.sqrt( COULOMB_CONST/(r2 + w2 ) )

    plt.figure()
    g,dg,ddg,ddgT = Gauss( r, r2=None, s=0.6 )
    plt.plot(r,g,'b')
    plt.plot(r,dg,'g')
    plt.plot(r,ddg,'r')
    rd  = r[1:-1]
    rdd = r[2:-2]
    dg_  = numDeriv( r,  g  );   plt.plot(rd ,dg_ , 'g:',lw=4)
    ddg_ = numDeriv( rd, dg_ );  plt.plot(rdd,ddg_, 'r:',lw=4)
    plt.plot(r,V,'k')
    plt.grid()
    plt.title( "check derivs of Gaussian (numeric/analytic) " )
    #plt.show()
    #exit()

    colors = ['r' ,'g','b','m','c','y']
    sigmas = [0.25,0.5,1.0,1.5,2.0]
    nbas   = len(sigmas)

    # ======== plot bais
    basis = makeBasis( r, sigmas=sigmas )
    plt.figure()
    for i,bas in enumerate(basis):
        #plt.plot(r,basf[0])
        name = "basis[%i]" %i
        f = bas[0]
        rho = f*f
        print name,".norm() : ", np.trapz( rho*S, r )
        c = colors[i]
        plt.plot(r,bas[0]  ,c=c,lw=2.,ls='-', label=name)
        plt.plot(r,bas[0]*S,c=c,lw=1., ls=':' )
        plt.plot(r,bas[2]*S,c=c,lw=1., ls='--' )
        plt.plot(r,bas[3]*S,c=c,lw=1., ls='-.' )
    plt.plot(r,V,'k')
    plt.xlim(0,xmax)
    plt.legend()
    plt.grid()
    plt.title( "Gaussian Basis" )
    #plt.show()
    #exit()

    b3D = True

    k_h2m=0.0500001
    if b3D:
        Hchis = Hbasis_3D( basis, V, k_h2m=k_h2m, bDebug=True )
    else:
        Hchis = Hbasis_1D( basis, V, k_h2m=k_h2m )
    Hmat  = np.zeros( (nbas,nbas) )
    Smat  = np.zeros( (nbas,nbas) )
    chis     = [ bas[0] for bas in basis ]
    chis_ddR = [ bas[2] for bas in basis ]
    chis_ddT = [ bas[3] for bas in basis ]
    for i in xrange(nbas):
        for j in xrange(nbas):
            #Hmat[i,j] =  np.trapz( Hchis[i] * basis[j][0]    , r )
            #Hmat[i,j] =  np.trapz( Hchis[i] * basis[j][0] * S, r )
            if b3D:
                Hmat[i,j]  =  np.trapz( chis[j] * Hchis[i] * S  , r )
                Smat[i,j]  =  np.trapz( chis[i] *  chis[j] * S  , r )
            else:
                Hmat[i,j]  =  np.trapz( chis[j] * Hchis[i], r )
                Smat[i,j]  =  np.trapz( chis[i] *  chis[j], r )
    print " Hmat \n", Hmat
    print " Smat \n", Smat

    # Generalized eigenproblem    Hmat*Cs = Es*S*Cs
    # Result should be B-orthogonal


    '''
    #Es,Cs = np.linalg.eig( Hmat, b=Smat )
    Es,Cs  = linalg.eig( Hmat, b=Smat )
    #Cs=Cs.T
    #Cs = np.dot( Smat, Cs )
    CC  = np.dot( Cs.T, Cs )
    print " CC \n", CC
    CSC = np.dot( Cs.T, np.dot( Smat, Cs ) )
    print " CSC \n", CSC
    '''

    # ================= Lowdin

    # 1) ------ Orthogonalize Basis Set
    Ses,SVs = np.linalg.eig( Smat )
    SVs=SVs.T
    print "eigval(S) ",   Ses
    print "SVs \n", SVs

    SVVS = np.dot(SVs.T,SVs)
    print "SVVS \n", SVVS

    sSes = 1.0/np.sqrt( Ses )
    D = np.diag(sSes)

    print "sqrt(eigval(S)) ", sSes
    for i,e in enumerate(sSes):
        SVs[i,:]*=e

    Uchis     = np.dot( SVs, chis     )
    Uchis_ddR = np.dot( SVs, chis_ddR )
    Uchis_ddT = np.dot( SVs, chis_ddT ) 
    print "Uchis.shape ", Uchis.shape

    xUUx = np.zeros((nbas,nbas))
    for i in xrange(nbas):
        for j in xrange(nbas):
            xUUx[i,j] = np.trapz( Uchis[i] * Uchis[j] * S, r )
    #xUUx = np.dot( Uchis, Uchis.T )
    print "xUUx \n", xUUx

    # ...... plot S eigstates
    plt.figure()
    for i in range(nbas):
        plt.plot( r, Uchis[i], label="Uchi[%i]" %i )
    plt.plot(r,V,'k')
    plt.xlim(0,xmax)
    plt.grid()
    plt.title( "S-eigenstates" )
    plt.legend()
    #plt.show()

    # 2) ------ Solve in Orthogonal Basis Set

    UHU = np.zeros((nbas,nbas))
    for i in xrange(nbas):
        HUchi_i   =  applyH( Uchis[i], Uchis_ddR[i], Uchis_ddT[i], V, k_h2m=k_h2m, bDebug=True  )
        for j in xrange(nbas):
            UHU[i,j]  =  np.trapz( HUchi_i * Uchis[j] * S, r )
    print "UHU\n", UHU

    Es,Clow = np.linalg.eig( UHU )
    print "Es ", Es
    print "Clow \n", Clow

    CUchis = np.dot(Clow.T, Uchis   )
    #print "CUchis \n", CUchis

    #UCCU = np.dot( CUchis, CUchis.T )
    UCCU = np.zeros((nbas,nbas))
    for i in xrange(nbas):
        for j in xrange(nbas):
            UCCU[i,j] = np.trapz( CUchis[i] * CUchis[j] * S, r )

    print "UCCU \n", UCCU

    # ...... plot H eigstates
    plt.figure()
    for i in range(nbas):
        plt.plot( r, CUchis[i], label="CUchi[%i]" %i )
    plt.plot(r,V,'k')
    plt.xlim(0,xmax)
    plt.grid()
    plt.title( "H-eigenstates" )
    plt.legend()
    plt.show()

    '''
    sSHSs = np.dot( SVs, np.dot(Hmat,SVs.T) )
    print "sSHSs \n",   sSHSs

    Es,Vs = np.linalg.eig ( sSHSs )
    print "eigenval(Hlow) ",   Es
    print "eigenvec(Hlow) \n", Vs

    VV = np.dot( Vs,Vs.T )
    print "<Psi|Psi>\n" , VV 

    Hdiag = np.dot( Vs, np.dot( sSHSs,Vs.T ) )

    print "Hdiag \n", Hdiag

    Psis     =  np.dot(Vs,chis    )
    Psis_ddR =  np.dot(Vs,chis_ddR)
    Psis_ddT =  np.dot(Vs,chis_ddT)
    SSmat  = np.zeros( (nbas,nbas) )
    for i in xrange(nbas):
        for j in xrange(nbas):
            SSmat[i,j]  =  np.trapz( Psis[i] * Psis[j] * S, r )
    print " SSmat \n", SSmat

    exit()

    for i in xrange(nbas):
        applyH( Psis[i], Psis_ddR[i], Psis_ddT[i], V, k_h2m=k_h2m, bDebug=True )

    for i in xrange(nbas):
        print "eig [%i] ei"%i , Es[i]," vi ",Cs[i]

    # ======== plot eigstates
    plt.figure()
    for i in range(nbas):
        Ci = Cs[i]
        plt.plot( r, np.dot(Ci,chis), label="Psi[%i]" %i )
    plt.plot(r,V,'k')
    plt.xlim(0,xmax)
    plt.grid()
    plt.legend()
    plt.show()
    '''




