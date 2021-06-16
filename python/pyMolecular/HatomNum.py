
'''

E = <psi|   Nabla^2 + V |psi>

'''

const_K_eVA  =  7.61996440364  # [eV]
const_El_eVA = 14.3996453606  # [eV]

import numpy as np

def LaplaceSpherical( r, f ):
    '''
    https://en.wikipedia.org/wiki/Laplace_operator#Three_dimensions
    Lf = (1/r^2) * d_r{ r^2 * d_r{f} } = (1/r^2) ( dd_r{f} + 2r*d_r{f} ) 
    for function constant in phi and theta
    for n-dimensional spherically symmetric function : https://en.wikipedia.org/wiki/Laplace_operator#N_dimensions
    L{f} = d_r{d_r{f}} - (n-1)/r * d_r{f}
    '''
    dr    = r[1]-r[0]
    df_r  = (f[2:]-f[:-2])/(dr*2)
    ddf_r = (f[2:]+f[:-2]-2*f[1:-1])/(dr*dr)
    #return 2*r[1:-1]*df_r  + ddf_r
    r_ = r[1:-1]
    return (df_r*(2/r_)  + ddf_r)  # s=0.5

def Gauss(r,s=1.0):
    '''
    https://homepages.inf.ed.ac.uk/rbf/HIPR2/log.htm
    Wolpram Alpha
    Laplace{(e^(-(x^2 + y^2 + z^2)/(2 s^2))} = (e^(-(x^2 + y^2 + z^2)/(2 s^2)) (-3 s^2 + x^2 + y^2 + z^2))/s^4
    Laplace{(exp(-r^2/(2s^2))} = (exp(-r^2/(2s^2)) (-3 s^2 + r^2))/s^4
    '''
    s2=s*s
    r2=r*r
    g   = np.exp( -r2/(2*s2) )
    Lg  = g * (-3*s2 + r2)/(s2*s2)
    return g,Lg

def integrateRadial( r, f ):
    S = 4*np.pi*(r**2)
    return np.trapz(f*S, x=r)


def test_NumLaplace():
    g,Lg = Gauss(rs,s=1.0)
    LgNum = LaplaceSpherical( rs, g )
    plt.plot(rs,g,  label='g')
    plt.plot(rs,Lg, label='Lg')
    plt.plot(rs[1:-1],LgNum, label="LgNum")
    plt.legend()
    plt.grid()
    #plt.show()

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    roff = 1e-4;
    Rmax = 10.0
    dr   = 0.05
    rs = np.arange(roff,Rmax+roff,dr)

    V    = 1/rs 
    ss = np.arange(0.5,4.0,0.1)
    Es = np.zeros((3,len(ss)))
    for i,s in enumerate(ss):
        g,Lg = Gauss(rs,s=s)
        S    = integrateRadial( rs, g*g   )
        Ekin = integrateRadial( rs, g*Lg  )
        Eae  = integrateRadial( rs, g*V*g )

        Ekin *= (-const_K_eVA /(S))
        Eae  *= (-const_El_eVA/(S))
        Etot = Ekin+Eae

        Es[:,i] = (Etot,Ekin,Eae)

        print Ekin, Eae, Etot
    
    plt.plot(ss,Es[0,:],'k',label="Etot")
    plt.plot(ss,Es[1,:],'r',label="Ekin")
    plt.plot(ss,Es[2,:],'g',label="Eae")

    import eFF_terms as eff
    EkAn,EaeAn = eff.Hatom(ss); EaeAn*=np.sqrt(0.5); EtotAn = EkAn+EaeAn

    plt.plot(ss,EtotAn,"k:",lw=2,label="EtotAn")
    plt.plot(ss,EkAn  ,"r:",lw=2,label="EkinAn")
    plt.plot(ss,EaeAn ,"g:",lw=2,label="EaeAn")

    plt.grid()

    #plt.figure()
    #plt.plot(ss, Es[1,:],label="Ekin")
    #plt.plot(ss,-Es[2,:],label="Eae")
    #plt.yscale('log')

    plt.show()
