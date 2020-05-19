
import numpy as np

'''

Product of Gaussians 1D

gij(x) =  exp( wi*(x-xi)^2 ) * exp( wj*(x-xj)^2 )

wi*( x^2 + 2*xi*x + xi^2 ) +  wj*( x^2 + 2*xj*x + xj^2 )
(wi+wj)*x^2  +   (wi*xi + wj*xj)*2*x   +   wi*x^2 + wj*x^2

wij*( x^2   +  2*x*xij +   xij^2 )   +    C

wij = (wi+wj)
wij*xij =   (wi*xi + wj*xj)

xij =    (wi*xi + wj*xj) / (wi + wj)

C =   wi*xi^2 + wj*xj^2     -    wij*xij^2
C =   wi*xi^2 + wj*xj^2     -   (wi*xi + wj*xj)^2 / (wi + wj)
C = wi*wj*(xi**2 - 2*xi*xj + xj**2)/(wi + wj)
C = wi*wj*(xi-xj)**2/(wi+wj)
'''

'''
# ======== Sympy Derivation
import sympy as sy
x,wi,xi,Ci,wj,xj,Cj = sy.symbols('x wi xi Ci wj xj Cj')
wij = (wi+wj)
xij = (wi*xi + wj*xj)/wij
C = wi*xi**2 + wj*xj**2   -  wij*xij**2

P1  = wi*(x-xi)**2 +  wj*(x-xj)**2 
P1cx = sy.collect(P1.expand(),x,evaluate=False)
print P1cx

P2 = wij*(x-xij)**2   +   C
P2=P2.expand()
P1=P1.expand()
Err = P1  - P2
print Err.simplify()
exit()
'''

'''
# ======== Sympy Check
import sympy as sy
x,wi,xi,Ci,wj,xj,Cj = sy.symbols('x wi xi Ci wj xj Cj')
wxi =  wi*xi
wxj =  wj*xj
W   = (wi+wj)
wx  = (wxi + wxj)
X   = wx/W
C   = wxi*xi + wxj*xj - wx*X

P1  = wi*(x-xi)**2 + wj*(x-xj)**2 
P2  = W*(x-X)**2   +   C
P2=P2.expand()
P1=P1.expand()
Err = P1  - P2
print Err.simplify()
exit()
'''

def gaussProduct1D( wi=1.0,xi=0.0,Ci=1, wj=1.0,xj=0.0,Cj=1 ):
    wxi  =  wi*xi
    wxj  =  wj*xj
    W    = (wi+wj)
    wx   = (wxi + wxj)
    X    = wx/W
    logC = wxi*xi + wxj*xj - wx*X
    C   = np.exp(-logC) * Ci * Cj
    return (W,X,C)

def evalGauss1D( xs, w=1.0, x0=0.0, C=1.0 ):
    return np.exp( -w*(xs-x0)**2 )*C

# ================== Main

if __name__ == "__main__":
    xs = np.linspace(-5,5,1000)
    import matplotlib.pyplot as plt
    
    g1 = ( 1.2,-0.8,1.0 )
    g2 = ( 0.6, 1.2,1.0 )

    g12 =  gaussProduct1D( wi=g1[0],xi=g1[1],Ci=g1[2],  wj=g2[0],xj=g2[1],Cj=g2[2] )

    y1     = evalGauss1D(xs,w=g1 [0],x0=g1 [1],C=g1 [2])
    y2     = evalGauss1D(xs,w=g2 [0],x0=g2 [1],C=g2 [2])
    y12    = evalGauss1D(xs,w=g12[0],x0=g12[1],C=g12[2])
    y12num = y1*y2

    plt.plot( xs, y1 , label='g1')
    plt.plot( xs, y2 , label='g1')
    plt.plot( xs, y12,    label='(g1*g2)ana')
    plt.plot( xs, y12num,':', label='(g1*g2)num')
    
    plt.legend()
    plt.grid()
    #plt.ylim(0,1.5)
    plt.show()


