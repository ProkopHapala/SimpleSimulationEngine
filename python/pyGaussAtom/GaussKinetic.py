
import numpy as np

'''

'''


# ======== Sympy Derivation
import sympy as sy
x,y,z,wi,xi,wj,xj = sy.symbols('x y z wi xi wj xj')



print "sy.erf( 1000000) ", sy.erf( 1000000).evalf()
print "sy.erf(-1000000) ", sy.erf(-1000000).evalf()


def Laplace(f):
    return sy.diff( f, x,x) + sy.diff( f, y,y) + sy.diff( f, z,z)

f2  = sy.exp( -wj*( (x-xj)**2  + y**2 + z**2 ) )
f1  = sy.exp( -wi*(  x**2      + y**2 + z**2 ) ) ; print "f"  , f1
Lf1 = Laplace(f1).simplify()                     ; print "Lf ", Lf1
t   = (Lf1*f2).simplify()                        ; print "t", t
#T   = sy.integrate( t, x )                       ; print "T", t


# t = f2*Lf1 = (px+pr)*gx*gr
gr = sy.exp( -(wi+wj)*(y**2 + z**2)  )
gx = sy.exp( -wi*x**2 - wj*(x-xj)**2 )
px =  4*(wi**2)*(x**2)
pr =  2*wi*( 2*wi*y**2 + 2*wi*z**2 - 3 )

print "check decomp t = f2*Lf1 = (px+pr)*gx*gr : ", (t - (px+pr)*gx*gr).simplify()


# 2*wi*(2*wi*y**2 + 2*wi*z**2 + 2*wi*x**2 - 3)*exp(-wi*(y**2 + z**2 + x**2) - wj*((x-xj)**2 + y**2 + z**2) )
# x:   exp( -wi*x**2 - wj*(x-xj)**2 )                 exp( -(wi+wj)*(y**2 + z**2) )
#   2*wi*x**2      +   C;   
#                      C = ( 2*wi*(2*wi*y**2 + 2*wi*z**2 - 3 )



# integal{     exp( -wi*x**2 - wj*(x-xj)**2 )  }
#     = 0.5* sqrt(pi/(wi+wj)) * exp( -wi*wj*xj**2/(wi+wj) )   *    erf( (wi*x+wj*(x-s) )/(sqrt(wi+wj)) )
#     

gx_chk = 0.5* sy.sqrt(sy.pi/(wi+wj)) * sy.exp( -wi*wj*xj**2/(wi+wj) )  *  sy.erf( (wi*x+wj*(x-xj) )/(sy.sqrt(wi+wj)) )
gx_chk = sy.diff( gx_chk, x )
gx_chk.simplify()
print " gx_chk"     , gx_chk
print " gx_chk - gx", (gx_chk - gx).expand().simplify()

exit()



Igx  = sy.integrate( gx   , x ); print Igx
Ipgx = sy.integrate( px*gx, x ); print Ipgx

# ================== Main

if __name__ == "__main__":
    xs = np.linspace(-5,5,1000)
    import matplotlib.pyplot as plt

    #plt.show()


