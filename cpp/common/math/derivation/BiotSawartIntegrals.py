#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

#VortexQuantum = 1.0/(4*np.pi)
VortexQuantum = 1.0

# === Functions

def IwireElement(y,x):
    return VortexQuantum*   y/np.sqrt(x**2+y**2)**3

def IwireFinite(y,x1,x2):
    # finite wire with constant distribution of current
    # integral_a sin(a)*(cos(a)/sin(a)) = integral_a cos(a) = sin a
    c1 = x1/np.sqrt( x1**2+y**2 )
    c2 = x2/np.sqrt( x2**2+y**2 )
    return VortexQuantum*  (c2-c1)/y

def IwireSemiInfinite( y, x ):
    # semi infinite wire with constant distribution of current
    c   = x/np.sqrt( x**2+y**2 )         # cos(theta)
    out =  (1.0-c)/y
    #print "IwireSemiInfinite : ", x,y,c, out
    return VortexQuantum*out

def IwireSlope( y, x1, x2 ):
    # finite wire with linear distribution of current
    # integral_a sin(a)*(cos(a)/sin(a)) = integral_a cos(a) = sin a
    s1 = y/np.sqrt( x1**2+y**2 );
    s2 = y/np.sqrt( x2**2+y**2 );
    return VortexQuantum*  (s1-s2)/y

def IwireQuad( y, x1, x2 ):
    # finite wire with quadratic distribution of current
    # integral_a sin(a)*(cos(a)/sin(a))^2 = cos(a) + log( (1-cos(a))/(1+cos(a)) ) 
    # https://www.wolframalpha.com/input/?i=integrate+sin(a)*(cos(a)%2Fsin(a))%5E2+da
    c1  = x1/np.sqrt( x1**2+y**2 )
    c2  = x2/np.sqrt( x2**2+y**2 )
    I1  = c1 + np.log( np.sqrt( (1-c1)/(1+c1) ) )
    I2  = c2 + np.log( np.sqrt( (1-c2)/(1+c2) ) )
    return VortexQuantum*  (I1-I2)/y

def IwireCub( y, x1, x2 ): 
    # finite wire with cubic distribution of current
    # integral_a sin(a)*( cos(a)/sin(a) )^3 = (cos(2a)-3)/(2sin(a)) 
    # https://www.wolframalpha.com/input/?i=integrate+sin(a)*(cos(a)%2Fsin(a))%5E3+da
    r1  = np.sqrt( x1**2+y**2 )
    r2  = np.sqrt( x2**2+y**2 )
    c1  = x1/r1;  s1 = y/r1;
    c2  = x2/r2;  s2 = y/r2;
    I1  = 0.5*( c1**2 - s1**2 - 3 )/s1
    I2  = 0.5*( c2**2 - s2**2 - 3 )/s2
    return VortexQuantum*  (I1-I2)/y

def IsheetSemiInfinite( y, x1, x2 ):
    # sheet of semiinfinite vortex filements with constant density
    # Integral_y  (1-cos(a))/y =  Integral_y ( 1 - x/sqrt(x^2 + y^2))/y dy  = log(sqrt(x^+y^2)+x)
    #I1 = np.log( np.sqrt(x1**2+y**2)+y )
    #I2 = np.log( np.sqrt(x2**2+y**2)+y )
    r1  = np.sqrt(x1**2+y**2)
    r2  = np.sqrt(x2**2+y**2)
    out = np.log( (r2+y)/(r1+y) ) 
    #print  "IsheetSemiInfinite : ", y,x1,x2,r1,r2, out
    return  VortexQuantum*out

def IsheetSlopeSemiInfinite( y, x1, x2 ):
    # sheet of semiinfinite vortex filements with linearly increasing constant density
    # Integral_y  (1-cos(a)) =  Integral_y ( 1 - x/sqrt(x^2 + y^2))*y dy = y^2/2 - x*sqrt(x^2+y^2)
    r1  = np.sqrt(x1**2+y**2)
    r2  = np.sqrt(x2**2+y**2)
    out = (x2-x1) + y*np.log( (r1+x1)/(r2+x2) ) 
    return  VortexQuantum*out

def IsheetQuatSemiInfinite( y, x1, x2 ):
    # sheet of semiinfinite vortex filements with linearly increasing constant density
    # Integral_y  (1-cos(a)) =  Integral_y ( 1 - x/sqrt(x^2 + y^2))*y dy = y^2 - x sqrt(x^2+y^2)
    r1  = np.sqrt(x1**2+y**2)
    r2  = np.sqrt(x2**2+y**2)
    out = (x2*x2-x1*x1)*0.5 - y*(r2-r1)
    return  VortexQuantum*out

def integral_check_1D( f, F, xs  ):
    ys  = f(xs)
    Ys  = F(xs)
    Ys_ = np.cumsum(ys)*(xs[1]-xs[0])
    print xs,ys,Ys_
    return Ys,Ys_

# === setup

xmin=0.5
xmax=5.0
y0 = 1.0
dx  =0.001; hdx = dx*0.5;

# === main

#IwireSemiInfinite( 1.0, 1.0 ); exit()
#IwireSemiInfinite( 1.5, 0.5 ); exit()
#IsheetSemiInfinite( 1.0, 1.0, 5.0 ); exit()

xs = np.arange( xmin, xmax, dx ) + hdx

#f  = lambda x: IwireElement(y0,x)
#F  = lambda x: IwireFinite (y0,xmin,x+hdx)

#f  = lambda x: IwireElement(y0,x)*x
#F  = lambda x: IwireSlope  (y0,xmin,x+hdx)

#f  = lambda x: IwireElement(y0,x)*x*x
#F  = lambda x: IwireQuad  (y0,xmin,x+hdx)

#f  = lambda x: IwireElement(y0,x)*x*x*x
#F  = lambda x: IwireCub    (y0,xmin,x+hdx)

#f  = lambda x: IwireSemiInfinite (x,y0)
#F  = lambda x: IsheetSemiInfinite(y0,xmin,x+hdx)

#f  = lambda x: IwireSemiInfinite      (x,y0)*x
#F  = lambda x: IsheetSlopeSemiInfinite(y0,xmin,x+hdx)

f  = lambda x: IwireSemiInfinite      (x,y0)*x*x
F  = lambda x: IsheetQuatSemiInfinite(y0,xmin,x+hdx)

Ys,Ys_=integral_check_1D( f, F, xs )

plt.plot(xs+hdx,Ys , label='ana')
plt.plot(xs+hdx,Ys_, label='num')

plt.legend()
plt.show()


