#!/usr/bin/env python

# Computes the Clenshaw Curtis quadrature 
#http://www.scientificpython.net/1/post/2012/04/clenshaw-curtis-quadrature.html

from pylab import *
from scipy.fftpack import ifft

def ClenshawCurtisCoefs(n1):
  if n1 == 1:
    x = 0
    w = 2
  else:
    n = n1 - 1
    C = zeros((n1,2))
    k = 2*(1+arange(np.floor(n/2)))
    C[::2,0] = 2/hstack((1, 1-k*k))
    C[1,1] = -n
    V = vstack((C,flipud(C[1:n,:])))
    F = real(ifft(V, n=None, axis=0))
    x = F[0:n1,1]
    w = hstack((F[0,0],2*F[1:n,0],F[n,0]))
  return x,w

CCX=[[]]
CCW=[[]]
for i in range(1,15): 
	CCXs,CCWs = ClenshawCurtisCoefs(i)
	CCX.append(CCXs)
	CCW.append(CCWs)

'''
np.set_printoptions(precision=16)
print " xs = ",xs
print " ws = ",ws
'''



def f1(t): 	return (  (0.138-t)*(0.5454-t)*(0.8354-t) )**2.0
def f1(t): 	return (  t*(t-1)*(0.138-t)*(0.5454-t)*(0.8354-t) )**2.0

ts = arange(0,1,0.001);


def ClenshawCurtisQuadrature( a, b, f, degree=1):
	xs    = ( CCX[degree] + 1.0) / 2.0
	coefs =   CCW[degree]
	print len(xs)
	#print xs
	#print coefs 
	ys = f(xs)
	return 0.5*(b-a)*sum(ys*coefs)

for i in range(2,12+1):
	integral= ClenshawCurtisQuadrature( 0.0,1.0, f1, degree=i)
	print " CCQ[",i,"] Integral = ",integral

print "check = ", sum(f1(ts))/1000;


plot( ts, f1(ts) );

show()
