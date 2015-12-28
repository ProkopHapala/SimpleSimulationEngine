#!/usr/bin/python

from pylab import *

'''

 ( (1-a)*yAL + a*yBL )^2   +     ( (1-a)*yAR + a*yBR )^2   minimize
 ( yAL + a*(yBL-yAL) )^2   +     ( yAR + a*(yBR-yAR) )^2   minimize
 ( yAL + a*yABL      )^2   +     ( yAR + a*yABR      )^2   minimize

  yAL*yAL + 2*a*yABL*yAL + a*a*yABL*yABL   +    yAR*yAR + 2*a*yABR*yAR + a*a*yABR*yABR

  2*yABL*   +   2*a*yABL*yAL            2*yABR*yAR   +   2*a*yABR*yABR = 0 

  a * ( yABL^2 + yABR^2  )   =   yABL*yAL + yABR*yABR

  a   =  ( yABL*yAL + yABR*yABR ) / ( yABL^2 + yABR^2  ) 

'''

[yAL,yBL,yAR,yBR] = random(4) * 100

yABL = yBL - yAL;
yABR = yBR - yAR;
a_best    =   -( yABL*yAL + yABR*yAR ) / ( yABL*yABL + yABR*yABR ) ; 

print " 1-a, a : ",  1-a_best, a_best

a = linspace( -10, 10.0, 1000 )

dens_tot = ((1-a)*yAL + a*yBL)**2  +   ((1-a)*yAR + a*yBR)**2

plot( a, dens_tot )
axvline( a_best,ls='--' )

show()
