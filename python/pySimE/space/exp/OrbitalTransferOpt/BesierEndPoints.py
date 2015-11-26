#!/usr/bin/env python

from sympy import *
from sympy.matrices import Matrix

'''

We have C2 continuous besier curve defined by

position y(t)
t  t2 t3 t4
-------------
Y0	+1 -3 +3 -1
Y1	+4    -6 +3
Y2	+1 +3 +3 -3
Y3	         +1 

first derivative y(t)
t  t2 t3 t4
-------------
Y0	-3  +6 -3
Y1	   -12 +9
Y2	+3  +6 -9
Y3	       +3
       
second derivative y(t)
t  t2 t3 t4
-------------
Y0	  6  -6
Y1	-12  18
Y2	 +6 -18
Y3	     +6

We want to find such controlpoints to set specific position and derivative at endpoints


'''

Y0, Y1, Y2, Y3, y0, y1, dy0, dy1, t = symbols('Y0 Y1 Y2 Y3 y0 y1 dy0 dy1 t')

eqs1 = (
     (Y0 + 4*Y1 + Y2)/6 - y0,
  (-3*Y0 +      3*Y2)/6 - dy0
)

out = solve( eqs1, Y0, Y1 )
print out

eqs2 = (
	     (    Y1 + 4*Y2  +   Y3)/6  -  y1   ,
           (- 3*Y1         + 3*Y3)/6  - dy1
)

out = solve( eqs2, Y2, Y3 )
print out