# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 21:53:33 2013

@author: asiJa
"""

# http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
# http://mathworld.wolfram.com/CubicSpline.html
# http://people.math.sfu.ca/~stockie/teaching/macm316/notes/splines.pdf


# clamped spline is defined by derivatives at start =  A and end = B 
def Coefs_clamped( y ,h, A, B ):
    n = len(ys)
    b[0]=
    d[0]=
    for i in range(n-1):
        b[i+1] = 2*( h[i] + h[i+1] )
        d[i+1] = 6*( (y[i+2]-y[i+1])/h[i+1] -  (y[i+1]-y[i])/h[i]  )
    TDMASolve( hs, b, hs, d):
    



# Tridiagonal Matrix Solver, 
# a - subdiagonal
# b - diagonal
# c - superdiagonal
# d - left
# note: function also modifies b[] and d[] params while solving
def TDMASolve(a, b, c, d):
    n = len(d) # n is the numbers of rows, a and c has length n-1
    for i in xrange(n-1):
        d[i+1] -= d[i] * a[i] / b[i]
        b[i+1] -= c[i] * a[i] / b[i]
    for i in reversed(xrange(n-1)):
        d[i] -= d[i+1] * c[i] / b[i+1]
    return [d[i] / b[i] for i in xrange(n)] # return the solution
