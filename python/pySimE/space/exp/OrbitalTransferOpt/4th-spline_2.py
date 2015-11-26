#!/usr/bin/env python

from sympy import *
from sympy.matrices import Matrix

'''

We want a set of localized basis fuctions to describe trajectory x(t),  v(t)=x'(t), a(t)=x''(t)
We have 4 boundary conditions
x(0   )=x_start,  x'(0   )=v_start
x(tend)=x_end  ,  x'(tend)=v_end


In between these points we can choose arbitary two parameters of the three:
xi=x(ti), vi=x'(ti), ai=x''(ti) 
There are there posibilities how to do that

'''

x0, x1, v0, v1, a0, a1, t = symbols('x0, x1, v0, v1, a0, a1, t')



def evalMat( Mlist,  blist  ):
	M = Matrix(Mlist)
	b = Matrix(blist)
	c = M.LUsolve( b )
	print 
	print " ================= Case: ",b.transpose()
	print 
	print c
	return c


ux0 =[ 1, 0, 0, 0 ]
ux1 =[ 1, 1, 1, 1 ]
uv0 =[ 0, 1, 0, 0 ]
uv1 =[ 0, 1, 2, 3 ]
ua0 =[ 0, 0, 2, 0 ]
ua1 =[ 0, 0, 2, 6 ]

S='''
================= CASE I =====================
0      T1     t1  T2  t2     T3    t3  T4 tend
|-------------|-------|------------|------|
x0            x1      x2           x3     xend
x0            v1      v2           v3     vend
==============================================
'''
print S

evalMat( [ ux0, ux1,   uv0, uv1  ],  [x0,x1,v0,v1]  )


S='''
================= CASE II ====================
0      T1     t1  T2  t2     T3    t3  T4 tend
|-------------|-------|------------|------|
x0            x1      x2           x3     xend
v0            a1      a2           a3     vend
==============================================
'''
print S

evalMat( [ ux0, ux1,   ua0, ua1  ],  [x0,x1,a0,a1]  )
evalMat( [ ux0, ux1,   uv0, ua1  ],  [x0,x1,v0,a1]  )
evalMat( [ ux0, ux1,   ua0, uv1  ],  [x0,x1,a0,v1]  )

S='''
================= CASE III ===================
0      T1     t1  T2  t2     T3    t3  T4 tend
|-------------|-------|------------|------|
x0            v1      v2           v3     xend
v0            a1      a2           a3     vend
==============================================
'''
print S

print " case ", [v0,v1,a0,a1] , " Is singular"
#evalMat( [ uv0, uv1,   ua0, ua1  ],  [v0,v1,a0,a1]  )
evalMat( [ ux0, uv1,   uv0, ua1  ],  [x0,v1,v0,a1]  )
evalMat( [ uv0, ux1,   ua0, uv1  ],  [v0,x1,a0,v1]  )


