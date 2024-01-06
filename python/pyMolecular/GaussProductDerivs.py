
import sympy as sy




si, sj, s,  x, xi, xj, X = sy.symbols(  "si sj s  x xi xj X" )

gi = sy.exp( (x-xi)**2/(2*si**2))
gj = sy.exp( (x-xj)**2/(2*sj**2))

g = gi*gj

#lg = sy.simplify(sy.log(g))

lgi = -(x-xi)**2/(2*si**2)
lgj = -(x-xj)**2/(2*sj**2)

lg = lgi + lgj

#print sy.factor(sy.simplify( sy.expand(lg) ))
col = sy.collect(  sy.expand(lg), x, evaluate=False )

A = col[x**2]
B = col[x   ]
C = col[1   ]


lg = sy.collect(  sy.expand(lg), x )

print "A : ", A
print "B : ", B
print "C : ", C


X = sy.simplify( -B/A )
W = sy.simplify( A)

print "X ", sy.factor(sy.simplify(X))
print "W ", sy.factor(sy.simplify(W))

lg_ = sy.collect( W*((x-X)**2), x )

#print sy.simplify( sy.expand(lg_ - lg-C) )

print "lg_ ", lg_
print "lg  ", sy.simplify(lg - C )



