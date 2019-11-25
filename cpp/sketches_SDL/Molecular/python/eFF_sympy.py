#!/usr/bin/python

import sympy as sy


r, si, sj = sy.symbols( "r si sj" )

e1 = (2*si*sj/(si**2 + sj**2))**5
e2 = sy.exp( -2*r**2/(si**2+sj**2) )
e3 = (r/si)**2

E = e1*e2*e3


def makeDeriv( fx, dx, label=None ):
    dfx = sy.diff( fx, dx )
    de1si=sy.simplify(dfx) 
    if label is not None:
        #print label, dfx
        print label+"/e ", sy.simplify(dfx/fx) 
    return dfx

def error(a,b, label=None ):
    err = sy.simplify(a-b)
    if label is not None:
        print label, err
    return  err

de1si = makeDeriv( e1, si, "e1_si" )
de1sj = makeDeriv( e1, sj, "e1_sj" )
de1r  = makeDeriv( e1, r , "e1_r" )

de2si = makeDeriv( e2, si, "e2_si" )
de2sj = makeDeriv( e2, sj, "e2_sj" )
de2r  = makeDeriv( e2, r , "e2_r" )

de3si = makeDeriv( e3, si, "e3_si" )
de3sj = makeDeriv( e3, sj, "e3_sj" )
de3r  = makeDeriv( e3, r , "e3_r"  )

dEsi = sy.diff( E, si )
dEsj = sy.diff( E, sj )
dEr  = sy.diff( E, r  )

dEsi_ = e1*e2*de3si + e1*e3*de2si + e2*e3*de1si
dEsj_ = e1*e2*de3sj + e1*e3*de2sj + e2*e3*de1sj
dEr_  = e1*e2*de3r  + e1*e3*de2r  + e2*e3*de1r


error( dEsi, dEsi_, "dEsi " )
error( dEsj, dEsj_, "dEsj " )
error( dEr , dEr_ , "dEr  " )

