#!/usr/bin/python

import sympy as sy

#print sy.__atribute__

t, F, ve, m0, fmass = sy.symbols( 't F ve m0 fmass' )

# ============ Linear Solution

'''

(specific impulse)
K = mass_flow = dv/dt = F/ve
 
P = F * ve = (F/ve) * ve * ve

s = Int_t{ Int_t{ a(t)} }

a = F/m = F/(  m0 - K*t ) = F/(  m0 - (F/ve)*t ) = F*ve/( m0*ve - F*t )

a = (P/ve)/(  m0 -  (P/(ve*ve)) * t )



tend:
 fmass =  m0/(m0 - (F/ve)*t)
 fmass = 1/(1 - (F/(ve*m0))*t )
 fmass*(1 - (F/(ve*m0))*t ) = 1
 fmass - fmass*F/(ve*m0) * t = 1
 fmass-1 = fmass*F/(ve*m0) * t
 (fmass-1)*(ve*m0)/(fmass*F) = tend


 fmass =  (m0 - (F/ve)*t)/m0
 fmass = 1 - (F/(ve*m0))*t
 1-fmass = (F/(ve*m0))*t
 (1-fmass)*ve*m0/F = tend
 



v = 

( -t*ve*log(F*t - m0*ve) + ve*t + 
    m0*ve*ve*log(F*t - m0*ve)/F

( -t*ve*v + ve*t +  m0*ve*ve*v/F


'''

#tend = ((fmass-1)/fmass) * ((ve*m0)/F)
tend = (1-fmass)*ve*m0/F 
fmass_ = m0/(m0 - (F/ve)*tend)

print "tend ", tend 

print " fmass_ ", sy.simplify( fmass )

a = F/(m0 - (F/ve)*t)

print "a  ", a

v = sy.integrate( a, t ).simplify()

print "v ", v

s = sy.integrate( v, t ).simplify()

#s_err = (  s -  (-t*-v + ve*t +  m0*ve*-v/F) ).simplify()
s_err = ( s - ( ve*t + v*( t -  m0*ve/F ) ) ).simplify()

print " s_err ", s_err

print "s ", s

print sy.cse(s)

s_1 = s.subs( t, tend ).simplify()

print "  s[fmass] : ", s_1



# ============ mono-Exponential Solution





# ============ Multi Exponential-Solution
