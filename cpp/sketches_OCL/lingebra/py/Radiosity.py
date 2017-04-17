#!/usr/bin/python

'''
https://en.wikipedia.org/wiki/Radiosity_(radiometry)
for each surface

Q_i    ... heat stored in element
J_i    ... radiance of element (power emmited from the surface) (there may be more surfaces than elements - two surfaces per element)
P_i    ... external heat influx (e.g. reactor, or irradiation by sun)
A_i    ... area of surface
T_i    ... temperature at surface 1
e_i    ... emmisivity
F_{ij} ... voiew factor, radiation coupling between the two surfaces (depends on angle and distance between surfaces)
sigma  ... stefan boltzman constant


Pure absorption-emmision (black body)
0 = dQ_i/dt =   P_i    -        sigma*e_i*A_i*T_i**4    +     Sum_j{ sigma F_{ij} e_i e_j A_i A_j T_j**4 } 
#           external heating      termo-emmision                      absorbtion         

Pure scattering (problem we have 2N unknown J_i, and T_i^4 )
we need more conditition   
1) coupling between radiation over surfaces ( conservation of energy in radiation transfer)
2) stationarity conditions - that temperature of each element is not changing any more

J_i =  sigma*e_i*A_i*T_i**4 + (1-e_i)*A_i Sum_j{ sigma F_{ij} J_i } 
#       termo-emmision            diffuse reflection


0 = dQ_i/dt = P_i - sigma*(e1_i+e2_i)*A_i*T_i**4 + sigma*A_i*Sum_j{ e1_i F_{ij} + e2_i F_{ij} J_j } 


solution:

# termo-emmision  is  total emmision minus reflected inflix 
Mi = sigma*e_i*A_i*T_i**4  =  J_i - (1-e_i)*A_i Sum_j{ sigma F_{ij} J_i }                    ...     Mi = Cij Jj


Cij = -1    ...  sigma (1-e_i) A_i F_{ij} ...
Kij =  0    ...  sigma    e_i  A_i F_{ij} 


0 = Pi - Cij Ji + e_i * Ji




Martinovy rovnice:

Pin_i  = Sum_j ( t1_ij PoutL_j + t2_ij Pout1_j )
PoutL_i = (1 - AL_i ) Pin_i + AL_i.S_i.sigma.T_i^4
PoutR_i = (1 - AR_i ) Pin_i + AR_i.S_i.sigma.T_i^4
PoutL_i + PoutR_i = Pin_i + Bi

4N rovnic, pro 4N neznamych (Pin_i,PoutL_i,PoutR_i,PoutL_i,T_i)


Martinovy muj mod 5N rovnic for 5N neznamych:   (PinL_i,PinR_i,PoutL_i,PoutR_i,PoutL_i,T_i)

PinL_i  = Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j )
PinR_i  = Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j )
PoutL_i = (1 - AL_i ) PinL_i     + sigma.AL_i.S_i.T_i^4
PoutR_i = (1 - AR_i ) PinR_i     + sigma.AR_i.S_i.T_i^4
PoutL_i + PoutR_i = PinL_i + PinR_i + Bi


simplyfy:


PoutL_i = (1 - AL_i ) Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j )   + sigma.AL_i.S_i.T_i^4
PoutR_i = (1 - AR_i ) Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j )   + sigma.AR_i.S_i.T_i^4
PoutL_i + PoutR_i = Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) + Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) + Bi

(1 - AL_i -1 ) Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j )   + sigma.AL_i.S_i.T_i^4
(1 - AR_i -1 ) Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j )   + sigma.AR_i.S_i.T_i^4
= Bi

S_i.sigma.T_i^4  = ( PoutL_i - (1 - AL_i ) Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) ) /AL_i
                 = ( PoutR_i - (1 - AR_i ) Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) ) /AR_i

AL_i ( sigma.S_i.T_i^4  - Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) + 
AR_i ( sigma.S_i.T_i^4  - Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) 
= Bi



AL_i ( ( PoutL_i - (1 - AL_i ) Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) ) /AL_i  - Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) + 
AR_i ( ( PoutR_i - (1 - AR_i ) Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) ) /AR_i  - Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) 
= Bi

AL_i ( ( PoutL_i + ( AL_i - 1 ) Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) - AL_i Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) ) /AL_i  + 
AR_i ( ( PoutR_i + ( AR_i - 1 ) Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) - AR_i Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) ) /AR_i
= Bi

| AL_i ( ( PoutL_i - Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) )/AL_i   
| AR_i ( ( PoutR_i - Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) )/AR_i   
| = Bi    


Radiance equation for two-side elements
==========================================================================
|   ( ( PoutL_i - Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) )            |
|   ( ( PoutR_i - Sum_j ( tRL_ij PoutL_j + tRR_ij PoutR_j ) )            |
|   = Bi                                                                 |
==========================================================================

then we can obtain T^4 from equation 
T_i^4  = ( PoutL_i - (1 - AL_i ) Sum_j ( tLL_ij PoutL_j + tLR_ij PoutR_j ) ) / ( AL_i.S_i.sigma )

moje rovnice:

PM_i ( ALi + ARi ) = ALi.PLi + ARi.PRi + Bi
PR_i ( CRi + Sum_j tRij  + ARi ) = Sum_j ( tRij PR_j  ) + ARi
PL_i ( CLi + Sum_j tLij  + ALi ) = Sum_j ( tLij PL_j  ) + ALi

'''
