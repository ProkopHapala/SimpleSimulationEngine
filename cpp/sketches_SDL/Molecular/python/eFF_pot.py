#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

R2SAFE = 1.0e-8


def expQ( x, qq, w2, bExp, aExp ):
    r2     = x*x
    invrw2 = 1./( r2 + w2 );
    invrw  = np.sqrt(invrw2);
    E      =  qq*invrw;
    fr     = -qq*invrw2*invrw;
    if( bExp<0 ):
        r     = np.sqrt( r2+R2SAFE );
        Eexp  = aExp*np.exp( bExp*r );
        fr   += bExp*Eexp/r;
        E    += Eexp;
    f = x * fr
    return E, f

def combineAbW( abwi, abwj ):
    a =  abwi[0]*abwj[0];
    b = (abwi[1]+abwj[1])*0.5;
    w =  abwi[2]+abwj[2];
    return ( a, b, w )

default_eAbWs = [
[  0.0,  0.0, 0.00 ],  # Q = 0 //
[  0.0,  0.0, 0.7 ],  # Q = 1 // H
[ 10.0, -3.0, 0.7 ],  # Q = 2 // Be?
[ 10.0, -3.0, 0.7 ],  # Q = 3 // B
[ 10.0, -3.0, 0.7 ],  # Q = 4 // C
];

default_aAbWs = [
[  0.0,  0.0, 0.02 ],  # Q = 0 //
[  0.0,  0.0, 0.10 ],  # Q = 1 // H
[ 20.0, -5.0, 0.25 ],  # Q = 2 // Be?
[ 20.0, -5.0, 0.25 ],  # Q = 3 // B
[ 20.0, -5.0, 0.25 ],  # Q = 4 // C
];


default_eeAbw = [  2.0, -3.0, 2.0 ]
default_epAbw = [  0.0,  0.0, 2.0 ]

def eval_H2( x, y, aAbw=default_aAbWs[1], aeAbw=default_eAbWs[1], epAbw=default_epAbw, QQaa=+1., QQae=-1., QQee=+1. ):

    aaAbw = combineAbW( aAbw, aAbw )

    print "aaAbw, aeAbw, epAbw ", aaAbw, aeAbw, epAbw

    raa = x*2 + y*0
    rae = np.sqrt( x**2 + y**2 )
    ree = y*2 + x*0

    #print "raa.shape, rae.shape, ree.shape ", raa.shape, rae.shape, ree.shape

    Eaa,Faa  = expQ( raa, QQaa, aaAbw[2],  aaAbw[1],  aaAbw[0] )
    Eae,Fae  = expQ( rae, QQae, aeAbw[2],  aeAbw[1],  aeAbw[0] )
    Eee,Fee  = expQ( ree, QQee, epAbw[2],  epAbw[1],  epAbw[0] )

    return Eaa, Eae*2.0, Eee



'''

Note: It seems that H2 Molecule cannot be sable without varying Kinetic Energy

see:  
[1] http://aip.scitation.org/doi/10.1063/1.3272671
The dynamics of highly excited electronic systems: Applications of the electron force field
Julius T. Su, William A. Goddard

[2] http://dx.doi.org/10.1016/j.mechmat.2015.02.008
Non-adiabatic dynamics modeling framework for materials in extreme conditions
Hai Xiao, Andrés Jaramillo-Botero, Patrick L. Theofanis, William A. Goddard,
'''


if __name__ == "__main__":
    xs = np.arange( 0.0, 6.0, 0.05,  )

    #bEE     = -1.0;
    #aEE     =  2.0;
    #bEEpair = -1.0;
    #aEEpair =  0.1;

    ielem = 1;
    QQae = -1.0;
    QQaa = +1.0;
    QQee =  1.0;
    w2ee = (1.0)**2;
    eeAbw = default_eeAbw 
    epAbw = default_epAbw 
    eAbw  = default_eAbWs[ielem]
    aAbw  = combineAbW( default_aAbWs[ielem] , default_aAbWs[ielem] )

    print  "eAbw ", eAbw
    print  "aAbw ", aAbw

    Eee,Fee   = expQ( xs, QQee, eeAbw[2], eeAbw[1], eeAbw[0] )
    EeeP,FeeP = expQ( xs, QQee, epAbw[2], epAbw[1], epAbw[0] )
    Eae,Fae   = expQ( xs, QQae, eAbw[2], eAbw[1], eAbw[0] )
    Eaa,Faa   = expQ( xs, QQaa, aAbw[2], aAbw[1], aAbw[0] )

    #print "xs  ", xs
    #print "Eee ", Eee

    plt.plot(xs, Eee  ,label='Eee' )
    plt.plot(xs, EeeP ,label='EeeP')
    plt.plot(xs, Eae  ,label='Eae' )
    plt.plot(xs, Eaa  ,label='Eaa' )

    '''
    Eaa, Eae, Eee = eval_H2( xs, y=1.0 )
    plt.plot(xs, Eaa+Eae+Eee, 'k'  ,label='Eee' )
    plt.plot(xs, Eaa  ,label='Eaa' )
    plt.plot(xs, Eae  ,label='Eea ')
    plt.plot(xs, Eee  ,label='Eee' )
    '''

    '''
    eAbw      = default_eAbWs[1]
    Eae,Fae   = expQ( xs, -1.0, eAbw[2], eAbw[1], eAbw[0] )
    plt.plot(xs, Eae, label='H-e' )

    ys = [0.0,0.25,0.5,1.0,2.0,4.0]
    for y in ys:
        Eaa, Eae, Eee = eval_H2( xs, y=y )
        plt.plot(xs, Eaa+Eae+Eee, label='y %g' %y )
    '''

    plt.legend()
    plt.grid()
    plt.show()





