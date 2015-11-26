#!/usr/bin/python

from pylab import *

def centrifugal( r, v, mass=1 ):
	return mass*(v**2)/r

def xExamples( examples, x0=1.0, ls='--', clr='k' ):
	for example in examples:
		axhline( example[1], ls=ls, color=clr )
		text   ( x0 ,example[1], example[0]    )
		#text   ( x0 ,x0, 'Hey'   )

r = 10**linspace( 1, 7 )


speeds=[   1310,   4300, 7900 ]
labels=[ 'Moon', 'Mars', 'LEO' ]
for i,speed in enumerate(speeds):
	g=centrifugal( r,  	speed )
	plot( r/1000, g/9.81, label=labels[i] )
legend()


grid()
xscale('log'); xlabel( 'length [km]' );
yscale('log'); ylabel( 'accelartion [G]' );

accel_examples=[
('Earth', 1.0 ),
('Houman limit',18.0),
('ICBM missile',100.0),
('smart shell',15500.0 )
]

xExamples( accel_examples )


show()

