# -*- coding: utf-8 -*-
"""
http://adg.stanford.edu/aa241/drag/volumedrag.html
http://www.desktop.aero/appliedaero/compress3d/ssdragest.html
http://adg.stanford.edu/aa241/drag/SSDragCalc.html
http://aerorocket.com/HTV-3X.html
http://what-when-how.com/space-science-and-technology/rocket-propulsion-theory/
"""

from pylab import *;

mach = arange( 0,6, 0.01 );

#print mach

beta = sqrt((mach**2)-1.0)

cl  = 0.1
#ar  = 2.5
arl = 1
tc  = 1

'''
for arl in range(1,7):
	x   = pi*arl/4
	#cdwl = cl*cl*x/4*(sqrt(1+beta**2/x**2)-1) 
	#plot( mach, cdwl, label = 'Lift')
	cdwv = 4*(tc**2)*(beta**2+2*x**2)/(beta**2+x**2)**1.5
	plot( mach, cdwv, label = 'Volume')
'''

'''
c0 = 1

cdwv = 1/sqrt( mach**2 - 1.0 + 0.2 )
cdwv[mach<1]=0
plot( mach, c0 + cdwv, label = 'Volume')
'''

def dragFunc(  mach, cdf = 0.12, cds =0.05, cdw = 0.40, M0=0.70, Mtrans = 0.15  ):
	cd  = cdw*sqrt( 1/(mach**2 - M0**2) )
	dM = (mach -1 + Mtrans) /(2*Mtrans)
	w  =  3*dM**2 - 2*dM**3 
	cdtrans      = cds*(1-w) + w * cd
	indtrans     = (dM<1) & (dM>0) 
	cd[indtrans] = cdtrans[indtrans]
	cd[dM<0    ] = cds
	return cd + cdf
	
plot(mach, dragFunc(mach))

xticks(arange(0,6,0.5))
ylim(0,1)

grid()
#legend()

show()