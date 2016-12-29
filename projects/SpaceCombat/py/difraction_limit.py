#!/usr/bin/python

import numpy as np
import matplotlib.pylab as plt

# ========== Functions

eV2nm = 1239.84193

wavelengths = [ ("CO2",10.6e-6), ("Nd:YAG",1024e-9), ("Ar8+",47e-9), ("1keV",eV2nm*1e-12), ("XFEL",0.05e-9) ]
distances   = [ (1e+6,10e+6,100e+6)]

# z * lambda =  A * d

def LaserRange( A, lamb=1e-6, d=1.0 ):
    return A*d/lamb
    
def LaserAperture( z, lamb=1e-6, A=1.0, d=1.0 ):
    return z*lamb/d

def main():
    zs = np.arange( 3.0, 9.0+0.001, 0.5 )
    zs = np.power( 10.0, zs ) 
    #Amin=+1e+100; Amax=-1e+100
    #print( zs )
    ax = plt.gca()
    Atxt = 0.013
    for i,(label,lamb) in enumerate(wavelengths):
        As = LaserAperture( zs, lamb=lamb )
        #print( As )
        lamb_ = lamb*1e+9
        label=label+(" %3.3gnm %3.3geV" %(lamb_,eV2nm/lamb_))
        plt.plot( zs/1000.0, As, label=label )
        #Amin = min( Amin, As.min() )
        #Amax = max( Amax, As.max() )
        ax.text( LaserRange(Atxt, lamb=lamb)/1000.0, Atxt, label, rotation=50, color='k', horizontalalignment='left', rotation_mode='anchor', verticalalignment='baseline')
    plt.xscale('log')
    plt.yscale('log')
    #print( zs.min(),zs.max(), Amax,Amin )
    #plt.xlim(zs.min(),zs.max())
    #plt.ylim(Amin,Amax)
    plt.ylim(0.01,100.0)
    plt.xlabel(r"Distance [ km ]")
    plt.ylabel(r"Aperture [ m  ]")
    #plt.legend(loc=2)
    plt.grid()
    #plt.savefig( "timeToDistance.png", bbox_inches='tight' )	
    plt.show()

if __name__ == "__main__":
    main()

