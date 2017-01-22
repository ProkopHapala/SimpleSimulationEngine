#!/usr/bin/python

# compute distance by integration of tsielkovsky equation
# s2 = 0.5 a t^2 
# s1 = v t
# s2 = 0.5 a (s1/v)^2
# a  = 2 s2 / (s1/v)^2 
# v  = s1   / (2*s2/a)
# s2


import numpy as np
import matplotlib.pylab as plt

const_G = 9.81;

# ========= setup

dists     = 10.0**np.linspace( 4,9,40)
a_targets = 10.0**np.linspace(-3,1,5)
v_shots   = 10.0**np.linspace( 4,6,3)

clrs = ['b', 'c', 'm', 'y', 'k']

r_target_min = 1.0
r_target_max = 1000.0

# ========= Functions

def ManueverRad( dist, v_shot, a_target ):
    t = dist / v_shot
    return 0.5*a_target*(t**2)
    
def hitProb( dist, v_shot, a_target, r_target ):
    t      = dist / v_shot
    r_circ = 0.5*a_target*(t**2)
    return 1/(1 + (r_circ/r_target)**2)
    
def ManueverAcc( dist, r_target, v_shot ):
    t = dist / v_shot
    return 2.0*r_target/(t**2)    
    
def ManueverSpeed( dist, r_target, a_target ):
    t = np.sqrt( 2.0 * r_target / a_target )
    return dist / t    
    
def ManueverDist( r_target, a_target, v_shot ):    
    t = np.sqrt( 2.0 * r_target / a_target )
    return t * v_shot
    
def plotHitProb( r_target = 10.0, v_shot=100.0e+3 ):
    ax = plt.gca()
    for i,a in enumerate(a_targets):   
        #for v in v_shots: 
        #rs   = ManueverRad( dists, v, a )
        #plt.plot( dists*1e-3, rs, label='v=%f[km/s] a=%f[g]' %(v,a/const_G) )
        ps = hitProb( dists, v_shot, a, r_target )
        plt.plot( dists*1e-3, ps, label='v=%f[km/s] a=%f[g]' %(v_shot*1e-3,a/const_G), c=clrs[i] )
        ind = np.searchsorted( 1.0/ps, 1.0/1e-1 )
        ax.text( dists[ind]*1e-3, ps[ind], 'a=%1.2g [g]' %(a/const_G), rotation=-75, color=clrs[i], horizontalalignment='left', rotation_mode='anchor')
    plt.xscale('log'); plt.yscale('log');
    plt.ylim( 1e-6, 1.0 )
    plt.xlabel('distance [km]'); plt.ylabel('hit probablility [1]')
    plt.title('Hit probability to %3.0f [m] target by %3.0f [km/s] projectile' %(r_target,v_shot*1e-3) ); 
    plt.grid()
    
def plotAccs( ):
    for i,v in enumerate(v_shots): 
        plt.plot( dists*1e-3, ManueverAcc( dists, r_target_max, v )/const_G, c=clrs[i], lw=2, label='v=%f[km/s]' %(v*1e-3) )
        plt.plot( dists*1e-3, ManueverAcc( dists, r_target_min, v )/const_G, c=clrs[i], ls='--' )
    plt.xscale('log'); plt.yscale('log');
    plt.xlabel('dist [km]'); plt.ylabel('target acclearation [g]')
    plt.ylim( a_targets[0]/const_G, a_targets[-1]/const_G )
    #plt.title('Shooting on %f [m] target' %r_target ); 
    plt.legend(); plt.grid()
    
def plotSpeeds( ):
    ax = plt.gca()
    for i,a in enumerate(a_targets): 
        vs = ManueverSpeed( dists, r_target_max, a )
        plt.plot( dists*1e-3, vs*1e-3, c=clrs[i], lw=2, label='a=%1.4f[g]' %(a/const_G) )
        ax.text( ManueverDist( r_target_max, a, v_shots[0] )*1e-3, v_shots[0]*1e-3, 'a=%1.2g [g] r=%3.0f [m]' %(a/const_G,r_target_max), rotation=70, color=clrs[i], horizontalalignment='left', rotation_mode='anchor')
        plt.plot( dists*1e-3, ManueverSpeed( dists, r_target_min, a )*1e-3, c=clrs[i], ls='--' )
        #ax.text(txt_x, txt_y, txt, rotation=-40, color=clrs[i], horizontalalignment='left', rotation_mode='anchor', verticalalignment='top')
    plt.xscale('log'); plt.yscale('log');
    plt.xlabel('distance [km]'); plt.ylabel('projectile velocity [km/s]')
    plt.ylim( v_shots[0]*1e-3, v_shots[-1]*1e-3 )
    plt.title('Shooting on %3.0f [m] target ( -- %3.0f [m] )' %(r_target_max,r_target_min) ); 
    #plt.legend(); 
    plt.grid()

def main():
    #plt.figure(figsize=(12,6))
    #dists     = 10.0**np.linspace(3,9,7)
    #plt.subplot(1,2,1); plotAccs  (  )
    plt.subplot(1,2,1); plotSpeeds(  )
    plt.subplot(1,2,2); plotHitProb()
    plt.show()

if __name__ == "__main__":
    main()
    

