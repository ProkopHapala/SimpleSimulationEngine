#!/usr/bin/python

import os
import re
import numpy    as np 
import matplotlib.pyplot as plt 

def getDiffuse( hitn, light_dir ):
    return hitn[:,:,0]*light_dir[0] + hitn[:,:,1]*light_dir[1] + hitn[:,:,2]*light_dir[2]

def getSpecular( hitn, light_dir, rd, gloss=256.0, power=2 ):
    slr       = light_dir[None,None,:] - rd[:,:,0:3];  
    specular  = ( ( slr[:,:,0]*hitn[:,:,0] + slr[:,:,1]*hitn[:,:,1] +slr[:,:,2]*hitn[:,:,2] )
          /np.sqrt( slr[:,:,0]**2          + slr[:,:,1]**2         + slr[:,:,2]**2) )
    #return 1.0/(1.0+gloss*(1-specular)**power )
    return specular**16

def blur( F ):
    return ( F[1:-1,1:-1] + F[:-2,1:-1] + F[2:,1:-1] + F[1:-1,:-2] + F[1:-1,2:] )*(1.0/5.0)
    
def min_blur( F ):
    F_ = F[1:-1,1:-1]
    F_ = np.minimum( F[1:-1,1:-1], F_ )
    F_ = np.minimum( F[ :-2,1:-1], F_ )
    F_ = np.minimum( F[2:  ,1:-1], F_ )
    F_ = np.minimum( F[1:-1, :-2], F_ )
    return F_
    
def getAO( hitp ):
    ts0 = hitp[:,:,3]
    ts  = blur( ts0 ); ts  = blur( ts ); ts  = blur( ts ); ts  = blur( ts )
    #ts  = min_blur( ts0 ); ts  = min_blur( ts ); ts  = min_blur( ts ); ts  = min_blur( ts )
    ts -= ts0[4:-4,4:-4]
    ts  = np.minimum(ts,0)
    ts  = 1.0/(1.0+2000.0*ts*ts)-1.0 
    return ts   

def plot_depth( hitp, hitn ):
    #plt.imshow( hitp.x )
    plt.imshow( hitn[:,:,0:3] )
