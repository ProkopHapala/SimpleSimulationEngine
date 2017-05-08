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

def plot_depth( hitp, hitn ):
    #plt.imshow( hitp.x )
    plt.imshow( hitn[:,:,0:3] )
