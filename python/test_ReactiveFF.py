#!/usr/bin/python

import sys
import os
import numpy as np

import pyMolecular.ReactiveFF  as rff
import pyMolecular.atomicUtils as au


#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import common    as PPU

def h2bonds( itypes, poss, hbonds, bsc=1.0 ):
    natom = len(poss)
    #ntot  = natom + natom*4
    xyzs    = np.zeros( (natom,5,3) )
    itypes_ = np.zeros( (natom,5), dtype=np.int32 )

    mask = itypes>0

    xyzs[:,0,:] = poss[:,:]
    xyzs[:,1,:] = poss[:,:] + hbonds[:,0,:]*bsc
    xyzs[:,2,:] = poss[:,:] + hbonds[:,1,:]*bsc
    xyzs[:,3,:] = poss[:,:] + hbonds[:,2,:]*bsc
    xyzs[:,4,:] = poss[:,:] + hbonds[:,3,:]*bsc
    itypes_[:,0 ] = 5 + itypes[:]
    itypes_[:   ,1:4] = 1
    itypes_[mask,4  ] = 1
    return xyzs.reshape(-1,3), itypes_.reshape(-1)

if __name__ == "__main__":
    print " ================ START "
    print " ================ START "

    #os.chdir( "/u/25/prokoph1/unix/git/SimpleSimulationEngine/cpp/Build/apps/MolecularEditor2" )

    c6    = -15.0
    R2vdW = 8.0
    rff.insertAtomType( 3, 1, 0.65,  1.0, -0.7, c6, R2vdW )
    rff.insertAtomType( 4, 2, 0.8 , 1.0, -0.7, c6, R2vdW )

    natom = 50
    rff.ralloc(natom)

    types  = rff.getTypes(natom)
    poss   = rff.getPoss(natom)
    qrots  = rff.getQrots(natom)
    hbonds = rff.getHbonds(natom)

    #itypes  = np.random.randint( 2, size=natom, dtype=np.int32 ); print "itypes", itypes
    itypes  = (np.random.rand( natom )*2.0 ).astype(np.int32); print "itypes", itypes

    rff.setTypes( natom, itypes )

    poss [:,:]  = ( np.random.rand(natom,3)-0.5 ) * 20.0
    qrots[:,:]  = np.random.rand(natom,4)-0.5
    rs          = np.sum(qrots**2, axis=1 )
    qrots      /= rs[:,None]

    #rff.relaxNsteps( nsteps=2000, F2conf=0.0, dt=0.05, damp=0.9 )

    fout = open( "rff_movie.xyz",'w')
    for itr in range(50):
        F2 = rff.relaxNsteps( nsteps=100, F2conf=0.0, dt=0.05, damp=0.9 )
        print ">> itr ", itr," F2 ", F2
        #au.writeToXYZ( fout, itypes, poss  )
        xyzs, itypes_ = h2bonds( itypes, poss, hbonds, bsc=1.0 )
        #print "itypes_,xyzs shapes : ", itypes_.shape,xyzs.shape
        au.writeToXYZ( fout, itypes_, xyzs  )
    fout.close()

    #print "rots_", rots_
    print ">>>> ALL DONE <<<<"
