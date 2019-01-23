#!/usr/bin/python

import sys
import os
import numpy as np
import time

import pyMolecular.ReactiveFF  as rff
import pyMolecular.atomicUtils as au


#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#import common    as PPU

def h2bonds( itypes, poss, hbonds, bsc=1.1 ):
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

def ebond2caps( ebonds, Ecut=-0.1 ):
    caps = np.zeros(ebonds.shape,dtype=np.int32) - 1
    caps[ebonds>Ecut] = 1;
    return caps

def removeSaturatedBonds(caps, itypes, xyzs,  ):
    itypes = itypes.reshape(-1,5)
    xyzs   = xyzs  .reshape(-1,5,3)
    #print ebonds
    #mask  = ebonds > Ecut
    mask   = caps >= 0
    mask[ itypes[:,4]==0,3] = False
    #print mask
    xyzs_   = [ xyzs  [:,0,:], xyzs  [mask[:,0],1,:], xyzs  [mask[:,1],2,:], xyzs  [mask[:,2],3,:], xyzs  [mask[:,3],4,:] ]
    itypes_ = [ itypes[:,0  ], itypes[mask[:,0],1  ], itypes[mask[:,1],2  ], itypes[mask[:,2],3  ], itypes[mask[:,3],4  ] ]
    return np.concatenate(xyzs_), np.concatenate(itypes_)

if __name__ == "__main__":
    print " ================ START "
    print " ================ START "

    #os.chdir( "/u/25/prokoph1/unix/git/SimpleSimulationEngine/cpp/Build/apps/MolecularEditor2" )

    c6    = -15.0
    R2vdW = 8.0
    rff.insertAtomType( 3, 1, 0.65,  1.0, -0.7, c6, R2vdW, 0.2 )
    rff.insertAtomType( 4, 2, 0.8 , 1.0, -0.7, c6, R2vdW, 0.0 )


    '''
    natom = 2
    rff.ralloc(natom)
    types  = rff.getTypes(natom)
    poss   = rff.getPoss(natom)
    qrots  = rff.getQrots(natom)
    hbonds = rff.getHbonds(natom)
    ebonds = rff.getEbonds(natom)
    itypes  = np.zeros(natom).astype(np.int32); print "itypes", itypes
    rff.setTypes( natom, itypes )
    rff.setSurf(K=1.0, x0=0.0, h=np.array([0.0,0.0,1.0]) )

    poss [:,:]  = np.array( [[-1.0,0.0,0.0],[1.0,0.0,0.0]] )
    qrots[:,:]  = np.random.rand(natom,4)-0.5
    rs          = np.sum(qrots**2, axis=1 )
    qrots      /= rs[:,None]
    '''

    natom = 20
    rff.ralloc(natom)
    types  = rff.getTypes(natom)
    poss   = rff.getPoss(natom)
    qrots  = rff.getQrots(natom)
    hbonds = rff.getHbonds(natom)
    ebonds = rff.getEbonds(natom)
    caps   = rff.getBondCaps(natom)
    #itypes  = np.random.randint( 2, size=natom, dtype=np.int32 ); print "itypes", itypes
    itypes  = (np.random.rand( natom )*1.3 ).astype(np.int32); print "itypes", itypes
    rff.setTypes( natom, itypes )
    poss [:,:]  = ( np.random.rand(natom,3)-0.5 ) * 10.0
    poss [:,2]  = 0.15
    qrots[:,:]  = np.random.rand(natom,4)-0.5
    rs          = np.sum(qrots**2, axis=1 )
    qrots      /= rs[:,None]

    rff.setBox( p0=np.array([-5.0,-5.0,-1.0]), p1=np.array([5.0,5.0,1.0]), K=-1.0, fmax=1.0  )
    rff.setSurf(K=-0.2, x0=0.0, h=np.array([0.0,0.0,1.0]) )

    #rff.relaxNsteps( nsteps=2000, F2conf=0.0, dt=0.05, damp=0.9 )

    '''
    fout = open( "rff_movie.xyz",'w')
    for itr in range(50):
        F2 = rff.relaxNsteps( nsteps=50, F2conf=0.0, dt=0.15, damp=0.9 )
        print ">> itr ", itr," F2 ", F2
        #au.writeToXYZ( fout, itypes, poss  )
        xyzs, itypes_ = h2bonds( itypes, poss, hbonds, bsc=1.1 )
        #print "itypes_,xyzs shapes : ", itypes_.shape,xyzs.shape
        xyzs, itypes_ = removeSaturatedBonds(ebonds, itypes_, xyzs, Ecut=-0.1 )
        #print ebonds
        au.writeToXYZ( fout, itypes_, xyzs  )
    fout.close()
    t2 = time.clock();
    '''
    
    t1 = time.clock();
    fout = open( "rff_movie.xyz",'w')
    for itr in range(10):
        F2 = rff.relaxNsteps( nsteps=50, F2conf=0.0, dt=0.15, damp=0.9 )
        print ">> itr ", itr," F2 ", F2 #, caps
        xyzs, itypes_ = h2bonds( itypes, poss, hbonds, bsc=1.1 )
        xyzs, itypes_ = removeSaturatedBonds(caps, itypes_, xyzs )
        au.writeToXYZ( fout, itypes_, xyzs  )
    rff.passivateBonds( -0.1 );
    print "passivation ", caps
    for itr in range(30):
        F2 = rff.relaxNsteps( nsteps=50, F2conf=0.0, dt=0.05, damp=0.9 )
        print ">> itr ", itr," F2 ", F2 #, caps
        xyzs, itypes_ = h2bonds( itypes, poss, hbonds, bsc=1.1 )
        xyzs, itypes_ = removeSaturatedBonds(caps, itypes_, xyzs )
        au.writeToXYZ( fout, itypes_, xyzs  )
    fout.close()
    t2 = time.clock();
    print "Relaxation time ", t2-t1

    #print ebonds
    au.saveXYZ( itypes+5, poss, "rff_skelet.xyz" )

    #print "rots_", rots_
    print ">>>> ALL DONE <<<<"
