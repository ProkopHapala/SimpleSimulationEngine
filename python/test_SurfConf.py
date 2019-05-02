#!/usr/bin/python

import sys
import os
import numpy as np

import pyMolecular.RigidMol  as rmol
import pyMolecular.atomicUtils as au


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import common    as PPU

def combineGeoms(mol,surf):
    es   = mol[0] + surf[0]
    xyzs = np.hstack( [np.array( mol[1:4] ), np.array(  surf[1:4] )] ).transpose().copy()
    return es, xyzs

def convGeom(mol):
    es   = mol[0]
    xyzs = np.array( mol[1:4] ).transpose().copy()
    return es, xyzs

def assignREQ( es, atomTypes_dct ):
    n = len(es)
    REQs = np.zeros( (n,3) )
    for i,e in enumerate(es):
        t=atomTypes_dct[ e ]
        REQs[i,0] = t[0];
        REQs[i,1] = t[1];
    return REQs

def rotDistPivot( rot1, rot2, ipiv=2):
    #print rot1[:,ipiv], rot2[:,ipiv]
    d = rot1[:,ipiv]-rot2[:,ipiv]
    return np.sqrt(np.dot(d,d))
    #return 1-np.dot(rot1[:,ipiv],rot2[:,ipiv])

def getSimilarRotations( rots, rot0, dcut=0.1, ipiv=2):
    rots_ = []
    for rot in rots:
        dist = rotDistPivot( rot, rot0, ipiv=ipiv )
        if(dist<dcut):
            rots_.append(rot)
    return rots_

def findSimilarRot(rot,rot0s, dcut=0.1, ipiv=2):
    for i,rot0 in enumerate(rot0s):
        dist = rotDistPivot( rot, rot0, ipiv=ipiv )
        if(dist<dcut):
            return i
    return -1

def pruneRots(rots, dcut=0.1, ipiv=2):
    rots_ = []
    for i,rot in enumerate(rots):
        isim = findSimilarRot(rot,rots_,dcut=dcut, ipiv=ipiv)
        print "==== rot ", i, isim #, rot
        if isim < 0:
            rots_.append(rot)
    return rots_

def sphereTangentSpace(n=100):
    golden_angle = np.pi * ( 3.0 - np.sqrt(5.0) )
    theta  = golden_angle * np.arange(n)
    z      = np.linspace(1.0 - 1.0/n, 1.0/n - 1.0, n)
    radius = np.sqrt( 1.0 - z*z )
    cas  = np.cos(theta)
    sas  = np.sin(theta)
    rots = np.zeros( (n,3,3) )
    rots[:,2,0] = radius * cas
    rots[:,2,1] = radius * sas
    rots[:,2,2] = z
    rots[:,0,0] = -sas
    rots[:,0,1] =  cas
    rots[:,1,:] =  np.cross( rots[:,2,:], rots[:,0,:] )
    return rots

def drawTangetSpace(fws,ups,lfs, lenght=0.1):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    #print fws
    ax.quiver( fws[:,0], fws[:,1], fws[:,2], fws[:,0], fws[:,1], fws[:,2], length=lenght, color='r' )
    ax.quiver( fws[:,0], fws[:,1], fws[:,2], lfs[:,0], lfs[:,1], lfs[:,2], length=lenght, color='g' )
    ax.quiver( fws[:,0], fws[:,1], fws[:,2], ups[:,0], ups[:,1], ups[:,2], length=lenght, color='b' )
    ax.set_xlim3d(-1.0, 1.0)
    ax.set_ylim3d(-1.0, 1.0)
    ax.set_zlim3d(-1.0, 1.0)
    #plt.axis('equal')

def drawRots( rots, lenght=0.1):
    rots = np.array(rots)
    drawTangetSpace(rots[:,2,:],rots[:,1,:],rots[:,0,:], lenght=lenght)

def randomRotations( n=100 ):
    '''
     http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
     http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
     http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.53.1357&rep=rep1&type=pdf
       RAND _ ROTATION   Author: Jim Arvo, 1991
       This routine maps three values (x[0], x[1], x[2]) in the range [0,1]
       into a 3x3 rotation matrix, M.  Uniformly distributed random variables
       x0, x1, and x2 create uniformly distributed random rotation matrices.
       To create small uniformly distributed "perturbations", supply
       samples in the following ranges
           x[0] in [ 0, d ]
           x[1] in [ 0, 1 ]
           x[2] in [ 0, d ]
      where 0 < d < 1 controls the size of the perturbation.  Any of the
      random variables may be stratified (or "jittered") for a slightly more
      even distribution.
    '''
    vrand = np.random.rand(n,3)
    theta = vrand[:,0] * 2.0*np.pi
    phi   = vrand[:,1] * 2.0*np.pi
    z     = vrand[:,2] * 2.0;     
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    r  = np.sqrt( z );
    Vx = np.sin ( phi ) * r;
    Vy = np.cos ( phi ) * r;
    Vz = np.sqrt( 2.0 - z );
    # Compute the row vector S = Transpose(V) * R, where R is a simple
    # rotation by theta about the z-axis.  No need to compute Sz since
    # it's just Vz.
    st = np.sin( theta );
    ct = np.cos( theta );
    Sx = Vx * ct - Vy * st;
    Sy = Vx * st + Vy * ct;
    # Construct the rotation matrix  ( V Transpose(V) - I ) R, which
    # is equivalent to V S - R.
    rots = np.zeros( (n,3,3) )
    rots[:,0,0]= Vx * Sx - ct
    rots[:,0,1]= Vx * Sy - st
    rots[:,0,2]= Vx * Vz
    rots[:,1,0]= Vy * Sx + st
    rots[:,1,1]= Vy * Sy - ct
    rots[:,1,2]= Vy * Vz
    rots[:,2,0]= Vz * Sx
    rots[:,2,1]= Vz * Sy
    rots[:,2,2]= 1.0 - z
    return rots
    '''
    return [
    [ Vx * Sx - ct,   Vx * Sy - st,   Vx * Vz],
    [ Vy * Sx + st,   Vy * Sy - ct,   Vy * Vz],
    [ Vz * Sx     ,   Vz * Sy     ,   1.0 - z]]
    '''

def quat2mat(q):
    x=q[0]; y=q[1]; z=q[2]; w=q[3];
    r2 = x*x + y*y + z*z + w*w;
    s  = 2 / r2;
    xs = x * s;  ys = y * s;  zs = z * s;
    xx = x * xs; xy = x * ys; xz = x * zs;
    xw = w * xs; yy = y * ys; yz = y * zs;
    yw = w * ys; zz = z * zs; zw = w * zs;
    return np.array(  [
        [1 - (yy + zz),     (xy - zw),     (xz + yw) ],
        [    (xy + zw), 1 - (xx + zz),     (yz - xw) ],
        [    (xz - yw),     (yz + xw), 1 - (xx + yy) ]
    ] )

def mat2quat(m):
    t = m[0,0] + m[1,1] + m[2,2];
    if (t >= 0):
        s = np.sqrt(t + 1);
        w = 0.5 * s;
        s = 0.5 / s;
        x =  (m[2,1] - m[1,2]) * s;
        y =  (m[0,2] - m[2,0]) * s;
        z =  (m[1,0] - m[0,1]) * s;
    elif ((m[0,0] > m[1,1]) and (m[0,0] > m[2,2])):
        s = np.sqrt(1 + m[0,0] - m[1,1] - m[2,2]);
        x = s * 0.5;
        s = 0.5 / s;
        y = (m[1,0] + m[0,1]) * s;
        z = (m[0,2] + m[2,0]) * s;
        w = (m[2,1] - m[1,2]) * s;
    elif (m[1,1] > m[2,2]):
        s = np.sqrt(1 + m[1,1] - m[0,0] - m[2,2]);
        y = s * 0.5;
        s = 0.5 / s;
        x = (m[1,0] + m[0,1]) * s;
        z = (m[2,1] + m[1,2]) * s;
        w = (m[0,2] - m[2,0]) * s;
    else:
        s = np.sqrt(1 + m[2,2] - m[0,0] - m[1,1]);
        z = s * 0.5;
        s = 0.5 / s;
        x = (m[0,2] + m[2,0]) * s;
        y = (m[2,1] + m[1,2]) * s;
        w = (m[1,0] - m[0,1]) * s;
    return np.array([x,y,z,w])

def initSurf( surfFile, cell, ns=[60,60,100] ):
    rmol.initRigidSubstrate ( surfFile, np.array(ns,dtype=np.int32), np.array([0.0,0.0,0.0]), np.array(cell) )
    if os.path.isfile("data/FFPauli.bin"):
        print "gridFF found on disk => loading "
        rmol.loadGridFF()
    else:
        print "gridFF not found on disk => recalc "
        rmol.recalcGridFF( np.array([1,1,1],dtype=np.int32) )
        rmol.saveGridFF()
    #rmol.debugSaveGridFF( "FFtot_z_Na.xsf", np.array([1.3,0.0447214,0.0]) )

def getSurfConfs( rots, molFile, pos=[ 5.78, 6.7, 12.24 ], nMaxIter=200, Fconv=0.01 ):

    #print "DEBUG 0"

    rmol.clear()

    #mol   = rmol.loadMolType( molFile )                              ; #print "DEBUG 0.1"

    atomTypes =  np.array(  [ atomTypeNames[e] for e in es[:nAtomMol] ], dtype=np.int32 )
    print atomTypes
    apos = xyzs[:nAtomMol,:]
    #REQs = np.array( [ [1.5, 0.01, 0.0], ]*len(apos) )                  # REQs are not set
    REQs = assignREQ( es[:nAtomMol], atomTypes_dct )
    #print "REQs ======= ", REQs, " apos.shape ", apos.shape, " qs.shape ", qs.shape
    REQs[:,2] = qs[:]
    #print "REQs ======= ", REQs, " apos.shape ", apos.shape, " qs.shape ", qs.shape
    #exit()
    #print REQs; exit()

    mol = rmol.registerRigidMolType( apos, REQs, atomTypes )

    rot0  = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])    ; #print "DEBUG 0.2"
    rmol.insertMolecule( mol, np.array(pos), rot0, True )            ; #print "DEBUG 0.3"

    #exit()

    #print "DEBUG 1"

    #rmol.save2xyz( "world_debug_00.xyz" )

    # ========= Relaxation

    rmol.bakeMMFF();   #print "DEBUG 1.1"
    rmol.prepareOpt(); #print "DEBUG 1.2"
    rmol.setOptFIRE( dt_max=0.2, dt_min=0.01, damp_max=0.1, minLastNeg=5, finc=1.1, fdec=0.5, falpha=0.98, kickStart=1.0 ); #print "DEBUG 1.3"

    #print "DEBUG 3"

    poses = rmol.getPoses();    #print "rmol.getPoses() ", poses_
    apos  = rmol.getAtomPos();  #print "rmol.getAtomPos() ", apos

    mol_name = molFile.split("/")[1].split(".")[0]

    '''
    rots_ = []
    for irot,rot in enumerate(rots):
        q = mat2quat(rot)
        poses[0,4:8] = q
        F2 = rmol.relaxNsteps( nMaxIter, Fconv**2 ); 
        print "irot, F2 ", irot, F2
        rot_ = quat2mat(poses[0,4:8])
        rots_.append(rot_)
        #print "rot ", irot, rot,"\n -> ", rot_
    '''

    '''
    rots_ = []
    fout = open( "movie_%s_rots.xyz" %(mol_name) ,'w')
    for irot,rot in enumerate(rots):
        q = mat2quat(rot)
        poses[0,4:8] = q
        F2 = rmol.relaxNsteps( 100, Fconv**2 ); 
        rot_ = quat2mat(poses[0,4:8])
        isim = findSimilarRot(rot_,rots_, dcut=0.1, ipiv=2)
        print "irot, F2 ", irot, F2, isim
        if(isim<0):
            xyzs[:nAtomMol,:] = apos[:,:]
            au.writeToXYZ( fout, es, xyzs )
            rots_.append(rot_)
    fout.close()
    print "len(rots_)", len(rots_)
    '''

    
    rots_ = []
    for irot,rot in enumerate(rots):
        #fout = rmol.openf( "movie.xyz", -1, "w" )
        mol_name = molFile.split("/")[1].split(".")[0]
        print mol_name
        fout = open( "movie_%s_%03i.xyz" %(mol_name,irot) ,'w')
        q = mat2quat(rot)
        print "q ", q
        poses[0,4:8] = q
        for i in range(nMaxIter):
            #F2 = rmol.relaxNsteps( nMaxIter, Fconv**2 ); 
            F2 = rmol.relaxNsteps( 1, 0.0 ); 
            #rot_ = quat2mat(poses[0,4:8])
            #rots_.append(rot_)
            print ">>> i ", i, poses
            #print "|F| ", np.sqrt(F2)
            xyzs[:nAtomMol,:] = apos[:,:]
            au.writeToXYZ( fout, es, xyzs )
            #rmol.write2xyz( fout )
            #rmol.save2xyz( "world_debug_%03i.xyz" %i )
        #print "rot  ", rot
        #print "rot_ ", rot_
        #print "irot  -", irot
        fout.close()
    

    #del  poses
    #del  apos
    return rots_

#  >> itr 0 F2 0.557349 dt 0.05 qrot (-0.353364,-0.352836,-0.612781,0.612486) int 139984312000528 




if __name__ == "__main__":
    print " ================ START "
    print " ================ START "
    print " ================ START "
    print " ================ START "
    print " ================ START "

    #os.chdir( "/u/25/prokoph1/unix/git/SimpleSimulationEngine/cpp/Build/apps/MolecularEditor2" )
    os.chdir( "/home/prokop/git/SimpleSimulationEngine/cpp/Build/apps/MolecularEditor2")

    #FFparams=np.genfromtxt("inputs/atomtypes.ini",dtype=[('rmin',np.float64),('epsilon',np.float64),('alpha',np.float64),('atom',np.int),('symbol', '|S10')],usecols=(0,1,2,3,4))
    #print "FFparams[:,3]", FFparams[:,3]
    #print "FFparams ", FFparams;  

    #"H   1 1 1 0     1.487   0.0006808   0xFFFFFF"
    atomTypes_dat =np.genfromtxt("common_resources/AtomTypes.dat",dtype=[('symbol', '|S4'),('i1',np.int32),('i2',np.int32),('i3',np.int32),('i4',np.int32),('R',np.float64),('eps',np.float64),('alpha',np.int32)],usecols=(0,1,2,3,4,5,6))
    print "atomTypes_dat ", atomTypes_dat

    #for a in atomTypes_dat:
    #    print a
    #    print a[0]
    #    print a[1]
    atomTypes_dct = { a[0]:(a[5],a[6]) for a in atomTypes_dat }
    print "atomTypes_dct ", atomTypes_dct

    water   = au.loadAtoms( "inputs/water_T5_ax.xyz" );  #print Campher
    campher = au.loadAtoms( "inputs/Campher.xyz" );      #print Campher
    surf    = au.loadAtoms( "inputs/Cu111_6x6_2L.xyz" ); #print Surf

    #print water
    #exit()

    cell = [[15.31593,0.0,0.0],[0.0,13.26399,0.0],[0.0,0.0,20.0]]

    '''
    #rots  = sphereTangentSpace(n=100)
    rots  = randomRotations(n=10000)
    rots = pruneRots(rots, dcut=0.3, ipiv=2)
    #rots = getSimilarRotations( rots, rots[0], dcut=0.1, ipiv=2)
    print "len(rots)", len(rots)
    #print rots
    drawRots( rots, lenght=0.1)
    plt.show()
    #exit()
    '''
    rots  = randomRotations(n=5)

    print " rmol.initParams( ) "
    rmol.initParams( "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat" )
    atomTypeNames = rmol.loadAtomTypeNames( "common_resources/AtomTypes.dat" )
    print "atomTypeNames", atomTypeNames
    initSurf( "inputs/Cu111_6x6_2L.xyz", cell, ns=[60,60,100] )
    
    rmol.setCoulombMirror( np.array([0.0,0.0,1.0]), np.array([0.0,0.0,7.0]) )
    
    '''
    print "========== water_T5_ax.xyz ==========="
    print "========== water_T5_ax.xyz ==========="
    print "========== water_T5_ax.xyz ==========="
    nAtomMol = len(water[0])
    es, xyzs = combineGeoms(water,surf)
    #es, xyzs = convGeom(water)
    rots_ = getSurfConfs( rots, "inputs/water_T5_ax.xyz", pos=[ 5.78, 6.7, 10.00 ],  nMaxIter=100, Fconv=0 )
    rots_ = pruneRots(rots_, dcut=0.01, ipiv=2)
    print "len(rots_)", len(rots_)
    '''

    #print "========== Campher.xyz ==========="
    #print "========== Campher.xyz ==========="
    #print "========== Campher.xyz ==========="
    nAtomMol = len(campher[0])
    qs       = np.array(campher[4])
    print "qs ", qs
    es, xyzs = combineGeoms(campher,surf)
    #es, xyzs = convGeom(campher)
    rots_ = getSurfConfs( rots, "inputs/Campher.xyz", pos=[ 5.78, 6.7, 12.24 ], nMaxIter=100, Fconv=0 )
    #rots_ = pruneRots(rots_, dcut=0.01, ipiv=2)
    #drawRots( rots, lenght=0.1)
    #plt.show()
    print "len(rots_)", len(rots_)

    print rots_

    #print "rots_", rots_
    print ">>>> ALL DONE <<<<"




