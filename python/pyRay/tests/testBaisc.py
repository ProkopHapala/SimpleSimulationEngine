#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt 
sys.path.append("../../")
import pyRay       as ra
import pyRay.scene as scn 
import pyRay.ocl   as ocl 
import pyRay.image as img 

src_scene = '''
    dist = sdSphere( pos, 0.5f );   
    dist = opU( dist, sdSphere( pos-(float3)(0.7f,0.0f,0.2f), 0.7f ) ); 
'''

src_scene = ''' 
    dist = opU( dist, sdSphere  ( pos- (float3)(0.0f,0.5f,0.0), 0.5f )  );
    dist = opU( dist, sdCone    ( pos, (float3)(0.8f,0.6f,0.3f ) ) );
    dist = opU( dist, sdCylinder( pos, (float2)(0.1f,0.2f)       ) );
    dist = opS( dist, sdTorus   ( pos- (float3)(0.0f,0.1f,0.0), (float2)(0.25f,0.1f)     ) );
    dist = opS( dist, sdPlane( pos, (float3)(1.0,0.0,0.0), 0.0 ) ); 
    dist = opS( dist, sdCylinder( pos, (float2)(0.05f,0.9f)       ) );
    dist = opS( dist, sdSphere  ( pos- (float3)(0.0f,0.5f,0.0), 0.3f )  );
'''

#ra.screen_resolution[0] = 16
#ra.screen_resolution[1] = 16
#ra.screen_resolution[0] = 1024
#ra.screen_resolution[1] = 1024
#ra.defines["HIT_PREC"]  = 0.01
#ra.defines["MAX_STEPS"] = 16

prog_scr = scn.makeProgram( src_scene )
ocl.make_program(prog_scr)

fw = np.array([0.5,0.5,1.0])
light_dir = np.array([0.5,-1.0,0.0]); light_dir /=np.sqrt(np.dot(light_dir,light_dir));
spec_vec  = light_dir+fw

cam       = ra.getCamMat          ( fw )
rd, ro    = ocl.getRays           ( cam, fw0=-4.0, t0=8.0, tg=(0.25,0.25) )
print "rd.nbytes", rd.nbytes
kargs     = ocl.prep_rayTraceBasic( )
hitp,hitn = ocl.run_rayTraceBasic ( rd, ro, kargs )

#print hitp[:,:,3]
#print hitn

mask = hitp[:,:,3] > 20.0 
hitp[mask,3] = np.NaN

plt.figure(figsize=(10,5))
plt.subplot(1,2,1); plt.imshow( hitn[:,:,3] , interpolation='nearest' ); plt.colorbar()
#plt.subplot(1,2,1); plt.imshow( np.log( hitp[:,:,3] ), interpolation='nearest' )
#plt.subplot(1,2,2); plt.imshow( hitn[:,:,0:3]*0.5+0.5, interpolation='nearest' )

AO = 1/(1+0.01*hitn[:,:,3]**2) 

diffuse  = img.getDiffuse ( hitn, light_dir )
specular = img.getSpecular( hitn, light_dir, rd, gloss=256.0, power=2 )
#plt.subplot(1,2,2); plt.imshow( diffuse*0.5+0.5, interpolation='nearest', cmap='gray' )
#plt.subplot(1,2,2); plt.imshow( specular, interpolation='nearest', cmap='gray' )
#plt.subplot(1,2,2); plt.imshow( (diffuse+specular)*0.5+0.5, interpolation='nearest', cmap='gray' )
#plt.subplot(1,2,2); plt.imshow( (diffuse+specular)*0.5+0.5+0.5*AO, interpolation='nearest', cmap='gray' )
plt.subplot(1,2,2); plt.imshow( (diffuse+AO)*0.5+0.5, interpolation='nearest', cmap='gray' )
#plt.subplot(1,2,2); plt.imshow( +0.5*AO, interpolation='nearest', cmap='gray' )


#AO = img.getAO ( hitp )
#plt.subplot(1,2,2); plt.imshow( AO, interpolation='nearest', cmap='gray' )

plt.show()
