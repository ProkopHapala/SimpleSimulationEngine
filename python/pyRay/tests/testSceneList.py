#!/usr/bin/python

import sys
sys.path.append("../../")
#import pyRay       as ra
import pyRay.scene as scn 

# TODO : how to pass arguments from function header?
object1 = ("obj1",(),        [( "U","sdBox"   ,"%s",((1.0,1.0,1.0),) ),( "S","sdSphere","%s",(1.2,) )])
object2 = ("obj1",("f","f3"),[( "U","sdBox"   ,"%s",("2",)           ),( "S","sdSphere","%s",("1",) )])
object3 = ("obj2",("f","f2"),[( "U","sdBox"   ,"%s",(("2",1.0),)     ),( "S","sdSphere","%s",("1",) )])

scene   = [
( "U","sdBox"   ,"%s",((1.0,1.0,1.0),) ),
( "S","sdSphere","%s",(1.2,) ),
]

scene_src = scn.parseSceneList(scene)

print scene_src




