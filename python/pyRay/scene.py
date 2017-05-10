#!/usr/bin/python

import os
import re
import numpy    as np 

from common import *

def makeProgram( src_scene, src_user_func="" ):
    src_primitives = loadToString(CL_PATH+'/primitives.cl')
    src_operations = loadToString(CL_PATH+'/operations.cl')
    src_raytrace   = loadToString(CL_PATH+'/rayScene_basic.cl', dct=defines )
    src_raytrace   = re.sub( "//===SCENE"         , src_scene,     src_raytrace, count=1, flags=0)
    src_raytrace   = re.sub( "//===USER_FUNCTIONS", src_user_func, src_raytrace, count=1, flags=0)
    return src_primitives + src_operations + src_raytrace

def setDefines( src, dct ):
    for il,line in enumerate(src):
        if line.startswith( '#define' ):
            for key in dct:
                if key in line:
                    words    = line.split()
                    words[2] = str( dct[key] )+"\n"
                    src[il]  = " ".join( words )
    return src            

def loadToString( fname, dct = None ):
    with open(fname, 'r') as f:
        lines = f.readlines()
        if dct is not None:
            lines=setDefines(lines, dct)
    return "".join(lines)
    
def parseArg(arg):
    typ = type(arg)
    #print typ
    if typ   is float:
        return "%ff" %arg
    elif typ is int:
        return "%i" %arg
    elif (typ is tuple) or (typ is list) or (typ is np.ndarray):
        n = len(arg)
        #print n
        if n == 2:
            return "(float2)(%ff,%ff)" %arg
        elif n == 3:
            return "(float3)(%ff,%ff,%ff)" %arg
        elif n == 4:
            return "(float4)(%ff,%ff,%ff,%ff)" %arg
        elif n == 8:
            return "(float8)(%ff,%ff,%ff,%ff,%ff,%ff,%ff,%ff)" %arg
        else:
            return "ERROR1"
    else:
        return "ERROR0"

def parseObject( obj ):
    sargs = ",".join( [ parseArg(arg) for arg in obj[3] ] )
    spos  = obj[2] %"pos"
    #return "%s(%s(%s,%s));" %(obj[0],obj[1],spos, sargs)
    return "dist=op%s(dist,%s(%s,%s));" %(obj[0],obj[1],spos,sargs)

def parseSceneList( lst ):
    return "\n".join([ parseObject( obj ) for obj in lst] )


