#!/usr/bin/python

import os
import re
import numpy    as np 

from common import *

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

def makeProgram( src_scene ):
    src_primitives = loadToString(CL_PATH+'/primitives.cl')
    src_operations = loadToString(CL_PATH+'/operations.cl')
    src_raytrace   = loadToString(CL_PATH+'/rayScene_basic.cl', dct=defines )
    src_raytrace   = re.sub( "//===SCENE", src_scene, src_raytrace, count=1, flags=0)
    return src_primitives + src_operations + src_raytrace
