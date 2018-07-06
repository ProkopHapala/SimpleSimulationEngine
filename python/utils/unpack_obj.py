#!/usr/bin/python

fin =  open( "lightTransportTest.obj", "r" )

def line2vert( line ):
    wds = line.split()
    return [float(wds[1]),float(wds[2]),float(wds[3])]

def line2face( line ):
    wds = line.split()
    fc = []
    for wd in wds[1:]:
        #fc.append( wd.split("/")[0] )
        fc.append( int(wd.split("/")[0])-1 )
    return  fc

def faces2tris(faces):
    tris = []
    for f in faces:
        for i in range(len(f)-2):
            #tris.append(f[i:i+3])
            tris.append( [f[0],f[i+1],f[i+2]] )
    return tris

def formatCarray( lst ):
    s = "{"
    i = 0
    n1=len(lst)-1
    for item in lst:
        if type(item) is list:
            s +=  formatCarray( item )
        else:
            s += str(item)
        if i<n1: s+= ","
        i+=1
    return s+"}"


verts = []
faces = []

for line in fin:
    if line [:2]== "v ":
        verts.append( line2vert( line ) )
    if line [:2]== "f ":
        faces.append( line2face( line ) )


'''
print verts

print faces
print "=============="
print faces2tris(faces)


'''



print "verts="
print formatCarray( verts )

tris=faces2tris(faces)
print "tris="
print formatCarray( tris )

print "nVerts",len(verts) 
print "nTris",len(tris)
print "nFaces",len(faces)

