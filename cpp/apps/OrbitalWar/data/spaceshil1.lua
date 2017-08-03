--- https://www.lua.org/pil/5.3.html

newSpaceShip( "Ship1" )

origin = {0.0,0.0,0.0}
xvec   = {1.0,0.0,0.0} 
yvec   = {1.0,0.0,0.0}
zvec   = {1.0,0.0,0.0}

n0 = Node( origin )
n1 = Node( {-100.0,   0.0,    0.0} )
n2 = Node( {+100.0,   0.0,    0.0} )
n3 = Node( {   0.0,-200.0,    0.0} )
n4 = Node( {   0.0,+200.0,    0.0} )
n5 = Node( {   0.0,   0.0, -300.0} )
n6 = Node( {   0.0,   0.0, +800.0} )

--         from to  Up nseg,mseg   width, hegith 
g1 = Girder( n0, n1, zvec, 10, 2, {10.0,8.0} )
g2 = Girder( n0, n2, zvec, 10, 2, {10.0,8.0} )
g3 = Girder( n0, n3, xvec, 20, 2, {10.0,8.0} )
g4 = Girder( n0, n4, xvec, 20, 2, {10.0,8.0} )
g5 = Girder( n0, n5, xvec, 30, 2, {10.0,8.0} )
g6 = Girder( n0, n6, xvec, 80, 2, {10.0,8.0} )

--          type     thick[mm]
M_rope1 = {"kevalr", 10   }

Rope(n5,n1,M_rope1); Rope(n6,n1,M_rope1)
Rope(n5,n2,M_rope1); Rope(n6,n2,M_rope1)
Rope(n5,n3,M_rope1); Rope(n6,n3,M_rope1)
Rope(n5,n4,M_rope1); Rope(n6,n4,M_rope1)

--                type          T[K]
M_radiator = {"LithiumHeatPipe", 1280.0 }

Radiator( g5, g1 )
Radiator( g5, g2 )
Radiator( g5, g3 )
Radiator( g5, g4 )
Radiator( g6, g1 )
Radiator( g6, g2 )
Radiator( g6, g3 )
Radiator( g6, g4 )

Ring  ( n0, n6, xvec, 80, 2, {10.0,8.0},  nodes={} )

--                   path    side offset
Gun ( "railgun", {n5,n0,n6}, { 1.0,0.0} )
Gun ( "XFEL",    {n5,n0,n6}, {-1.0,0.0} )

Tank  ( "", R=16 )
Tank  ( "water", R=16 )

Engine( n5, maxPower=7.0e+9, maxTrust=10.0e+3, maxVexh=50e+3, nodes={} )
