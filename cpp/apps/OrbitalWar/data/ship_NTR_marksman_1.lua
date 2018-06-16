
require( "data/lua/utils" )
print("BEGIN spaceship");

origin = {0.0,0.0,0.0}
xvec   = {1.0,0.0,0.0}
yvec   = {0.0,1.0,0.0}
zvec   = {0.0,0.0,1.0}

n0 = Node( origin );
n1 = Node( {-50.0,   0.0,    0.0} )
n2 = Node( { 50.0,   0.0,    0.0} )
n3 = Node( {   0.0,-100.0,    0.0} )
n4 = Node( {   0.0, 100.0,    0.0} )
n5 = Node( {   0.0,   0.0, -800.0} )
n6 = Node( {   0.0,   0.0,  200.0} )

defGirderWidth = 1.5

--         from to  Up nseg,mseg   width, hegith
g1 = Girder( n0, n1, zvec, 10, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g2 = Girder( n0, n2, zvec, 10, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g3 = Girder( n0, n3, xvec, 20, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g4 = Girder( n0, n4, xvec, 20, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g5 = Girder( n0, n5, xvec, 30, 2, {8.0,8.0}, "Steel" )
g6 = Girder( n0, n6, xvec, 20, 2, {8.0,8.0}, "Steel" )

g7 = Girder( n3, n5, xvec, 70, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g8 = Girder( n4, n5, xvec, 70, 2, {defGirderWidth,defGirderWidth}, "Steel" )

--          type     thick[mm]
Rope(n5,n1, 25, "Kevlar");
Rope(n6,n1, 25, "Kevlar" )
Rope(n5,n2, 25, "Kevlar");
Rope(n6,n2, 25, "Kevlar" )
Rope(n6,n3, 25, "Kevlar");
Rope(n6,n4, 25, "Kevlar");

Rope(n1,n3, 25, "Kevlar");
Rope(n1,n4, 25, "Kevlar");
Rope(n2,n3, 25, "Kevlar");
Rope(n2,n4, 25, "Kevlar");

tanks( 6,0.0, 10.0, 5.0, -50.0, -20.0 )
girderFan( 4, 0.0, 20.0, 40.0,-50.0, -250, 10, 1.0, true )

girderFan( 4, 0.0, 10.0, 20.0,-650.0, -750, 10, 1.0, true )
girderFan( 4, 0.0, 15.0, 30.0,-350.0, -500, 10, 1.0, true )

Tank( {0.0,0,-800}, zvec, {8.0,1.5,40.0}, "H2")

--Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
Radiator( g5,0.1,0.8, g7,0.02,0.8, 1280.0 )
Radiator( g5,0.1,0.8, g8,0.02,0.8, 1280.0 )

--      node1,2, up,  nseg    R     {width,height}
Ring( {0.0,0.0,16.0}, zvec, xvec, 16, 20.0, {1.0,1.0}, "Steel" )
Ring( {0.0,0.0,-300.0}, xvec, yvec, 16, 60.0, {1.0,1.0}, "Steel" )
-- Ring( {0.0,0.0,0.0}, yvec, xvec, 16, 160.0, {8.0,5.0}, "Steel" )

--  There should be mechanism how to generate nodes on-top of ship components (anchor points)

Gun( g6, 0.1, 0.8, "XFEL" )

print("END spaceship");
