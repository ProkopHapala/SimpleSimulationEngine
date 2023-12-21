--- https://www.lua.org/pil/5.3.html

require( "data/lua/utils" )

print("BEGIN spaceship");

-- http://www.mse.mtu.edu/~drjohn/my4150/props.html

--- newSpaceShip( "Ship1" )

--- long, perp, zigIn, zigOut
st1  = StickMaterial( "GS1_long", "Steel", 0.1,  0.005 )
st2  = StickMaterial( "GS1_perp", "Steel", 0.05, 0.003 )
st3  = StickMaterial( "GS1_in",   "Steel", 0.04, 0.002 )
st4  = StickMaterial( "GS1_out",  "Steel", 0.04, 0.002 )
-- stK1 = StickMaterial( "RK_1","Kevalr", 0.01, 0.001 )

n0 = Node( origin );
n1 = Node( {-100.0,   0.0,    0.0} )
n2 = Node( { 100.0,   0.0,    0.0} )
n3 = Node( {   0.0,-200.0,    0.0} )
n4 = Node( {   0.0, 200.0,    0.0} )
n5 = Node( {   0.0,   0.0, -300.0} )
n6 = Node( {   0.0,   0.0,  800.0} )

print( "Lua:Nodes:" , n0,n1,n2,n3,n4,n5,n6 )

defGirderWidth = 4.0

g1_st = { st1, st2, st3, st4 };

print( "Lua:g1_st:" , st1, st2, st3, st4 )

--         from to  Up nseg,mseg   width, hegith
g1 = Girder( n0, n1, zvec, 10, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
g2 = Girder( n0, n2, zvec, 10, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
g3 = Girder( n0, n3, xvec, 20, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
g4 = Girder( n0, n4, xvec, 20, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
g5 = Girder( n0, n5, xvec, 30, 2, {8.0,8.0}, "Steel", g1_st )
g6 = Girder( n0, n6, xvec, 80, 2, {8.0,8.0}, "Steel", g1_st )

g7 = Girder( n3, n6, xvec, 70, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
g8 = Girder( n4, n6, xvec, 70, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )

--          type     thick[mm]
Rope(n5,n1, 25, "Kevlar" );
Rope(n6,n1, 25, "Kevlar" )
Rope(n5,n2, 25, "Kevlar" );
Rope(n6,n2, 25, "Kevlar" )
Rope(n5,n3, 25, "Kevlar" );
--Rope(n6,n3, 25, "Kevlar" )
Rope(n5,n4, 25, "Kevlar");
--Rope(n6,n4, 25, "Kevlar" )

Rope(n1,n3, 25, "Kevlar");
Rope(n1,n4, 25, "Kevlar");
Rope(n2,n3, 25, "Kevlar");
Rope(n2,n4, 25, "Kevlar");

tanks( 6,0.0, 10.0, 5.0, -50.0, -20.0 )

print( "#------ Lua:girderFans:" )

girderFan( 4, 0.0, 20.0, 40.0,-50.0, -250, 10, 1.0, true, g1_st  )
girderFan( 4, 0.0, 20.0, 40.0, 20.0,  300, 10, 1.0, true, g1_st  )
girderFan( 8, 0.0, 80.0, 80.0, -380, -200, 10, 1.0, false, g1_st )

-- =RadiatorType = {"LithiumHeatPipe", 1280.0 }

--Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
Radiator( g6,0.15,0.8, g7,0.02,0.8, 1280.0 )
Radiator( g6,0.15,0.8, g8,0.02,0.8, 1280.0 )

--      node1,2, up,  nseg    R     {width,height}
Ring( {0.0,0.0,4.0}, zvec, xvec, 16, 100.0, {4.0,4.0}, "Steel", g1_st )
Ring( {20.0,0.0,0.0}, xvec, yvec, 16, 108.0, {4.0,4.0}, "Steel", g1_st )
-- Ring( {0.0,0.0,0.0}, yvec, xvec, 16, 160.0, {8.0,5.0}, "Steel" )

--  There should be mechanism how to generate nodes on-top of ship components (anchor points)

Thruster( {0.0,0.0,-300.0}, zvec, {5.0,100.0,120.0}, "ICF_Ebeam_magNozzle" )
-- Thruster( {-16,-16,16}, {1.0,2.0,3.0}, {5.0,100.0,50.0}, "ICF_Ebeam_magNozzle" )


Gun( g6, 0.1, 0.8, "XFEL" )


print("END spaceship");

--[=====[
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

Ring  ( n0, n6, xvec, 80, 2, {10.0,8.0} )

--                   path    side offset
Gun ( "railgun", {n5,n0,n6}, { 1.0,0.0} )
Gun ( "XFEL",    {n5,n0,n6}, {-1.0,0.0} )

-- Tank  ( "", R=16 )
-- Tank  ( "water", R=16 )

Thruster{ n5, maxPower=7.0e+9, maxTrust=10.0e+3, maxVexh=50e+3, nodes={} }

--]=====]

