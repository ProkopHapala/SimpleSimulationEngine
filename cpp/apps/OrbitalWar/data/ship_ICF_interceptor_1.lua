require( "data/lua/utils" )

--- long, perp, zigIn, zigOut
st1  = StickMaterial( "GS1_long", "Steel", 0.1,  0.005 )
st2  = StickMaterial( "GS1_perp", "Steel", 0.05, 0.003 )
st3  = StickMaterial( "GS1_in",   "Steel", 0.04, 0.002 )
st4  = StickMaterial( "GS1_out",  "Steel", 0.04, 0.002 )
-- stK1 = StickMaterial( "RK_1","Kevalr", 0.01, 0.001 )

g1_st = { st1, st2, st3, st4 };

n0  = Node( origin );
n1 = Node( { 0.0,0.0,300.0} )


g1 = Girder( n0, n1, xvec, 30, 2, {8.0,8.0}, "Steel", g1_st )
Gun( g1, 0.0, 1.0, "XFEL" )
tanks( 6,0.0, 10.0, 5.0, 100.0, 220.0 )

-- girderFan( 4, 0.0, 20.0, 40.0,200.0, 20, 10, 1.0, true )
girderFan( 6, 0.0, 20.0, 20.0,200.0, 20, 10, 1.0, true )

girderFan( 8, 0.0, 260.0, 260.0, -260.0, 200, 10, 4.0, false )

Thruster( {0.0,0.0,0.0}, zvec, {5.0,400.0,600.0}, "ICF_Ebeam_magNozzle" )
