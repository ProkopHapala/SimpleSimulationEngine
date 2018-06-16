require( "data/lua/utils" )

origin = {0.0,0.0,0.0}
xvec   = {1.0,0.0,0.0}
yvec   = {0.0,1.0,0.0}
zvec   = {0.0,0.0,1.0}

n0  = Node( origin );
n1 = Node( { 0.0,0.0,300.0} )

g1 = Girder( n0, n1, xvec, 30, 2, {8.0,8.0}, "Steel" )
Gun( g1, 0.0, 1.0, "XFEL" )
tanks( 6,0.0, 10.0, 5.0, 100.0, 220.0 )

-- girderFan( 4, 0.0, 20.0, 40.0,200.0, 20, 10, 1.0, true )
girderFan( 6, 0.0, 20.0, 20.0,200.0, 20, 10, 1.0, true )

girderFan( 8, 0.0, 260.0, 260.0, -260.0, 200, 10, 4.0, false )

Thruster( {0.0,0.0,0.0}, zvec, {5.0,400.0,600.0}, "ICF_Ebeam_magNozzle" )
