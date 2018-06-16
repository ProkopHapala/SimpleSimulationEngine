require( "data/lua/utils" )

n0 = Node( origin );
n1 = Node( { 0.0,0.0,600.0} )
n2 = Node( { 0.0,0.0,-200.0} );
n3 = Node( { -200.0,0.0,0.0} );
n4 = Node( {  200.0,0.0,0.0} );

g1 = Girder( n2, n1, xvec, 30, 2, {8.0,8.0}, "Steel" )
g2 = Girder( n3, n1, xvec, 30, 2, {8.0,8.0}, "Steel" )
g3 = Girder( n4, n1, xvec, 30, 2, {8.0,8.0}, "Steel" )
Gun( g1, 0.0, 1.0, "XFEL" )
tanks( 6,0.0, 10.0, 5.0, 100.0, 500.0 )

-- girderFan( 6, 0.0, 20.0, 20.0,200.0, 20, 10, 1.0, true )
-- Ring( {0.0,0.0,4.0}, yvec, xvec, 16, 100.0, {4.0,4.0}, "Steel" )
Thruster( {0.0,0.0,0.0}, yvec, {5.0,300.0,20.0}, "ICF_Ebeam_magNozzle" )
