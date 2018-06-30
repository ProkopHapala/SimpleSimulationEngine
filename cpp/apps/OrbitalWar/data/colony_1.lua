--- https://www.lua.org/pil/5.3.html

require( "data/lua/utils" )

print("BEGIN spaceship");

-- http://www.mse.mtu.edu/~drjohn/my4150/props.html

--- newSpaceShip( "Ship1" )

n1 = Node( {-300.0, -300.0,     0.0} )
n2 = Node( { 300.0, -300.0,     0.0} )
n3 = Node( {   0.0,  300.0,  -300.0} )
n4 = Node( {   0.0,  300.0,   300.0} )

n5 = Node( {   800.0,  800.0,   800.0} )

n6 = Node( {   800.0,  800.0,   -800.0} )


--[[
n6 = Node( {   1000.0,  1000.0,   1000.0} )
n7 = Node( {   1000.0,  1000.0,    900.0} )
n8 = Node( {   1000.0,  1000.0,   2000.0} )
n9 = Node( {   1000.0,  1000.0,   0.0} )
]]

--         from to  Up nseg,mseg   width, hegith
g1 = Girder( n1, n2, zvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g2 = Girder( n1, n3, zvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g3 = Girder( n1, n4, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g4 = Girder( n2, n3, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g5 = Girder( n3, n4, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g6 = Girder( n4, n2, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )

g4 = Girder( n5, n3, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g5 = Girder( n5, n4, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g6 = Girder( n5, n2, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )

g6 = Girder( n5, n6, xvec, 100, 2, {defGirderWidth,defGirderWidth}, "Steel" )

--[[
g6 = Girder( n5, n2, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
g6 = Girder( n5, n2, xvec, 50, 2, {defGirderWidth,defGirderWidth}, "Steel" )
]]



Rock( {-300.0, -300.0,     0.0}, xvec, yvec, {130.0,120.0,140.0} );
Rock( { 300.0, -300.0,     0.0}, xvec, yvec, {140.0,130.0,145.0} );
Rock( {   0.0,  300.0,  -300.0}, xvec, yvec, {150.0,140.0,145.0} );
Rock( {   0.0,  300.0,   300.0}, xvec, yvec, {160.0,150.0,130.0} );


Rock( {   800.0,  800.0,   800.0} , xvec, yvec, {260.0,250.0,230.0} );


Balloon( {0.0,0.0,0.0}, xvec, yvec, {700.0,700.0,700.0} );
