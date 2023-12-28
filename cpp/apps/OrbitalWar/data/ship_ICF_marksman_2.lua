--- https://www.lua.org/pil/5.3.html

require( "data/lua/utils" )

print("BEGIN spaceship");

-- http://www.mse.mtu.edu/~drjohn/my4150/props.html

--- newSpaceShip( "Ship1" )

print( "#------ Lua:StickMaterials" )

--- long, perp, zigIn, zigOut
st1  = StickMaterial( "GS1_long", "Steel", 100, 5 )
st2  = StickMaterial( "GS1_perp", "Steel", 50,  3 )
st3  = StickMaterial( "GS1_in",   "Steel", 40,  2 )
st4  = StickMaterial( "GS1_out",  "Steel", 40,  2 )
-- stK1 = StickMaterial( "RK_1","Kevalr", 0.01, 0.001 )

print( "#------ Lua:Nodes" )

n0  = Node( origin );
nLf = Node( {-100.0,   0.0,    0.0} )
nRt = Node( { 100.0,   0.0,    0.0} )
nDw = Node( {   0.0,-200.0,    0.0} )
nUp = Node( {   0.0, 200.0,    0.0} )
nBk = Node( {   0.0,   0.0, -300.0} )
nFw = Node( {   0.0,   0.0,  800.0} )

print( "#------ Lua:Girders" )


print( "Lua:Nodes:" , n0,nLf,nRt,nDw,nUp,nBk,nFw )

defGirderWidth = 4.0

g1_st = { st1, st2, st3, st4 };

print( "Lua:g1_st:" , st1, st2, st3, st4 )

--         from to  Up nseg,mseg   width, hegith
gLf = Girder( n0, nLf, zvec, 10, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
gRt = Girder( n0, nRt, zvec, 10, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
gDw = Girder( n0, nDw, xvec, 20, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
gUp = Girder( n0, nUp, xvec, 20, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
gBk = Girder( n0, nBk, xvec, 30, 2, {8.0,8.0}, "Steel", g1_st )
gFw = Girder( n0, nFw, xvec, 80, 2, {8.0,8.0}, "Steel", g1_st )

gWdw = Girder( nDw, nFw, xvec, 70, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )
gWup = Girder( nUp, nFw, xvec, 70, 2, {defGirderWidth,defGirderWidth}, "Steel", g1_st )

print( "#------ Lua:Ropes" )

--          type     thick[mm]
Rope(nBk,nLf, 25, "Kevlar" );
Rope(nFw,nLf, 25, "Kevlar" )
Rope(nBk,nRt, 25, "Kevlar" );
Rope(nFw,nRt, 25, "Kevlar" )
Rope(nBk,nDw, 25, "Kevlar" );
--Rope(nFw,nDw, 25, "Kevlar" )
Rope(nBk,nUp, 25, "Kevlar");
--Rope(nFw,nUp, 25, "Kevlar" )

Rope(nLf,nDw, 25, "Kevlar");
Rope(nLf,nUp, 25, "Kevlar");
Rope(nRt,nDw, 25, "Kevlar");
Rope(nRt,nUp, 25, "Kevlar");

print( "#------ Lua:BoundNodes and Girders" )

--- Finished basic layout of the ship ( main girders and ropes )

nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {0.0,100.0,100.0} );
--nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {0.0,-100.0,0.0} );
--nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {100.0,0.0,0.0} );
--nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {-100.0,0.0,0.0} );
nB2 = BoundNode( gUp, ComponetKind.Girder, 0.5, {0.0,100.0,100.0} );
gB  = Girder( nB1, nB2, xvec, 10, 2, {4.0,4.0}, "Steel", g1_st )

print( "#------ Lua:Radiataors" )

-- =RadiatorType = {"LithiumHeatPipe", 1280.0 }
--Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
--Radiator( gFw,0.15,0.8, gWdw,0.02,0.8, 1280.0 )
--Radiator( gFw,0.15,0.8, gWup,0.02,0.8, 1280.0 )

print( "#------ Lua:Rings & Sliders" )

--      node1,2, up,  nseg    R     {width,height}
r_roll  = Ring( {0.0,0.0,8.0}, zvec, xvec, 64, 55.0, {4.0,4.0}, "Steel", g1_st )
r_pitch = Ring( {8.0,0.0,0.0}, xvec, yvec, 64, 80.0, {16.0,4.0}, "Steel", g1_st )
-- Ring( {0.0,0.0,0.0}, yvec, xvec, 16, 160.0, {8.0,5.0}, "Steel" )

--- print( "ComponetKind[Girder,Ring,Rope]", ComponetKind.Girder,ComponetKind.Ring,ComponetKind.Rope)
-- Slider( gFw,r_roll, ComponetKind.Girder,ComponetKind.Ring,   3,0, 0,0  )
-- Slider( r_pitch,gFw, ComponetKind.Ring,ComponetKind.Girder,   3,3, 0,0.1  )

--  There should be mechanism how to generate nodes on-top of ship components (anchor points)

print( "#------ Lua:Thrustrs,Tansk,Guns etc." )

Thruster( {0.0,0.0,-300.0}, zvec, {5.0,100.0,120.0}, "ICF_Ebeam_magNozzle" )
-- Thruster( {-16,-16,16}, {1.0,2.0,3.0}, {5.0,100.0,50.0}, "ICF_Ebeam_magNozzle" )
tanks( 6,0.0, 10.0, 5.0, -50.0, -20.0 )
Gun( gFw, 0.1, 0.8, "XFEL" )


print("END spaceship");
