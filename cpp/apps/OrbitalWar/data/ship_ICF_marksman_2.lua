--- https://www.lua.org/pil/5.3.html

require( "data/lua/utils" )

print("BEGIN spaceship");

-- http://www.mse.mtu.edu/~drjohn/my4150/props.html

--- newSpaceShip( "Ship1" )

print( "#------ Lua:StickMaterials" )

--- long, perp, zigIn, zigOut
st1  = StickMaterial( "GS1_long", "Steel", 100, 5, 0.0 )
st2  = StickMaterial( "GS1_perp", "Steel", 50,  3, 0.0 )
st3  = StickMaterial( "GS1_in",   "Steel", 40,  2, 0.0 )
st4  = StickMaterial( "GS1_out",  "Steel", 40,  2, 0.0 )
-- stK1 = StickMaterial( "RK_1","Kevalr", 0.01, 0.001, 0.0 )

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

g1_st = { st1, st2, st3, st4 };

print( "Lua:g1_st:" , st1, st2, st3, st4 )

Rg1 = 6.0
Rg2 = 3.0
Rg3 = 2.0

--         from to  Up nseg,mseg   width, hegith
gFw = Girder( n0, nFw, xvec, 80, 2, {Rg1,Rg1}, "Steel", g1_st )
gBk = Girder( n0, nBk, xvec, 30, 2, {Rg1,Rg1}, "Steel", g1_st )
gDw = Girder( n0, nDw, xvec, 20, 2, {Rg1,Rg1}, "Steel", g1_st )
gUp = Girder( n0, nUp, xvec, 20, 2, {Rg1,Rg1}, "Steel", g1_st )

gLf  = Girder( n0, nLf, zvec, 10, 2, {Rg2,Rg2}, "Steel", g1_st )
gRt  = Girder( n0, nRt, zvec, 10, 2, {Rg2,Rg2}, "Steel", g1_st )

gWdw = Girder( nDw, nFw, xvec, 70, 2, {Rg3,Rg3}, "Steel", g1_st )
gWup = Girder( nUp, nFw, xvec, 70, 2, {Rg3,Rg3}, "Steel", g1_st )

print( "#------ Lua:Ropes" )

-- ToDo: Ropes should be pre-strained (pre-tensioned) to avoid slack, we should set pre-strain force for each rope, the leght should be calculated from the rope material properties and pre-strain force


preS1 = 0.001
--          type     thick[mm]
Rope(nBk,nLf, 25.,4, preS1,  "Kevlar" );
Rope(nFw,nLf, 25.,4, preS1, "Kevlar" )
Rope(nBk,nRt, 25.,4, preS1, "Kevlar" );
Rope(nFw,nRt, 25.,4, preS1, "Kevlar" )
Rope(nBk,nDw, 25.,4, preS1, "Kevlar" );
--Rope(nFw,nDw, 25.,4, preS1, "Kevlar" )
Rope(nBk,nUp, 25.,4, preS1, "Kevlar");
--Rope(nFw,nUp, 25.,4, preS1, "Kevlar" )

Rope(nLf,nDw, 25.,4, preS1, "Kevlar");
Rope(nLf,nUp, 25.,4, preS1, "Kevlar");
Rope(nRt,nDw, 25.,4, preS1, "Kevlar");
Rope(nRt,nUp, 25.,4, preS1, "Kevlar");

preS2 = 0.003

Rope2( {gFw,gLf}, {0.5,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gRt}, {0.5,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gDw}, {0.5,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gUp}, {0.5,2.0}, 25.,4, preS2, "Kevlar");

Rope2( {gFw,gLf}, {0.25,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gRt}, {0.25,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gDw}, {0.25,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gUp}, {0.25,2.0}, 25.,4, preS2, "Kevlar");

Rope2( {gFw,gLf}, {0.75,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gRt}, {0.75,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gDw}, {0.75,2.0}, 25.,4, preS2, "Kevlar");
Rope2( {gFw,gUp}, {0.75,2.0}, 25.,4, preS2, "Kevlar");

print( "#------ Lua:BoundNodes and Girders" )

--- Finished basic layout of the ship ( main girders and ropes )

--nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {0.0,100.0,100.0} );
--nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {0.0,-100.0,0.0} );
--nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {100.0,0.0,0.0} );
--nB1 = BoundNode( gFw, ComponetKind.Girder, 0.5, {-100.0,0.0,0.0} );
--nB2 = BoundNode( gUp, ComponetKind.Girder, 0.5, {0.0,100.0,100.0} );
--gB  = Girder( nB1, nB2, xvec, 10, 2, {4.0,4.0}, "Steel", g1_st )

print( "#------ Lua:Radiataors" )

-- =RadiatorType = {"LithiumHeatPipe", 1280.0 }
--Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
--Radiator( gFw,0.15,0.8, gWdw,0.02,0.8, 1280.0 )
--Radiator( gFw,0.15,0.8, gWup,0.02,0.8, 1280.0 )

print( "#------ Lua:Rings & Sliders" )


-- r_pitch = Ring( {8.0,0.0,0.0}, xvec, yvec, 64, 80.0, {16.0,4.0}, "Steel", g1_st )

r_roll  = Ring2( {gLf,gUp,gDw,gRt}, {0.7,0.3,0.3,-1.0}, {0.0,0.0,1.0}, 64, { 4.0,4.0}, "Steel", g1_st, 0 )
-- r_roll  = Ring2( {gLf,gUp,gDw,-1}, {0.8,0.3,0.3,-1.0}, {0.0,0.0,10.0}, 64, { 4.0,4.0}, "Steel", g1_st, 0 )
r_pitch = Ring2( {gBk,gUp,gDw,gFw}, {0.2,0.4,0.4,-1.0}, {1.0,0.0,0.0}, 64, {16.0,4.0}, "Steel", g1_st, 1 )
-- r_pitch = Ring2( {gBk,gUp,gDw,-1}, {0.2,0.4,0.4,-1.0}, {1.0,0.0,0.0}, 64, {16.0,4.0}, "Steel", g1_st )

--      node1,2, up,  nseg    R     {width,height}
--r_roll  = Ring( {0.0,0.0,8.0}, zvec, xvec, 64, 55.0, {4.0,4.0}, "Steel", g1_st )
--r_pitch = Ring( {8.0,0.0,0.0}, xvec, yvec, 64, 80.0, {16.0,4.0}, "Steel", g1_st )
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
