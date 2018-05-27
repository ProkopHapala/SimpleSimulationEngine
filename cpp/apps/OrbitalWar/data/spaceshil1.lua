--- https://www.lua.org/pil/5.3.html

print("BEGIN spaceship");

-- http://www.mse.mtu.edu/~drjohn/my4150/props.html
Material{ name="Kevlar", density=1.44e+3, Spull=3.6e+9, Spush=0.0, Kpull=154.0e+9, Kpush=0.0, reflectivity=0.6,  Tmelt=350 }
Material{ name="Steel" , density=7.89e+3, Spull=1.2e+9, Spush=0.0, Kpull=200.0e+9, Kpush=0.0, reflectivity=0.85, Tmelt=800 }
--Material{ name="Titanium" , density=7.89e+3, Spull=3.6e+9, Spush=0.0, Kpull=154.0e+9, Kpush=0.0, reflectivity=0.7, Tmelt=450 }

--- newSpaceShip( "Ship1" )

origin = {0.0,0.0,0.0}
xvec   = {1.0,0.0,0.0} 
yvec   = {0.0,1.0,0.0}
zvec   = {0.0,0.0,1.0}

n0 = Node( origin );
n1 = Node( {-100.0,   0.0,    0.0} )
n2 = Node( { 100.0,   0.0,    0.0} )
n3 = Node( {   0.0,-200.0,    0.0} )
n4 = Node( {   0.0, 200.0,    0.0} )
n5 = Node( {   0.0,   0.0, -300.0} )
n6 = Node( {   0.0,   0.0,  800.0} )

print( "Lua:Nodes:" , n0,n1,n2,n3,n4,n5,n6 )

--         from to  Up nseg,mseg   width, hegith 
g1 = Girder( n0, n1, zvec, 10, 2, {10.0,8.0}, "Steel" )
g2 = Girder( n0, n2, zvec, 10, 2, {10.0,8.0}, "Steel" )
g3 = Girder( n0, n3, xvec, 20, 2, {10.0,8.0}, "Steel" )
g4 = Girder( n0, n4, xvec, 20, 2, {10.0,8.0}, "Steel" )
g5 = Girder( n0, n5, xvec, 30, 2, {10.0,8.0}, "Steel" )
g6 = Girder( n0, n6, xvec, 80, 2, {10.0,8.0}, "Steel" )

g7 = Girder( n3, n6, xvec, 70, 2, {10.0,8.0}, "Steel" )
g8 = Girder( n4, n6, xvec, 70, 2, {10.0,8.0}, "Steel" )

--          type     thick[mm]
Rope(n5,n1, 25, "Kevlar"); 
Rope(n6,n1, 25, "Kevlar" )
Rope(n5,n2, 25, "Kevlar"); 
Rope(n6,n2, 25, "Kevlar" )
Rope(n5,n3, 25, "Kevlar"); 
--Rope(n6,n3, 25, "Kevlar" )
Rope(n5,n4, 25, "Kevlar"); 
--Rope(n6,n4, 25, "Kevlar" )

Rope(n1,n3, 25, "Kevlar");
Rope(n1,n4, 25, "Kevlar");
Rope(n2,n3, 25, "Kevlar");
Rope(n2,n4, 25, "Kevlar");


-- =RadiatorType = {"LithiumHeatPipe", 1280.0 }

--Radiator( g5,0.2,0.8, g1,0.2,0.8, 1280.0 )
Radiator( g6,0.2,0.8, g7,0.1,0.8, 1280.0 )
Radiator( g6,0.2,0.8, g8,0.1,0.8, 1280.0 )



Tank( {16,16,16}, zvec, {50.0,10.0}, "H2");



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
