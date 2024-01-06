--- https://www.lua.org/pil/5.3.html

require( "data/lua/utils" )

print( "#------ Lua:StickMaterials START" )

--- long, perp, zigIn, zigOut
st1  = StickMaterial( "GS1_long", "Steel", 100, 5, 0.0 )
st2  = StickMaterial( "GS1_perp", "Steel", 50,  3, 0.0 )
st3  = StickMaterial( "GS1_in",   "Steel", 40,  2, 0.0 )
st4  = StickMaterial( "GS1_out",  "Steel", 40,  2, 0.0 )
-- stK1 = StickMaterial( "RK_1","Kevalr", 0.01, 0.001, 0.0 )

print( "#------ Lua:StickMaterials DONE" )
