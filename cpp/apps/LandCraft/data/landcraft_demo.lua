-- Minimal demo sequence for LandCraft Lua API
-- Assumes the app is built with -DWITH_LUA=ON and LandCraft registers its API on start.

-- Generate terrain
generate_terrain{ seed=16464, maxHeight=500.0 }

-- Rain accumulation and river extraction
local wmax = gather_rain{ minFlow=100.0 }
local nR = find_rivers{ minFlow=50.0 }
print("[Lua] Rivers found:", nR, " wmax:", wmax)

-- Build a straight road across interesting terrain
local rid = build_road{ ax=10, ay=15, bx=55, by=38 }
print("[Lua] Road id:", rid)

-- Query profile along the road
local g, w = road_profile{ id=rid }
print("[Lua] Road profile samples:", #g)

-- Optional: dump first few samples
for i=1,math.min(5,#g) do
    print(string.format("[Lua] i=%d ground=%.3f water=%.3f", i, g[i], w[i]))
end
