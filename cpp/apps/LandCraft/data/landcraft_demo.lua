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

-- Flooding/Draining example: set one inflow and one outflow then relax
inflow { ix=20, iy=20, delta=50.0 }
outflow{ ix=40, iy=40 }
for i=1,50 do relax_all{} end
print("[Lua] Flood/Drain relax done")

-- PathFinder: bind, set params, add centers, expand, find connections, make paths
pf_bind{}
pf_params{ ch2=1.0, chminus=0.0, chplus=0.0 }
pf_clear_centers{}
pf_add_center{ ix=10, iy=10 }
pf_add_center{ ix=60, iy=50 }
pf_prepare{}
for i=1,20 do pf_step{} end
local nconn = pf_find_connections{}
local npaths = pf_make_paths{}
print("[Lua] PF connections:", nconn, " paths:", npaths)
local paths = pf_paths{}
if #paths>0 then print("[Lua] First path length:", #paths[1]) end

-- Vehicles: create default type, spawn on road, step a few times, query state
local vtype = vehicle_type_default{}
local vid = vehicle_spawn{ road=rid, type=vtype }
for i=1,10 do vehicle_step_all{ dt=1.0 } end
local status, ipath, dir, onWay = vehicle_state{ id=vid }
print("[Lua] Vehicle:", status, ipath, dir, onWay)

-- Economy: load technologies, create factory, set tech, feed inputs, produce
-- Note: adjust file path if you have a specific technologies file
local ntech = econ_load{ file = "cpp/apps/LandCraft/data/Technologies.txt" }
print("[Lua] Technologies loaded:", ntech)
if ntech>0 then
    local fid = factory_new{}
    factory_set_tech{ id=fid, tech=0 }
    factory_set{ id=fid, commodity="ore", amount=100.0 }
    factory_set{ id=fid, commodity="fuel", amount=50.0 }
    local mass = factory_produce{ id=fid, N=1.0 }
    print("[Lua] Factory produced mass:", mass)
end
