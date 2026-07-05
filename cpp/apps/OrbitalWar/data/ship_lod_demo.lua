-- ship_lod_demo.lua — low-res OBJ (sketch_prism.obj) → fulfillTag → sketch + blocks export
-- Test: tests_bash/Orbital/spaceCraftLOD_demo.sh → review_lod_demo/*.obj + *.svg
-- Caveats: fulfillTag is girder-only; shields/radiators not in CLI blocks mesh; no sliders/rings in OBJ import.
Material({ name="Steel", density=7800, Spull=5e8, Spush=5e8, Kpull=2e11, Kpush=2e11, reflectivity=0.3, Tmelt=1800 })
StickMaterial("GS1_long", "Steel", 50, 2, 0.0)

local grp = fromObj("sketch_prism.obj", { scale=5.0, pos={0,0,0} })
print(string.format("fromObj: %i nodes, %i girders, %i shields, %i radiators",
    #grp.nodes, #grp.girders, #grp.shields, #grp.radiators))
printDeferred()

-- Fulfill all girders tagged GS1_long (workshop stick type from usemtl)
fulfillTag("girder:GS1_long", { up={0,0,1}, nseg=24, mseg=4, wh={0.05,0.05}, st={0,0,0,0} })
