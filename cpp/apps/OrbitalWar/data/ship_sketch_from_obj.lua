-- ship_sketch_from_obj.lua — minimal fromObj + fulfillTag smoke test (flat sketch_box.obj)
Material({ name="Steel", density=7800, Spull=5e8, Spush=5e8, Kpull=2e11, Kpush=2e11, reflectivity=0.3, Tmelt=1800 })
StickMaterial("GS1_long", "Steel", 50, 2, 0.0)

local grp = fromObj("sketch_box.obj", { scale=10.0, pos={0,0,0} })
print("fromObj nodes=", table.concat(grp.nodes, ","))
printDeferred()

fulfillTag("girder:GS1_long", { up={0,0,1}, nseg=40, mseg=4, wh={0.05,0.05}, st={0,0,0,0} })

exportSketchOBJ("sketch_box_tags.obj")
