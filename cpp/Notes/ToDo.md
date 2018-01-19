
## Maps

- Do path finding
- All functionality of SimplexGrid souch as hydraulic errosion, flooding and raymarching should be made also for SquareTerrain
- Make order in different types of terrains, rulers and related datastructures.
  - Each terrain or other datastructure should use ruler (not it's own ruler-like methods)
  - path-finding and hydraulic errosion should be shared on ruler level
  - Do ve need specialized Terran hydraulicas and path-finding class


### Maps unifications

- We have following typs of Rulers
    - 2D:
        SimplexRuler.h  (double Step)
        Ruler2DFast.h - (Vec2 Step)
        SquareRuler.h   (double step)
        Map2D.h         (double step)
    - 3D:
        CubicRuler.h  (Mat3 Step)
        CubeGridRuler (double Step)  (in grids3D.h    )
        GridShape     (Mat3d   cell) (in dataStructures)
- Grids:
    - Templated Oversophisticated
        Grid3D.h - Template const size
        GridMap2D.h GridMap2D_stub - Templates (double step)
        SimplexGrid.h templates
    - 3D:
        TileBuffer3D.h   (Templated)
    - 2D:
        HashMap2D.h
        TerrainCubic : public Map2D 
        TerrainHydraulics  (not simplex) ... agnostic to ruler but with hexagonal neighborhoos 
        TerrainSimplex    ( copied hydraulics, + simplex ruler )
        TerrainSquare     ( copied hydraulics, contains SquareRuler )
        TileBuffer2D.h   (Templated)
