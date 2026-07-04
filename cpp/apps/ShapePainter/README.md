# ShapePainter

Experimental drawing program where the brush creates useful geometric shapes rather than simple strokes. A quick hybrid between pixel and vector graphics — each brush stroke defines a parametric shape (corners, arcs, polygons) that can be edited after placement.

## Files

- **shapePainter_main.cpp** — main application: canvas, brush selection, shape editing
- **ShapePainter.h** — brush system: `Brush` base class, `CornerBrush` (parametric corner shape), shape manipulation
- **CMakeLists.txt** — build target: `ShapePainter_main`
