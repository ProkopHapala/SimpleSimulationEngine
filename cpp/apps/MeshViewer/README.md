# MeshViewer

Fast SDL2/OpenGL mesh visualizer with folder browser and animation playback.

## Features

- **Folder browser**: TreeView with single-click file loading. Mesh files (.obj, .obj+, .npz) highlighted
- **OBJ+ format**: standard OBJ with optional header `# nverts nedges nfaces` for single-pass loading
- **NPZ binary format**: packed binary with magic header + raw arrays for fastest loading
- **Animation**: multi-snapshot files (separated by `# --- snapshot N ---` in OBJ+, or concatenated in NPZ)
  - Each snapshot can have different vertex/edge/face counts (topology changes, like .xyz molecules)
  - Auto-play with Space, step with `[` / `]`, GUI slider for frame selection
- **Rendering**: solid (lit), wireframe (derived from triangles if no explicit edges), points — toggle with S/W/P or GUI checkboxes
- **Index labels**: billboarded vertex/edge/face numbers (keys 7/8/9 or GUI), constant screen size via Draw3D::drawInt
- **Face normals**: red arrows from face centroids (GUI checkbox)
- **Camera**: right-drag orbit (base class), wheel zoom, arrow keys rotate, numpad 0-6 preset views
- **SVG export**: faces/edges/points/normals/labels with backface culling, opacity, projection axis (key X)

## Build

```
cd cpp && mkdir -p build && cd build
cmake .. -DWITH_SDL=ON
make meshViewer
```

## Usage

```
./meshViewer [start_dir]
```

## File Formats

### OBJ+ (human-readable)
```
# 4 6 4          ← optional header: nverts nedges nfaces (enables single-pass load)
v 0 0 0
v 1 0 0
...
l 1 2           ← edges (1-indexed, same as OBJ)
...
f 1 2 3         ← faces (triangles)
# --- snapshot 1 ---   ← separator for multi-snapshot files
# 3 3 1
v ...
```

### NPZ (binary, packed)
```
magic[4] = "MNPZ"
nverts (int32), nedges (int32), nfaces (int32)
verts[nverts*3] (float64)
edges[nedges*2] (int32)
tris[nfaces*3] (int32)
[repeat for multi-snapshot]
```

## Controls

| Key | Browse mode | View mode |
|-----|-------------|-----------|
| Up/Down | Navigate files | Rotate camera (pitch) |
| Left/Right | — | Rotate camera (yaw) |
| Enter | Open dir/file | — |
| Backspace | Parent dir | — |
| Esc | Quit | Quit |
| Space | — | Play/pause animation |
| `[` / `]` | — | Prev/next frame |
| W/S/P | — | Toggle wire/solid/points |
| 7/8/9 | — | Toggle vert/edge/face numbers |
| F | — | Fit camera to mesh |
| R | — | Reload file |
| Numpad 0-6 | — | Preset views (front/back/top/bottom/left/right/default) |
| X | — | Export SVG |
| 1-6 | — | SVG options (faces/opacity/backface/edges/axis/points) |
| N/M | — | SVG vert/edge number labels |
| B | — | SVG normals |
| Right-drag | — | Orbit camera (base class) |
| Wheel | — | Zoom |
| Click (tree) | Load file / open dir | — |
