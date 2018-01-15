### Shaders

- make library of basic shaders in .h file file (easier code managemetn and name searching)

### Meshes

- make library of basic shapes in .h file

﻿### Terrain

 - Sampling height map in vertex shader
 - Vertex array composed of trinagles
 - Camera space / world space ?
    - Artifacts when camera move ?
 - Logaritmic depth buffer

﻿### Instances

﻿### Trees

 - interpolate texture atlas depending on ray angle from camera
 - Generate quads in geometry shader
    - shader can automatically do frustrum culling
 - signed distance in alpha channel to enable low-res texture atlas
    - texture atlas generated from hi-res renders by frag shader

### clouds

 - generate 3D texture by turbulent noise-wraping
