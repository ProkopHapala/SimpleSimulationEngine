# MusicVizualizer

Real-time music visualizer using SDL2_mixer for audio playback and OpenGL 3+ shaders for visualization. Features frequency spectrum analysis, various shader-based visual effects, and planned reaction-diffusion patterns.

## Planned Features (from source comments)

- Reaction-diffusion (Belousov-Zhabotinsky, Turing patterns)
- Fluid dynamics multi-pass rendering
- Kaleidoscope effects
- Body-of-rotation render for spectrum
- L-system branching driven by music spectrum
- Diffusion random walks with music spectrum
- Particle flow harmonics

## Files

- **MusicVisualizer_main.cpp** — main application: audio playback, spectrum analysis, shader visualization
- **MusicRendererOGL3.h** — OpenGL 3+ renderer: shader-based visual effects, spectrum rendering
- **MusicUtils.h** — audio utilities: FFT, spectrum extraction, beat detection
- **CMakeLists.txt** — build target: `MusicVisualizer_main` (requires SDL2_mixer, `WITH_MUSIC=ON`)
