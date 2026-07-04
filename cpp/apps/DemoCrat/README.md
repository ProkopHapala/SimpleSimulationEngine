# DemoCrat

Plugin-based demo framework using `dlopen`/`dlsym` to dynamically load and run visual demos as shared libraries (`.so`). Each plugin implements the `Demo` interface (`setup()`, `draw()`, `onMouse()`). The host application provides SDL2/OpenGL context, GUI, and hot-reloading of plugins.

## Files

- **DemoCrat_main.cpp** — host application: loads `.so` plugins, provides SDL2/OpenGL window, GUI for plugin selection, hot-reload support
- **Demo.h** — abstract `Demo` interface and C export function declarations for plugins
- **CMakeLists.txt** — build target: `DemoCrat_main`
- **data/** — demo plugin data files
