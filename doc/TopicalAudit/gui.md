---
type: TopicalAudit
title: GUI Components
tags: [topic, cpp, python, sdl2, opengl, gui, widgets, panels, input, rendering]
---

## Summary

GUI framework spanning C++ (SDL2/OpenGL) and Python (moderngl/pygame). C++ provides a widget-based UI system with panels, sliders, buttons, text input, and multi-panel layouts integrated with the SDL2 OpenGL application framework. Python provides `BaseGUI.py` as a lightweight base class for GUI applications, with specialized variants for GLCL and pySymGLSL. The C++ framework is the primary, feature-complete implementation.

## Implementations

| Language | Location | Status | Notes |
|----------|----------|--------|-------|
| C++ | `cpp/common_SDL/SDL2OGL/GUI.h` / `GUI.cpp` | active | Core GUI framework: `GUITextInput`, `GUIAbstractPanel`, `GUIPanel`, `MultiPanel`, `BoundGUI`, `GUI`. Widgets: sliders, buttons, text input, panels. Event handling, rendering, focus management. |
| C++ | `cpp/common_SDL/SDL2OGL/AppSDL2OGL.h` / `.cpp` | active | Base SDL2/OpenGL application: window creation, event loop, input handling, child windows, timing. |
| C++ | `cpp/common_SDL/SDL2OGL/AppSDL2OGL_3D.h` / `.cpp` | active | 3D extension: quaternion camera, trackball rotation, mouse/keyboard camera control, first-person mode toggle. |
| C++ | `cpp/common_SDL/SDL2OGL/Draw2D.h` | active | 2D drawing primitives for UI rendering. |
| C++ | `cpp/common_SDL/SDL2OGL/Draw3D.h` | active | 3D drawing primitives for scene rendering. |
| C++ | `cpp/common_SDL/SDL2OGL3/` | active | OpenGL 3+ specific rendering utilities. |
| Python | `python/BaseGUI.py` | active | Lightweight base class for Python GUI applications. Provides window setup, event loop skeleton, draw/update callbacks. |
| Python | `python/GLCL/GLGUI.py` | active | GUI for GLCL (OpenCL+OpenGL) applications. Extends BaseGUI. |
| Python | `python/pySymGLSL/GLSL_GUI.py` | active | GUI for pySymGLSL: playback control, viewport resize, shader graph rebuild, uniform updates. Extends BaseGUI. |
| Python | `python/GUI/GUITemplate.py` | active | Template for Python GUI applications. |
| Doc | `doc/Markdown/GUI.md` | doc | Auto-generated API documentation for C++ GUI classes. |

## Sub-topics

### C++ GUI Widget System

`GUI.h` defines the widget hierarchy:
- **`GUITextInput`**: Text editing with cursor, selection, keyboard input
- **`GUIAbstractPanel`**: Base panel with position, size, title bar, dragging
- **`GUIPanel`**: Container panel holding widgets (sliders, buttons, text)
- **`MultiPanel`**: Tabbed/stacked multi-panel container
- **`BoundGUI`**: Panel with data bindings
- **`GUI`**: Top-level manager — focus handling, event dispatch, drawing

### C++ Application Framework

- **`AppSDL2OGL`**: SDL2 window + OpenGL context, event loop, input state (keyboard/mouse), child window management, FPS timing. Virtual methods: `onEvent`, `onKey`, `onMouse`, `draw`, `update`
- **`AppSDL2OGL_3D`**: Adds 3D camera with quaternion rotation. Trackball-style mouse drag for rotation. Key controls: WASD movement, arrow keys for rotation, toggle perspective/first-person mode. `draw3D()` virtual method for scene rendering.

### Python GUI

- **`BaseGUI`**: Minimal base class — window creation (pygame or moderngl), event loop, `draw()` and `update()` virtual methods, `onKeyDown()` callback
- **`GLSL_GUI`**: Integrates with `GLSL_Simulation` — shader pipeline parsing, recompilation from disk, uniform value updates, playback start/stop, viewport resize
- **`GLGUI`**: OpenCL+OpenGL specific GUI integration

## Parity Status

- **C++ GUI ↔ Python BaseGUI**: Different paradigms. C++ has full widget system (panels, sliders, buttons, text input). Python has minimal base class with virtual callbacks. No widget library in Python.
- **C++ AppSDL2OGL_3D ↔ Python**: No Python equivalent of 3D camera framework. Python GUIs are 2D-focused or use moderngl context directly.
- No formal parity needed — C++ and Python GUIs serve different use cases (interactive 3D apps vs. lightweight simulation viewers).

## Open Issues

- C++ GUI has no layout manager — panels positioned manually or via hardcoded coordinates
- No C++ widget for color picker, file dialog, or dropdown menu
- Python GUI lacks widget library — each application reimplements UI controls
- `GUI.cpp` text input handling is complex with potential edge cases in cursor/selection logic
- No theming/styling system — colors and fonts hardcoded
- `AppSDL2OGL_3D` camera control is basic — no orbit target, no zoom-to-fit
- No automated tests for GUI event handling or widget behavior
