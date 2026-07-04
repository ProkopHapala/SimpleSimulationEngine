
/// @file meshViewer_main.cpp
/// @brief Entry point for MeshViewer — SDL2/OpenGL mesh visualizer.

#include "MeshViewer.h"

MeshViewerApp* thisApp = nullptr;

int main(int argc, char* argv[]) {
    SDL_DisplayMode dm = initSDLOGL();
    int junk;
    thisApp = new MeshViewerApp(junk, dm.w - 150, dm.h - 100);

    // Start in current directory or argument-provided path
    std::string startDir = ".";
    if (argc > 1) startDir = argv[1];
    thisApp->setDir(startDir);

    thisApp->loop(1000000);
    return 0;
}
