
/// @file meshViewer_main.cpp
/// @brief Entry point for MeshViewer — SDL2/OpenGL mesh visualizer.

#include "MeshViewer.h"
#include "argparse.h"

MeshViewerApp* thisApp = nullptr;

int main(int argc, char* argv[]) {
    SDL_DisplayMode dm = initSDLOGL();
    int junk;
    thisApp = new MeshViewerApp(junk, dm.w - 150, dm.h - 100);

    LambdaDict funcs;

    funcs["-dir"] = {1, [&](const char** ss){
        thisApp->setDir(ss[0], true);
    }};

    funcs["-file"] = {1, [&](const char** ss){
        std::string path = ss[0];
        size_t slash = path.find_last_of('/');
        std::string dir = (slash != std::string::npos) ? path.substr(0, slash) : ".";
        thisApp->setDir(dir, false);
        thisApp->loadFile(path);
    }};

    thisApp->setDir(".", true);
    process_args(argc, argv, funcs, false);

    thisApp->loop(1000000);
    return 0;
}
