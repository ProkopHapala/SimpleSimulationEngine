/// @file spaceCraftMeshExport.cpp
/// @brief Headless Lua ship → mesh export CLI — sketch or blocks LOD, optional OBJ/SVG/truss/tags outputs.
///
/// Loads a ship script (`-s`), builds **SpaceCraft** via **EditSpaceCraft** bindings, then:
/// - `-lod sketch` → **BuildCraft_sketch** (wireframe; no `.truss`)
/// - `-lod blocks` → **BuildCraft_blocks** from **SpaceCraft2Mesh_blocks.h** (FEM truss; default)
///
/// Sets `g_luaScriptDir` from the `-s` path so `fromObj("hull.obj")` resolves relative to the script folder.
/// `-svg` writes a 4-panel wireframe SVG; `-g` writes tagged sketch topology OBJ.
///
/// Open issues / caveats:
/// - **Uses SpaceCraft2Mesh_blocks.h**, not **SpaceCraft2Mesh2.h** — editor/dynamics apps may produce different
///   truss topology for the same Lua ship (duplicate **BuildCraft_blocks**; do not unify without review).
/// - **`-lod blocks` here omits sliders, shields, radiators, welds** — only nodes/girders/ropes/rings truss.
/// - **`-lod sketch` + `-t`**: truss export skipped with warning — sketch is not a simulation mesh.
/// - **Deferred girders/ropes** warn-skipped in blocks LOD — silent partial ships if `fulfillTag` forgotten.
/// - **`-svg` is edges-only** (Builder2::exportSVG sets `ntri=0`) — plate triangles from sketch LOD invisible.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "globals.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "LuaHelpers.h"
#include "IO_utils.h"
#include "argparse.h"

#include "MeshBuilder2.h"
#include "SpaceCraft.h"
#include "SpaceCraft2Mesh_blocks.h"
#include "SpaceCraftFromOBJ.h"
#include "EditSpaceCraft.h"

using namespace Mesh;
using namespace SpaceCrafting;

static void setLuaScriptDirFromPath( const char* shipPath ){
    char buf[512];
    strncpy(buf, shipPath, sizeof(buf)-1); buf[sizeof(buf)-1]=0;
    char* slash = strrchr(buf, '/');
    if(slash){ *slash=0; strncpy(g_luaScriptDir, buf, sizeof(g_luaScriptDir)-1); g_luaScriptDir[sizeof(g_luaScriptDir)-1]=0; }
    else { g_luaScriptDir[0]=0; }
}

static void print_usage(const char* prog){
    printf("Usage: %s -s <ship.lua> [-o <out.obj>] [-svg <out.svg>] [-t <out.truss>] [-lod sketch|blocks] [-g <tags.obj>] [-v <verbosity>]\n", prog);
}

static int reloadShip_noSim(const char* fname, Builder2& mesh, bool bSketchLOD){
    if(verbosity>0) printf("#### START reloadShip_noSim('%s') lod=%s\n", fname, bSketchLOD?"sketch":"blocks");
    if(!theSpaceCraft){
        theSpaceCraft = new SpaceCraft();
    }
    theSpaceCraft->clear();
    setLuaScriptDirFromPath(fname);
    if( Lua::dofile(theLua,fname) ){
        printf("ERROR in reloadShip_noSim() Lua::dofile(%s) failed\n", fname);
        return 1;
    }
    theSpaceCraft->checkIntegrity();
    mesh.clear();
    if(bSketchLOD) BuildCraft_sketch( mesh, *theSpaceCraft );
    else           BuildCraft_blocks( mesh, *theSpaceCraft, 30.0 );
    if(verbosity>0){
        mesh.printSizes();
        if(!bSketchLOD) checkSpaceCraftMesh(mesh, *theSpaceCraft, true, false);
    }
    if(verbosity>0) printf("#### END reloadShip_noSim('%s')\n", fname);
    return 0;
}

int main(int argc, char *argv[]){
    setbuf(stdout, NULL);

    const char* shipPath = nullptr;
    const char* objPath  = "ship.obj";
    const char* svgPath  = nullptr;
    const char* trussPath = nullptr;
    const char* tagsPath = nullptr;
    bool bSketchLOD = false;

    // Initialize Lua/SpaceCrafting
    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();

    LambdaDict funcs;
    funcs["-s"] = {1, [&](const char** ss){ shipPath  = ss[0]; if(verbosity>0) printf("ARG -s %s\n", shipPath); }};
    funcs["-o"] = {1, [&](const char** ss){ objPath   = ss[0]; if(verbosity>0) printf("ARG -o %s\n", objPath);  }};
    funcs["-svg"] = {1, [&](const char** ss){ svgPath   = ss[0]; if(verbosity>0) printf("ARG -svg %s\n", svgPath); }};
    funcs["-t"] = {1, [&](const char** ss){ trussPath = ss[0]; if(verbosity>0) printf("ARG -t %s\n", trussPath); }};
    funcs["-g"] = {1, [&](const char** ss){ tagsPath  = ss[0]; if(verbosity>0) printf("ARG -g %s\n", tagsPath);  }};
    funcs["-lod"] = {1, [&](const char** ss){
        if(strcmp(ss[0],"sketch")==0) bSketchLOD = true;
        else if(strcmp(ss[0],"blocks")==0) bSketchLOD = false;
        else { printf("ERROR: -lod must be 'sketch' or 'blocks' (got '%s')\n", ss[0]); exit(1); }
        if(verbosity>0) printf("ARG -lod %s\n", ss[0]);
    }};
    funcs["-v"] = {1, [&](const char** ss){ int v=1; sscanf(ss[0], "%d", &v); verbosity = v; printf("ARG -v %d\n", verbosity); }};

    process_args(argc, argv, funcs);

    if(!shipPath){
        print_usage(argv[0]);
        return 1;
    }

    if(verbosity>0){
        printf("spaceCraftMeshExport: ship='%s' obj='%s' truss='%s' lod=%s verbosity=%d\n",
            shipPath, objPath, trussPath ? trussPath : "(none)", bSketchLOD?"sketch":"blocks", verbosity);
    }

    Builder2 mesh;
    int err = reloadShip_noSim(shipPath, mesh, bSketchLOD);
    if(err){ return err; }

    if(verbosity>0){
        printf("Writing OBJ '%s' ...\n", objPath);
    }
    mesh.write_obj(objPath);
    if(verbosity>0){
        printf("DONE: OBJ '%s' written (nvert=%zu, nedge=%zu, ntri=%zu)\n", objPath, mesh.verts.size(), mesh.edges.size(), mesh.tris.size());
    }

    if(svgPath){
        if(verbosity>0) printf("Writing SVG '%s' (4 views) ...\n", svgPath);
        Mat3d rots[4] = { Mat3dIdentity, Mat3d{0,0,1, 0,1,0, -1,0,0}, Mat3d{1,0,0, 0,0,-1, 0,1,0}, Mat3d{0.7071,0,0.7071, 0,1,0, -0.7071,0,0.7071} };
        mesh.exportSVGmultiView(svgPath, 4, rots, 2);
        if(verbosity>0) printf("DONE: SVG '%s' written\n", svgPath);
    }

    if(tagsPath){
        writeSpaceCraftSketchOBJ(tagsPath, *theSpaceCraft);
    }

    if(trussPath){
        if(bSketchLOD){
            printf("WARNING: sketch LOD — .truss export skipped (no truss dynamics on sketch)\n");
        }else{
            if(verbosity>0){
                printf("Writing TRUSS '%s' ...\n", trussPath);
            }
            exportSimToFile(trussPath, mesh, workshop);
            if(verbosity>0){
                printf("DONE: TRUSS '%s' written\n", trussPath);
            }
        }
    }

    return 0;
}
