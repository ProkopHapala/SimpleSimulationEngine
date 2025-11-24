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
#include "EditSpaceCraft.h"

using namespace Mesh;
using namespace SpaceCrafting;

static void print_usage(const char* prog){
    printf("Usage: %s -s <ship.lua> [-o <out.obj>] [-t <out.truss>] [-v <verbosity>]\n", prog);
}

static int reloadShip_noSim(const char* fname, Builder2& mesh){
    if(verbosity>0) printf("#### START reloadShip_noSim('%s')\n", fname);
    if(!theSpaceCraft){
        theSpaceCraft = new SpaceCraft();
    }
    theSpaceCraft->clear();
    if( Lua::dofile(theLua,fname) ){
        printf("ERROR in reloadShip_noSim() Lua::dofile(%s) failed\n", fname);
        return 1;
    }
    theSpaceCraft->checkIntegrity();
    mesh.clear();
    BuildCraft_blocks( mesh, *theSpaceCraft, 30.0 );
    if(verbosity>0){
        mesh.printSizes();
        checkSpaceCraftMesh(mesh, *theSpaceCraft, true, false);
    }
    if(verbosity>0) printf("#### END reloadShip_noSim('%s')\n", fname);
    return 0;
}

int main(int argc, char *argv[]){
    setbuf(stdout, NULL);

    const char* shipPath = nullptr;
    const char* objPath  = "ship.obj";
    const char* trussPath = nullptr;

    // Initialize Lua/SpaceCrafting
    theSpaceCraft = new SpaceCraft();
    initSpaceCraftingLua();

    LambdaDict funcs;
    funcs["-s"] = {1, [&](const char** ss){ shipPath  = ss[0]; if(verbosity>0) printf("ARG -s %s\n", shipPath); }};
    funcs["-o"] = {1, [&](const char** ss){ objPath   = ss[0]; if(verbosity>0) printf("ARG -o %s\n", objPath);  }};
    funcs["-t"] = {1, [&](const char** ss){ trussPath = ss[0]; if(verbosity>0) printf("ARG -t %s\n", trussPath); }};
    funcs["-v"] = {1, [&](const char** ss){ int v=1; sscanf(ss[0], "%d", &v); verbosity = v; printf("ARG -v %d\n", verbosity); }};

    process_args(argc, argv, funcs);

    if(!shipPath){
        print_usage(argv[0]);
        return 1;
    }

    if(verbosity>0){
        printf("spaceCraftMeshExport: ship='%s' obj='%s' truss='%s' verbosity=%d\n", shipPath, objPath, trussPath ? trussPath : "(none)", verbosity);
    }

    Builder2 mesh;
    int err = reloadShip_noSim(shipPath, mesh);
    if(err){ return err; }

    if(verbosity>0){
        printf("Writing OBJ '%s' ...\n", objPath);
    }
    mesh.write_obj(objPath);
    if(verbosity>0){
        printf("DONE: OBJ '%s' written (nvert=%zu, nedge=%zu, ntri=%zu)\n", objPath, mesh.verts.size(), mesh.edges.size(), mesh.tris.size());
    }

    if(trussPath){
        if(verbosity>0){
            printf("Writing TRUSS '%s' ...\n", trussPath);
        }
        exportSimToFile(trussPath, mesh, workshop);
        if(verbosity>0){
            printf("DONE: TRUSS '%s' written\n", trussPath);
        }
    }

    return 0;
}
