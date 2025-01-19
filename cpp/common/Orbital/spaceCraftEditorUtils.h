
#ifndef spaceCraftEditorUtils_h
#define spaceCraftEditorUtils_h

#include "Tree.h"

//#include "Truss.h"
//#include "TriangleRayTracer.h"
//#include "Radiosity.h"

#include "SpaceCraft.h"
#include "TrussDynamics_d.h"

#include "GUI.h"

#ifdef Truss_h

int makeTruss( Truss& truss ){
    //truss.girder1( (Vec3d){-5.0,0.0,0.0}, (Vec3d){5.0,0.0,0.0}, (Vec3d){0.0,1.0,0.0}, 5, 1.0 );
    Truss trussPlan;
    trussPlan.loadXYZ(  "data/octShip.xyz" );
    //trussPlan.affineTransform( (Mat3d){5.5,0.0,0.0, 0.0,5.5,0.0, 0.0,0.0,5.5}, false );
    trussPlan.affineTransform( (Mat3d){6,0.0,0.0, 0.0,6,0.0, 0.0,0.0,6}, false );
    GirderParams * gpar = new GirderParams [trussPlan.edges.size()];
    //Vec3d ups  = new Vec3d[trussPlan.edges.size()];
    Vec3d ups[] = {
        (Vec3d){0.0,-1.0,0.0},
        (Vec3d){0.0,+1.0,0.0},
        (Vec3d){0.0,0.0,-1.0},
        (Vec3d){0.0,0.0,+1.0},
        (Vec3d){-1.0,0.0,0.0},
        (Vec3d){+1.0,0.0,0.0}
    };
    //truss.makeGriders( 6, &trussPlan.edges[0], &trussPlan.points[0], gpar, ups );
    //truss.makeGriders( trussPlan, gpar, ups, NULL );
    truss.wheel( {0.0,0.0,0.0}, {10.0,0.0,0.0}, {0.0,1.0,0.0}, 50, 0.5 );
    std::vector<Vec2i> ends;
    truss.makeGriders( trussPlan, gpar, ups, &ends );
    truss.autoBridge(ends.size(), &ends[0], 0.8, 0 );
    truss.panel( {0.0,0.0,0.0}, {5.0,0.0,0.0}, {0.0,5.0,0.0}, {5.0,5.0,0.0}, {5,5}, 1.0 );
    delete [] gpar;

    int glo = glGenLists(1);
    glNewList( glo, GL_COMPILE );
    SpaceCrafting::drawTruss( truss, true );
    glEndList();

    return glo;
}

#endif // Truss_h

#endif
