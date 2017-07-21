
#ifndef CMesh_h
#define CMesh_h

#include "Vec2.h"
#include "Vec3.h"

// CMesh ... C like "constant" or C-language (vs. C++)
class CMesh{ public:
    int nvert ;
    int nedge ;
    int ntri  ;
    int nfaces;
    Vec3d * verts;
    Vec2i * edges;
    Vec3i * tris ;  // later we may do polygon faces ?
    int   * ngons;
    int   * faces;
};

#endif






