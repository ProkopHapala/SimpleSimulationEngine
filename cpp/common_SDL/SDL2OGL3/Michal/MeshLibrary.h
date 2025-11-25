#ifndef MESH_LIBRARY_H
#define MESH_LIBRARY_H

#include "GLMesh.h"
#include "GLInstancedMesh.h"

namespace MeshLibrary {
    extern GLMesh<MPOS> point;
    extern GLMesh<MPOS> pointCross;
    
    extern GLMesh<MPOS> line;
    extern GLMesh<MPOS> line2D;

    extern GLMesh<MPOS> wireCube;
    extern GLMesh<MPOS,MNORMAL> cubeWithNormals;

    extern GLMeshBase<MPOS> sphere;
    extern GLInstancedMeshBase<GLvbo<MPOS>, MPOSOFFSET, MRADIUS, MCOLOR> sphereInstanced;
    extern GLMesh<MPOS> rect;
    extern GLMesh<MPOS> circle;
    extern GLMesh<MPOS> cross;
    extern GLMesh<MPOS> xmark;

    extern GLvbo<MPOS> octSphere;
    extern GLMesh<MPOS> octSphereMesh;
    extern GLInstancedMeshBase<GLvbo<MPOS>, MPOSOFFSET> octSphereInstanced;
}

#endif // MESH_LIBRARY_H
