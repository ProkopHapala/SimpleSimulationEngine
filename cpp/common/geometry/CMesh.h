
#ifndef CMesh_h
#define CMesh_h

/// @file CMesh.h
/// @brief Non-owning C-style mesh view — the bridge between Builder2's vector storage and C/OpenCL APIs.
///
/// CMesh exists because Builder2 uses std::vector (which owns and may reallocate its data),
/// but GPU upload and file I/O need raw pointers with explicit sizes. CMesh is a "view":
/// it borrows pointers from Builder2's vectors (or from a loaded file buffer) without copying.
///
/// The "C" in CMesh means "C-language style" — no methods, no ownership, no templates.
/// This makes it suitable for passing across the C/C++ boundary, to OpenCL kernels,
/// and for direct memcpy to/from disk. The caller is responsible for lifetime management.
///
/// Builder2::addCMesh() imports in the reverse direction: copies CMesh data into Builder2's
/// vectors, transforming from flat pointer arrays to managed std::vector storage.

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






