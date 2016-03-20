    #ifndef  Solids_h
    #define  Solids_h

    #include "Vec3.h"

    class Solid{

        int nvert;
        int nedge;
        int ntri;

        Vec3d * verts;
        int   * edges;
        Vec3i * tris;  // later we may do polygon faces ?

    };

    namespace Solids{

/*
        struct Tetrahedron{
            const static int    nVerts = 4;
            const static int    nEdges = 6;
            const static int    nTris  = 4;
            constexpr static Vec3d verts [nVerts]    = { {-1.0d,-1.0d,-1.0d}, {+1.0d,+1.0d,-1.0d}, {-1.0d,+1.0d,+1.0d}, {+1.0d,-1.0d,+1.0d} };
            constexpr static int   edges [nEdges][2] = { {0,1},{0,2},{0,3}, {1,2},{1,3},{2,3} };
            constexpr static int   tris  [nTris ][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3} };
        } tetrahedron;
*/

        // ==== Tetrahedron

        const static int Tetrahedron_nverts = 4;
        const static int Tetrahedron_nedges = 6;
        const static int Tetrahedron_ntris  = 4;
        const static int Tetrahedron_nfaces = 4;
        static Vec3d     Tetrahedron_verts   [Tetrahedron_nverts]   = { {-1.0d,-1.0d,-1.0d}, {+1.0d,+1.0d,-1.0d}, {-1.0d,+1.0d,+1.0d}, {+1.0d,-1.0d,+1.0d} };
        static int       Tetrahedron_edges   [Tetrahedron_nedges*2] = { 0,1,   0,2,   0,3,    1,2,  1,3,  2,3 };
        static int       Tetrahedron_tris    [Tetrahedron_ntris *3] = { 0,2,1, 0,1,3, 0,3,2,  1,2,3 };
        static int       Tetrahedron_ngons   [Tetrahedron_nfaces  ] = { 3,     3,     3,      3     };
        static int       Tetrahedron_faces   [Tetrahedron_nfaces*3] = { 0,2,1, 0,1,3, 0,3,2,  1,2,3 }; // order is important for normals

        // ==== Octahedron

        const static int Octahedron_nverts = 6;
        const static int Octahedron_nedges = 12;
        const static int Octahedron_nfaces = 8;
        const static int Octahedron_ntris  = 8;
        static Vec3d     Octahedron_verts   [Octahedron_nverts]   = { {-1.0d,0.0d,0.0d}, {+1.0d,0.0d,0.0d}, {0.0d,-1.0d,0.0d}, {0.0d,+1.0d,0.0d}, {0.0d,0.0d,-1.0d}, {0.0d,0.0d,+1.0d} };
        static int       Octahedron_edges   [Octahedron_nedges*2] = { 0,2, 0,3, 0,4, 0,5,   1,2, 1,3, 1,4, 1,5,   2,4, 2,5, 3,4, 3,5  };
        static int       Octahedron_ngons   [Octahedron_nfaces  ] = { 3,     3,     3,     3,       3,     3,     3,     3     };
        static int       Octahedron_faces   [Octahedron_nfaces*3] = { 0,4,2, 0,2,5, 0,3,4, 0,5,3,   1,2,4, 1,5,2, 1,4,3, 1,3,5 };

        // ==== Cube

        const static int  Cube_nverts = 8;
        const static int  Cube_nedges = 12;
        //const static int  nCube_tris  = 6;
        const static int  Cube_nfaces  = 6;
        static Vec3d      Cube_verts   [Cube_nverts]    = {
            {-1.0d,-1.0d,-1.0d},
            {-1.0d,-1.0d,+1.0d},
            {-1.0d,+1.0d,-1.0d},
            {-1.0d,+1.0d,+1.0d},
            {+1.0d,-1.0d,-1.0d},
            {+1.0d,-1.0d,+1.0d},
            {+1.0d,+1.0d,-1.0d},
            {+1.0d,+1.0d,+1.0d},
        };
        static int     Cube_edges [Cube_nedges*2] = { 0,1, 0,2, 0,4,  1,3,1,5, 2,3, 2,6,    7,5, 7,6, 7,3,  5,4, 6,4    };
        //static int   Cube_tris  [nCube_tris ][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3} };
        static int     Cube_ngons [Cube_nfaces  ] = { 4,        4,          4,         4,         4,         4        };
        static int     Cube_faces [Cube_nfaces*4] = { 0,1,3,2,  0,4,5,1,    0,2,6,4,   7,5,4,6,   7,3,1,5,   7,6,2,3  };


    };

    #endif﻿
