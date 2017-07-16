#ifndef  Solids_h
#define  Solids_h

#include "fastmath.h"
#include "Vec3.h"

    class Solid{ public:

        int nvert;
        int nedge;
        int ntri;

        Vec3d * verts;
        int   * edges;
        int   * tris;  // later we may do polygon faces ?

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

        //const static Solid Tetrahedron = (Solid){Tetrahedron_nverts,Tetrahedron_nedges,Tetrahedron_ntris, NULL,NULL,NULL};
        //const static Solid Tetrahedron = (Solid){Tetrahedron_nverts,Tetrahedron_nedges,Tetrahedron_ntris, &Tetrahedron_verts,&Tetrahedron_edges,&Tetrahedron_tris};

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

        // ==== RhombicDodecahedron

        const static int  RhombicDodecahedron_nverts = 14;
        const static int  RhombicDodecahedron_nedges = 24;
        const static int  RhombicDodecahedron_nfaces = 12;
        static Vec3d      RhombicDodecahedron_verts [RhombicDodecahedron_nverts] = {
            {-1.0d,-1.0d,-1.0d},
            {-1.0d,-1.0d,+1.0d},
            {-1.0d,+1.0d,-1.0d},
            {-1.0d,+1.0d,+1.0d},
            {+1.0d,-1.0d,-1.0d},
            {+1.0d,-1.0d,+1.0d},
            {+1.0d,+1.0d,-1.0d},
            {+1.0d,+1.0d,+1.0d},

            {-2.0d,0.0d,0.0d},
            {+2.0d,0.0d,0.0d},
            {0.0d,-2.0d,0.0d},
            {0.0d,+2.0d,0.0d},
            {0.0d,0.0d,-2.0d},
            {0.0d,0.0d,+2.0d}
        };
        static int     RhombicDodecahedron_edges [RhombicDodecahedron_nedges*2] = { 0,8, 0,10, 0,12,   1,8, 1,10, 1,13,   2,8, 2,11, 2,12,  3,8, 3,11, 3,13,  4,9, 4,10, 4,12,  5,9, 5,10, 5,13,  6,9, 6,11, 6,12,   7,9, 7,11, 7,13   };
        static int     RhombicDodecahedron_ngons [RhombicDodecahedron_nfaces  ] = { 4,4,4,4,  4,4,4,4,    4,4,4,4  };
        static int     RhombicDodecahedron_faces [RhombicDodecahedron_nfaces*4] = { 0,10,1,8,  0,8,2,12,    0,12,4,10,   7,11,3,13,   7,9,6,11,   7,13,5,9,   1,13,3,8,  1,10,5,13,   4,9,5,10,  4,12,6,9,  2,8,3,11,  2,11,6,12    };

        // ==== Icosahedron

        const static int  Icosahedron_nverts = 12;
        const static int  Icosahedron_nedges = 30;
        const static int  Icosahedron_nfaces = 20;
        const double a = 1.0d;
        const double b = GOLDEN_RATIO;
        static Vec3d      Icosahedron_verts [RhombicDodecahedron_nverts] = {
            {0.0d,-a,-b}, // 0
            {0.0d,-a,+b}, // 1
            {0.0d,+a,-b}, // 2
            {0.0d,+a,+b}, // 3
            {-a,-b,0.0d}, // 4
            {-a,+b,0.0d}, // 5
            {+a,-b,0.0d}, // 6
            {+a,+b,0.0d}, // 7
            {-b,0.0d,-a}, // 8
            {+b,0.0d,-a}, // 9
            {-b,0.0d,+a}, // 10
            {+b,0.0d,+a}, // 11
        };
        static int     Icosahedron_edges [Icosahedron_nedges*2] = {
            0,2,  0,8,  0,4, 0,6,  0,9,
            3,1,  3,10, 3,5, 3,7,  3,11,
            2,8,  8,4,  4,6, 6,9,  9,2,
            1,10, 10,5, 5,7, 7,11, 11,1,
            1,4,  10,8, 5,2, 7,9,  11,6,
            1,6,  10,4, 5,8, 7,2,  11,9
        };
        static int     Icosahedron_ngons [Icosahedron_nfaces  ] = { 3,3,3,3,3,  3,3,3,3,3,    3,3,3,3,3,    3,3,3,3,3  };
        static int     Icosahedron_faces [Icosahedron_nfaces*3] = {
            0,8,2,   0,4,8,   0,6,4,  0,9,6,   0,2,9,
            3,10,1,  3,5,10,  3,7,5,  3,11,7,  3,1,11,
            4,1,10,  8,10,5,  2,5,7,  9,7,11,  6,11,1,
            1,4,6,   10,8,4,  5,2,8,  7,9,2,   11,6,9
        };

        // ==== Icosahedron


    };

    #endif﻿
