#ifndef  Solids_h
#define  Solids_h

#include "Vec2.h"
#include "Vec3.h"
#include "CMesh.h"

namespace Solids{

    // ==== Tetrahedron

    const static int Tetrahedron_nverts = 4;
    const static int Tetrahedron_nedges = 6;
    const static int Tetrahedron_ntris  = 4;
    const static int Tetrahedron_nfaces = 4;
    static Vec3d     Tetrahedron_verts   [Tetrahedron_nverts  ] = { {-1.,-1.,-1.}, {+1.,+1.,-1.}, {-1.,+1.,+1.}, {+1.,-1.,+1.} };
    static Vec2i     Tetrahedron_edges   [Tetrahedron_nedges  ] = { 0,1,   0,2,   0,3,    1,2,  1,3,  2,3 };
    static Vec3i     Tetrahedron_tris    [Tetrahedron_ntris   ] = { 0,2,1, 0,1,3, 0,3,2,  1,2,3 };
    static int       Tetrahedron_ngons   [Tetrahedron_nfaces  ] = { 3,     3,     3,      3     };
    static int       Tetrahedron_faces   [Tetrahedron_nfaces*3] = { 0,2,1, 0,1,3, 0,3,2,  1,2,3 }; // order is important for normals

    //const static CMesh Tetrahedron = (CMesh){Tetrahedron_nverts,Tetrahedron_nedges,Tetrahedron_ntris, Tetrahedron_nfaces, NULL,NULL,NULL};
    const static CMesh Tetrahedron = (CMesh){Tetrahedron_nverts,Tetrahedron_nedges,Tetrahedron_ntris,Tetrahedron_nfaces, Tetrahedron_verts, Tetrahedron_edges, Tetrahedron_tris, Tetrahedron_ngons, Tetrahedron_faces};

    // ==== Octahedron

    const static int Octahedron_nverts = 6;
    const static int Octahedron_nedges = 12;
    const static int Octahedron_nfaces = 8;
    const static int Octahedron_ntris  = 8;
    static Vec3d     Octahedron_verts   [Octahedron_nverts  ] = { {-1.,0.,0.}, {+1.,0.,0.}, {0.,-1.,0.}, {0.,+1.,0.}, {0.,0.,-1.}, {0.,0.,+1.} };
    static Vec2i     Octahedron_edges   [Octahedron_nedges  ] = { 0,2, 0,3, 0,4, 0,5,   1,2, 1,3, 1,4, 1,5,   2,4, 2,5, 3,4, 3,5  };
    static Vec3i     Octahedron_tris    [Octahedron_ntris   ] = { 0,4,2, 0,2,5, 0,3,4, 0,5,3,   1,2,4, 1,5,2, 1,4,3, 1,3,5 };
    static int       Octahedron_ngons   [Octahedron_nfaces  ] = { 3,     3,     3,     3,       3,     3,     3,     3     };
    static int       Octahedron_faces   [Octahedron_nfaces*3] = { 0,4,2, 0,2,5, 0,3,4, 0,5,3,   1,2,4, 1,5,2, 1,4,3, 1,3,5 };

    const static CMesh Octahedron = (CMesh){Octahedron_nverts,Octahedron_nedges,Octahedron_ntris,Octahedron_nfaces, Octahedron_verts, Octahedron_edges, Octahedron_tris, Octahedron_ngons, Octahedron_faces};


    // ==== Cube

    const static int  Cube_nverts = 8;
    const static int  Cube_nedges = 12;
    const static int  Cube_ntris  = 12;
    const static int  Cube_nfaces = 6;
    static Vec3d      Cube_verts  [Cube_nverts] = {
        {-1.,-1.,-1.},
        {-1.,-1.,+1.},
        {-1.,+1.,-1.},
        {-1.,+1.,+1.},
        {+1.,-1.,-1.},
        {+1.,-1.,+1.},
        {+1.,+1.,-1.},
        {+1.,+1.,+1.},
    };
    static Vec2i   Cube_edges [Cube_nedges  ] = { 0,1, 0,2, 0,4,  1,3,1,5, 2,3, 2,6,    7,5, 7,6, 7,3,  5,4, 6,4    };
    static Vec3i   Cube_tris  [Cube_ntris   ] = { 0,1,3, 0,3,2,  0,4,5, 0,5,1,  0,2,6, 0,6,4,  7,5,4,  7,4,6,  7,3,1,   7,1,5,  7,6,2,  7,2,3  };
    static int     Cube_ngons [Cube_nfaces  ] = { 4,        4,          4,         4,         4,         4        };
    static int     Cube_faces [Cube_nfaces*4] = { 0,1,3,2,  0,4,5,1,    0,2,6,4,   7,5,4,6,   7,3,1,5,   7,6,2,3  };

        // drawFace( block, 0, block.pos+rot.a* L.a, Mat3d{ rot.b    ,rot.c    ,rot.a    }, {L.y,L.z} );
        // drawFace( block, 1, block.pos+rot.a*-L.a, Mat3d{ rot.b*-1.,rot.c*-1.,rot.a*-1.}, {L.y,L.z} );
        // drawFace( block, 2, block.pos+rot.b* L.b, Mat3d{ rot.c    ,rot.a    ,rot.b    }, {L.z,L.x} );
        // drawFace( block, 3, block.pos+rot.b*-L.b, Mat3d{ rot.c*-1.,rot.a*-1.,rot.b*-1.}, {L.z,L.x} );
        // drawFace( block, 4, block.pos+rot.c* L.c, Mat3d{ rot.a    ,rot.b    ,rot.c    }, {L.x,L.y} );
        // drawFace( block, 5, block.pos+rot.c*-L.c, Mat3d{ rot.a*-1.,rot.b*-1.,rot.c*-1.}, {L.x,L.y} );

    static Vec3d      Cube_normals  [Cube_nfaces] = { 
        { 1., 0., 0.},
        {-1., 0., 0.},
        { 0., 1., 0.},
        { 0.,-1., 0.},
        { 0., 0., 1.},
        { 0., 0.,-1.}
        
    };
    static Vec3d      Cube_ups      [Cube_nfaces] = { 
        //{ 0., 1., 0.},
        //{ 0.,-1., 0.},
        //{ 0., 0., 1.},
        //{ 0., 0.,-1.},
        //{ 1., 0., 0.},
        //{-1., 0., 0.}

        { 0., 1., 0.},
        { 0., 1., 0.},
        { 0., 0., 1.},
        { 0., 0., 1.},
        { 1., 0., 0.},
        { 1., 0., 0.}

        // { 0., 0.,-1.},
        // { 0., 0., 1.},
        // {-1., 0., 0.},
        // { 1., 0., 0.},
        // { 0.,-1., 0.},
        // { 0., 1., 0.}
    };



    const static CMesh Cube = (CMesh){Cube_nverts,Cube_nedges,Cube_ntris,Cube_nfaces, Cube_verts, Cube_edges, Cube_tris, Cube_ngons, Cube_faces};


    // ==== RhombicDodecahedron

    const static int  RhombicDodecahedron_nverts = 14;
    const static int  RhombicDodecahedron_nedges = 24;
    const static int  RhombicDodecahedron_ntris  = 24;
    const static int  RhombicDodecahedron_nfaces = 12;
    static Vec3d      RhombicDodecahedron_verts [RhombicDodecahedron_nverts] = {
        {-1.,-1.,-1.},
        {-1.,-1.,+1.},
        {-1.,+1.,-1.},
        {-1.,+1.,+1.},
        {+1.,-1.,-1.},
        {+1.,-1.,+1.},
        {+1.,+1.,-1.},
        {+1.,+1.,+1.},

        {-2.,0.,0.},
        {+2.,0.,0.},
        {0.,-2.,0.},
        {0.,+2.,0.},
        {0.,0.,-2.},
        {0.,0.,+2.}
    };
    static Vec2i   RhombicDodecahedron_edges [RhombicDodecahedron_nedges  ] = { 0,8, 0,10, 0,12,   1,8, 1,10, 1,13,   2,8, 2,11, 2,12,  3,8, 3,11, 3,13,  4,9, 4,10, 4,12,  5,9, 5,10, 5,13,  6,9, 6,11, 6,12,   7,9, 7,11, 7,13   };
    static Vec3i   RhombicDodecahedron_tris  [RhombicDodecahedron_nedges  ] = { 0,10,1, 0,1,8,   0,8,2, 0,2,12,   0,12,4, 0,4,10,    7,11,3, 7,3,13,   7,9,6, 7,6,11,   7,13,5, 7,5,9,   1,13,3, 1,3,8,  1,10,5, 1,5,13,   4,9,5, 4,5,10,  4,12,6, 4,6,9,  2,8,3, 2,3,11,  2,11,6, 2,6,12 };
    static int     RhombicDodecahedron_ngons [RhombicDodecahedron_nfaces  ] = { 4,4,4,4,  4,4,4,4,    4,4,4,4  };
    static int     RhombicDodecahedron_faces [RhombicDodecahedron_nfaces*4] = { 0,10,1,8,  0,8,2,12,    0,12,4,10,   7,11,3,13,   7,9,6,11,   7,13,5,9,   1,13,3,8,  1,10,5,13,   4,9,5,10,  4,12,6,9,  2,8,3,11,  2,11,6,12    };

    const static CMesh RhombicDodecahedron = (CMesh){RhombicDodecahedron_nverts,RhombicDodecahedron_nedges,RhombicDodecahedron_ntris,RhombicDodecahedron_nfaces, RhombicDodecahedron_verts, RhombicDodecahedron_edges, RhombicDodecahedron_tris, RhombicDodecahedron_ngons,RhombicDodecahedron_faces};


    // ==== Icosahedron

    const static int  Icosahedron_nverts = 12;
    const static int  Icosahedron_nedges = 30;
    const static int  Icosahedron_ntris  = 20;
    const static int  Icosahedron_nfaces = 20;
    const double a = 1.;
    const double b = GOLDEN_RATIO;
    static Vec3d      Icosahedron_verts [RhombicDodecahedron_nverts] = {
        {0.,-a,-b}, // 0
        {0.,-a,+b}, // 1
        {0.,+a,-b}, // 2
        {0.,+a,+b}, // 3
        {-a,-b,0.}, // 4
        {-a,+b,0.}, // 5
        {+a,-b,0.}, // 6
        {+a,+b,0.}, // 7
        {-b,0.,-a}, // 8
        {+b,0.,-a}, // 9
        {-b,0.,+a}, // 10
        {+b,0.,+a}, // 11
    };
    static Vec2i     Icosahedron_edges [Icosahedron_nedges] = {
        0,2,  0,8,  0,4, 0,6,  0,9,
        3,1,  3,10, 3,5, 3,7,  3,11,
        2,8,  8,4,  4,6, 6,9,  9,2,
        1,10, 10,5, 5,7, 7,11, 11,1,
        1,4,  10,8, 5,2, 7,9,  11,6,
        1,6,  10,4, 5,8, 7,2,  11,9
    };
    static Vec3i     Icosahedron_tris [Icosahedron_ntris] = {
        0,8,2,   0,4,8,   0,6,4,  0,9,6,   0,2,9,
        3,10,1,  3,5,10,  3,7,5,  3,11,7,  3,1,11,
        4,1,10,  8,10,5,  2,5,7,  9,7,11,  6,11,1,
        1,4,6,   10,8,4,  5,2,8,  7,9,2,   11,6,9
    };
    static int     Icosahedron_ngons [Icosahedron_nfaces  ] = { 3,3,3,3,3,  3,3,3,3,3,    3,3,3,3,3,    3,3,3,3,3  };
    static int     Icosahedron_faces [Icosahedron_nfaces*3] = {
        0,8,2,   0,4,8,   0,6,4,  0,9,6,   0,2,9,
        3,10,1,  3,5,10,  3,7,5,  3,11,7,  3,1,11,
        4,1,10,  8,10,5,  2,5,7,  9,7,11,  6,11,1,
        1,4,6,   10,8,4,  5,2,8,  7,9,2,   11,6,9
    };

    const static CMesh Icosahedron = (CMesh){Icosahedron_nverts,Icosahedron_nedges,Icosahedron_nfaces,Icosahedron_ntris, Icosahedron_verts, Icosahedron_edges, Icosahedron_tris, Icosahedron_ngons, Icosahedron_faces};


    // ==== Icosahedron


};

#endif
