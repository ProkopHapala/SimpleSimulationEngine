#ifndef  Solids_h
#define  Solids_h

class Solid{

	int nvert;
	int nedge;
	int ntri;

	Vec3d * verts;
	int   * edges;
	Vec3i * tris;  // later we may do polygon faces ? 

}

namespace Solids{
	const static int   nTetrahedron_verts = 4;
	const static int   nTetrahedron_edges = 6;
	const static int   nTetrahedron_tris  = 4;
	const static double Tetrahedron_verts [nTetrahedron_verts][3] = { {-1.0d,-1.0d,-1.0d}, {+1.0d,+1.0d,-1.0d}, {-1.0d,+1.0d,+1.0d}, {+1.0d,-1.0d,-1.0d} };
	const static int    Tetrahedron_edges [nTetrahedron_edges][2] = { {0,1},{0,2},{0,3}, {1,2},{1,3},{2,3} };
	const static int    Tetrahedron_tris  [nTetrahedron_tris ][3] = { {0,1,2},{0,1,3},{0,2,3},{1,2,3} };
};

#endif﻿
