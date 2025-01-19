add_library( Convex2d OBJECT Convex2d.cpp )
add_library( Voronoi  OBJECT Voronoi.cpp  )
add_library( voronoi  OBJECT voronoi.cpp  )
add_library( Mesh     OBJECT Mesh.cpp  )
add_library( MeshBuilder2 OBJECT MeshBuilder2.cpp )

file(GLOB geometry_HEADERS_ *.h )
add_custom_target(geometry_headers SOURCES ${math_HEADERS_} )
