add_library( Lingebra OBJECT Lingebra.cpp )
add_library( Convex2d OBJECT Convex2d.cpp )
add_library( Voronoi  OBJECT Voronoi.cpp  )
add_library( Noise    OBJECT Noise.cpp    )
add_library( Splines  OBJECT SplineManager.cpp  )

file(GLOB math_HEADERS_ *.h )
add_custom_target(math_headers SOURCES ${math_HEADERS_} )
