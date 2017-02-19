#add_library( Convex2d OBJECT Convex2d.cpp )

file(GLOB math_HEADERS_ *.h )
add_custom_target(math_headers SOURCES ${math_HEADERS_} )
