add_library( MeshBuilder2 OBJECT MeshBuilder2.cpp )

file(GLOB utils_HEADERS_ *.h )
add_custom_target(utils_headers SOURCES ${math_HEADERS_} )
