
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
)

add_library( TiledView      OBJECT TiledView.cpp    )
add_library( TerrainCubic   OBJECT TerrainCubic.cpp )
add_library( TerrainRBF     OBJECT TerrainRBF.cpp   )


file(GLOB maps_HEADERS_ *.h )
add_custom_target(maps_headers SOURCES ${maps_HEADERS_} )
