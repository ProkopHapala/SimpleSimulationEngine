
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
)

add_library( TiledView      OBJECT TiledView.cpp    )
add_library( TerrainCubic   OBJECT TerrainCubic.cpp )
