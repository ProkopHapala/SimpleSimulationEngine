
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
)

add_library( TerrainGrid2D  OBJECT  Grid2DAlgs.cpp TerrainHydraulics.cpp )
add_library( Fluid2D        OBJECT Fluid2D.cpp )
add_library( TerrainCubic   OBJECT TerrainCubic.cpp )
add_library( SimplexGrid    OBJECT SimplexGrid.cpp  )
add_library( TerrainSimplex      OBJECT TerrainSimplex.cpp )
add_library( TerrainRBF     OBJECT TerrainRBF.cpp   )
add_library( TrussBuilder   OBJECT TrussBuilder.cpp )

add_library( NBodyWorld2D   OBJECT NBodyWorld2D.cpp   )
#target_compile_options(mylib PRIVATE   -O3 --ffast-math )
#target_compile_options(NBodyWorld2D PRIVATE   -Ofast )

file(GLOB maps_HEADERS_ *.h )
add_custom_target(maps_headers SOURCES ${maps_HEADERS_} )
