
project(sketches_OCL)

find_package( OpenGL REQUIRED )
find_package( GLU    REQUIRED )
find_package( SDL2   REQUIRED )
find_package( OpenCL REQUIRED )

include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}/OCL
)

add_subdirectory (OCL)
add_subdirectory (lingebra)
add_subdirectory (NBody)
add_subdirectory (Grids)
