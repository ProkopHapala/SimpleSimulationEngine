find_package( OpenGL REQUIRED )
find_package( GLU    REQUIRED )
find_package( SDL2   REQUIRED )

project(sketches_SDL)

add_subdirectory (math)
add_subdirectory (2D)
add_subdirectory (3D)
add_subdirectory (Molecular)
add_subdirectory (Shooter)
add_subdirectory (OGL3)

if( WITH_MUSIC )
add_subdirectory (music)
endif() # WITH_MUSIC

if( WITH_NET )
add_subdirectory (network)
endif() # WITH_NET