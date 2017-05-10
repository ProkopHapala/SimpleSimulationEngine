find_package( OpenGL REQUIRED )
find_package( GLU    REQUIRED )
find_package( SDL2   REQUIRED )

#add_subdirectory (algorithms)
#add_subdirectory (math)
#add_subdirectory (dynamics)
#add_subdirectory (maps)
add_subdirectory (SDL2)
add_subdirectory (SDL2OGL)
add_subdirectory (SDL2OGL3)
#add_subdirectory (network)     # currently there are only headers
#add_subdirectory (engines)
