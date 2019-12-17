cmake_minimum_required ( VERSION 2.8 )
project ( SimpleSimulationEngine )

# ==== FLAGS

set(default_build_type Release)
if( UNIX )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++2a" )
endif()


# https://stackoverflow.com/questions/7724569/debug-vs-release-in-cmake
#message ( "CMAKE_CXX_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG}  " )
#message ( "CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE}" )
message ( "CMAKE_BUILD_TYPE   ${CMAKE_BUILD_TYPE}  " )
#if (  (CMAKE_BUILD_TYPE) AND (CMAKE_BUILD_TYPE MATCHES Debug))
#    message("debug mode FFLAGS : " )
#    SET( AXULIARY_COMPILE_FLAGS "-Wall -g -Og -fPIC ")
#else ()
#    message("release mode FFLAGS : " )
#    SET( AXULIARY_COMPILE_FLAGS "-w -O2 -fPIC -march=native ")
#endif ()
# -O3 and -Ofast does not work for SDL2_net
#SET( AXULIARY_COMPILE_FLAGS "-w -O3") 
#SET( AXULIARY_COMPILE_FLAGS "-w -O2 -fPIC -march=native ")
#SET( AXULIARY_COMPILE_FLAGS "-O2 -fPIC -march=native ")
SET( IGNORE_WARRNING_FLAGS "-Wno-maybe-uninitialized -Wno-char-subscripts -Wno-write-strings -Wno-format -Wno-parentheses -Wno-unused-but-set-variable  -Wno-narrowing -Wno-unused-result  -Wno-sign-compare -Wno-strict-aliasing -Wno-unused-variable -Wno-unused-value -Wno-comment -Wno-misleading-indentation " )
SET( AXULIARY_COMPILE_FLAGS "-O2 -fPIC -march=native -fno-strict-aliasing " )

#SET( AXULIARY_COMPILE_FLAGS "-w -g -Og -fPIC ")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS} ${IGNORE_WARRNING_FLAGS}" )
message ( "CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} )

# ==== DIRS

SET( MY_SRC_DIR    "${CMAKE_SOURCE_DIR}" )
SET( MY_BIN_SHARE  "${MY_SRC_DIR}/Build"       )


SET( COMMON_SRCS       "${MY_SRC_DIR}/common"     )
SET( COMMON_SDL_SRCS  "${MY_SRC_DIR}/common_SDL" )
SET( COMMON_RESOURCES "${MY_SRC_DIR}/common_resources"       )
SET( COMMON_SHADERS   "${MY_SRC_DIR}/common_resources/shaders" )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/algorithms
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/maps
        #${COMMON_SRCS}/optimization
        ${COMMON_SRCS}/utils
        ${COMMON_SRCS}/dataStructures
        ${COMMON_SRCS}/engine
        ${COMMON_SDL_SRCS}
        ${COMMON_SDL_SRCS}/network
        ${COMMON_SDL_SRCS}/SDL2
        ${COMMON_SDL_SRCS}/SDL2OGL
        ${COMMON_SDL_SRCS}/SDL2OGL3
        #${MY_SRC_DIR}/engines
)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_utils )

#file (GLOB_RECURSE _ALL_SRCS_ *.h *.cpp *.c )
file (GLOB_RECURSE _COMMON_SRCS_     ${COMMON_SRCS}/*.h     ${COMMON_SRCS}/*.cpp     ${COMMON_SRCS}/*.c     )
file (GLOB_RECURSE _COMMON_SDL_SRCS_ ${COMMON_SDL_SRCS}/*.h ${COMMON_SDL_SRCS}/*.cpp ${COMMON_SDL_SRCS}/*.c )
file (GLOB_RECURSE _COMMON_SHADERS_  ${COMMON_SHADERS}/*.glslv  ${COMMON_SHADERS}/*.glslf    )

# ==== BUILD TARGETS

add_subdirectory( ${MY_SRC_DIR}/common   ${MY_BIN_SHARE}/common    )
add_subdirectory( ${MY_SRC_DIR}/libs     )
add_subdirectory( ${MY_SRC_DIR}/tests    )   # currently empty

# ---- these need SDL
add_subdirectory( ${MY_SRC_DIR}/common_SDL   ${MY_BIN_SHARE}/common_SDL )
add_subdirectory( ${MY_SRC_DIR}/libs_SDL     )
add_subdirectory( ${MY_SRC_DIR}/sketches_SDL )
#add_subdirectory( ${MY_SRC_DIR}/engines      )
add_subdirectory( ${MY_SRC_DIR}/apps         )
