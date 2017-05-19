cmake_minimum_required ( VERSION 2.8 )
project ( SimpleSimulationEngine )

# ==== FLAGS

if( UNIX )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
endif()
#SET( AXULIARY_COMPILE_FLAGS "-w -O3") 
# -O3 and -Ofast does not work for SDL2_net
SET( AXULIARY_COMPILE_FLAGS "-w -O2 -fPIC")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS}" )

# ==== DIRS

SET( MY_SRC_DIR "${CMAKE_SOURCE_DIR}/../.." )
SET( MY_BIN_SHARE  "${MY_SRC_DIR}/Build"       )

SET( COMMON_SRCS     "${MY_SRC_DIR}/common"     )
SET( COMMON_SDL_SRCS "${MY_SRC_DIR}/common_SDL" )
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/algorithms
        ${COMMON_SRCS}/dataStructures
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/maps
        ${COMMON_SRCS}/optimization
	#${COMMON_SDL_SRCS}
        #${COMMON_SDL_SRCS}/network
        #${COMMON_SDL_SRCS}/SDL2
        #${COMMON_SDL_SRCS}/SDL2OGL
        #${COMMON_SDL_SRCS}/SDL2OGL3
        #${MY_SRC_DIR}/engines
)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MY_SRC_DIR}/cmake_utils )

# ==== BUILD TARGETS

add_subdirectory( ${MY_SRC_DIR}/common   ${MY_BIN_SHARE}/common    )
add_subdirectory( ${MY_SRC_DIR}/libs     ${MY_BIN_SHARE}/libs    )
#add_subdirectory( ${MY_SOURCE_DIR}/sketches    )   # currently empty

# ---- these need SDL
#add_subdirectory( ${MY_SRC_DIR}/common_SDL   ${MY_BIN_SHARE}/common_SDL   )
#add_subdirectory( ${MY_SRC_DIR}/libs_SDL     ${MY_BIN_SHARE}/libs_SDL     )
#add_subdirectory( ${MY_SRC_DIR}/sketches_SDL ${MY_BIN_SHARE}/sketches_SDL )
#add_subdirectory( ${MY_SRC_DIR}/engines      ${MY_BIN_SHARE}/engines      )
#add_subdirectory( ${MY_SRC_DIR}/apps         ${MY_BIN_SHARE}/apps         )
