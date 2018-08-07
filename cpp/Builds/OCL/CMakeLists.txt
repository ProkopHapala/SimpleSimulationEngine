cmake_minimum_required ( VERSION 2.8 )
project ( OCL_SimpleSimulationEngine )

# ==== FLAGS

if( UNIX )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
endif()
#SET( AXULIARY_COMPILE_FLAGS "-w -O3") 
# -O3 and -Ofast does not work for SDL2_net
#SET( AXULIARY_COMPILE_FLAGS "-w -O2 -fPIC")
SET( AXULIARY_COMPILE_FLAGS "-w -O2 -fPIC -march=native ")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS}" )

# ==== DIRS

SET( MY_SRC_DIR    "${CMAKE_SOURCE_DIR}/../.." )
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
        ${COMMON_SRCS}/OCL
        ${COMMON_SRCS}/utils
        ${COMMON_SDL_SRCS}
        ${COMMON_SDL_SRCS}/SDL2
        ${COMMON_SDL_SRCS}/SDL2OGL
        #${COMMON_SDL_SRCS}/SDL2OGL3
)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MY_SRC_DIR}/cmake_utils )

# ==== BUILD TARGETS

#add_subdirectory( ${MY_SRC_DIR}/common   ${MY_BIN_SHARE}/common    )
add_subdirectory( ${MY_SRC_DIR}/common/math           ${MY_BIN_SHARE}/common/math    )
add_subdirectory( ${MY_SRC_DIR}/common/algorithms     ${MY_BIN_SHARE}/common/algorithms    )
add_subdirectory( ${MY_SRC_DIR}/common/dataStructures ${MY_BIN_SHARE}/common/dataStructures    )
add_subdirectory( ${MY_SRC_DIR}/common/maps           ${MY_BIN_SHARE}/common/maps    )
add_subdirectory( ${MY_SRC_DIR}/common/dynamics       ${MY_BIN_SHARE}/common/dynamics    )
add_subdirectory( ${MY_SRC_DIR}/common/OCL            ${MY_BIN_SHARE}/common/OCL    )

#add_subdirectory( ${MY_SRC_DIR}/libs         ${MY_BIN_SHARE}/libs    )
#add_subdirectory( ${MY_SOURCE_DIR}/sketches    )   # currently empty

# ---- these need SDL
#add_subdirectory( ${MY_SRC_DIR}/common_SDL   ${MY_BIN_SHARE}/common_SDL   )
add_subdirectory( ${MY_SRC_DIR}/common_SDL/SDL2      ${MY_BIN_SHARE}/common_SDL/SDL2   )
add_subdirectory( ${MY_SRC_DIR}/common_SDL/SDL2OGL   ${MY_BIN_SHARE}/common_SDL/SDL2OGL   )

#add_subdirectory( ${MY_SRC_DIR}/libs_SDL     ${MY_BIN_SHARE}/libs_SDL    )
add_subdirectory( ${MY_SRC_DIR}/sketches_OCL ${MY_BIN_SHARE}/sketches_OCL )
add_subdirectory( ${MY_SRC_DIR}/apps_OCL     ${MY_BIN_SHARE}/apps_OCL     )
#add_subdirectory( ${MY_SRC_DIR}/engines      ${MY_BIN_SHARE}/engines     )
#add_subdirectory( ${MY_SRC_DIR}/apps         ${MY_BIN_SHARE}/apps        )
