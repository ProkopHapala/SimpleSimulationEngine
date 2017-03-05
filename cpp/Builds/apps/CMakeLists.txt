cmake_minimum_required ( VERSION 2.8 )
project ( SimpleSimulationEngine_libs )

# ==== FLAGS

if( UNIX )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
endif()
#SET( AXULIARY_COMPILE_FLAGS "-w -O3") 
# -O3 and -Ofast does not work for SDL2_net
SET( AXULIARY_COMPILE_FLAGS "-w -O2 -fPIC")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS}" )

# ==== FLAGS

SET( MY_SRC_DIR    "${CMAKE_SOURCE_DIR}/../../" )
#SET( MY_BIN_SHARE  "${MY_SRC_DIR}/Builds/shared_bins"       )
SET( MY_BIN_SHARE  "${MY_SRC_DIR}/Build"       )

SET( COMMON_SRCS "${MY_SRC_DIR}common" )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/algorithms
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/maps
        ${COMMON_SRCS}/optimization
        ${COMMON_SRCS}/network
	${COMMON_SRCS}/engines
        ${COMMON_SRCS}/SDL2OGL
        ${COMMON_SRCS}/SDL2OGL3
)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${MY_SRC_DIR}cmake_utils )

add_subdirectory( ${MY_SRC_DIR}/common ${MY_BIN_SHARE}/common )
add_subdirectory( ${MY_SRC_DIR}/apps   apps  )
