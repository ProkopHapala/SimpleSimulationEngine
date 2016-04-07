cmake_minimum_required ( VERSION 2.8 )
project ( SimpleSimulationEngine )

if( UNIX )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
endif()

#SET( AXULIARY_COMPILE_FLAGS "-w -O3") 
# -O3 and -Ofast does not work for SDL2_net

SET( AXULIARY_COMPILE_FLAGS "-w -O2")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS}" )

SET( COMMON_SRCS "${CMAKE_SOURCE_DIR}/common" )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/algorithms
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/maps
        ${COMMON_SRCS}/optimization
        ${COMMON_SRCS}/SDL2OGL
        ${COMMON_SRCS}/network
        ${COMMON_SRCS}/SDL2OGL3
)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_utils )

add_subdirectory( common )
add_subdirectory( apps )
add_subdirectory( tests )
