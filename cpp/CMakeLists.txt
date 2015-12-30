cmake_minimum_required ( VERSION 2.8 )
project ( SimpleSimulationEngine )

if( UNIX )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
endif()

SET( AXULIARY_COMPILE_FLAGS "-w -O3")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS}" )

SET( COMMON_SRCS "${CMAKE_SOURCE_DIR}/common" )

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_utils )

add_subdirectory( apps )
