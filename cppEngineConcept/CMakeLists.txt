cmake_minimum_required ( VERSION 2.8 )
project ( NewSimpleSimulationEngine )

if( UNIX )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
endif()

#SET( AXULIARY_COMPILE_FLAGS "-w -O3") 
# -O3 and -Ofast does not work for SDL2_net

SET( AXULIARY_COMPILE_FLAGS "-w -O2 -fPIC")

SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS}" )

SET( COMMON_SRCS "${CMAKE_SOURCE_DIR}" )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/physics
        ${COMMON_SRCS}/graphics
		${COMMON_SRCS}/engine
)

#set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_utils )

#add_subdirectory( engine )
add_subdirectory( tests  )
