
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
)

#add_library( MeshBuilder2 OBJECT MeshBuilder2.cpp )

file(GLOB dataStructures_ *.h )
add_custom_target(dataStructures SOURCES ${dataStructures_} )
