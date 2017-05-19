
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
)

file(GLOB dataStructures_ *.h )
add_custom_target(dataStructures SOURCES ${dataStructures_} )
