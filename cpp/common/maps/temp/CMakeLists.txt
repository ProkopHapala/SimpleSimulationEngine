
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
)

add_library( GridMap2D        OBJECT GridMap2D.cpp   )
