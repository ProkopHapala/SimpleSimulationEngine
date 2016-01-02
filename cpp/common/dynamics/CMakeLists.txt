
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
)

add_library( Body       OBJECT Body.cpp       )
add_library( Body2D     OBJECT Body2D.cpp     )
add_library( AeroSurf2D OBJECT AeroSurf2D.cpp )
add_library( Collisions OBJECT Collisions.cpp )
