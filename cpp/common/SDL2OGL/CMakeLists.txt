
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
)

add_library( drawMath        OBJECT drawMath.cpp   )
add_library( drawMath2D      OBJECT drawMath2D.cpp )
add_library( ScreenSDL2OGL   OBJECT ScreenSDL2OGL.cpp ScreenSDL2OGL_3D.cpp )
