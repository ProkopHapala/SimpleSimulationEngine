
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
)

add_library( SDL2OGL      OBJECT AppSDL2OGL.cpp ScreenSDL2OGL.cpp ScreenSDL2OGL_3D.cpp drawMath.cpp drawMath2D.cpp )
