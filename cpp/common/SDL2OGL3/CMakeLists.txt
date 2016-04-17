
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/SDL2OGL3
)

add_library( SDL2OGL3   OBJECT ScreenSDL2OGL3.cpp SceneNode.cpp GLObject.cpp Shader.cpp )

file(GLOB SDL2OGL3_HEADERS_ *.h )
add_custom_target(SDL2OGL3_headers SOURCES ${SDL2OGL3_HEADERS_} )
