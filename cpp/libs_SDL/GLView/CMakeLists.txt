
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SDL_SRCS}/SDL2OGL
)

#add_library( GLView SHARED GLView.cpp Plot2D.cpp )
#add_library(GLView OBJECT GLView.cpp Plot2D.cpp )

#add_library(GLView SHARED GLView.cpp ${COMMON_SDL_SRCS}/SDL2OGL/Plot2D.cpp )
#add_library(GLView SHARED GLView.cpp $<TARGET_OBJECTS:SDL2OGL> )
#add_library(Vis3D SHARED Vis3D.cpp $<TARGET_OBJECTS:SDL2OGL>)
#add_library(Vis3D SHARED Vis3D.cpp $<TARGET_OBJECTS:SDL2OGL>)

add_library( GLView SHARED GLView.cpp
        ${COMMON_SDL_SRCS}/SDL2OGL/Plot2D.cpp
        ${COMMON_SDL_SRCS}/SDL2OGL/Draw.cpp
        ${COMMON_SDL_SRCS}/SDL2OGL/Draw2D.cpp
        ${COMMON_SDL_SRCS}/SDL2OGL/Draw3D.cpp
        ${COMMON_SDL_SRCS}/SDL2OGL/GUI.cpp
        ${COMMON_SRCS}/math/Mesh.cpp
)

#add_library( GLView SHARED GLView.cpp $<TARGET_OBJECTS:SDL2OGL> ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )
