
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
)

add_library( TiledView      OBJECT TiledView.cpp    )

add_library( SDL2OGL      OBJECT AppSDL2OGL_3D.cpp AppSDL2OGL.cpp PlotScreen2D.cpp ScreenSDL2OGL_3D.cpp ScreenSDL2OGL.cpp  Plot2D.cpp Draw.cpp Draw2D.cpp Draw3D.cpp GUI.cpp 
#$<TARGET_OBJECTS:Mesh> 
${COMMON_SRCS}/math/Mesh.cpp
)

#add_library( GLView      OBJECT GLView.cpp Plot2D.cpp Draw.cpp Draw2D.cpp Draw3D.cpp GUI.cpp ${COMMON_SRCS}/math/Mesh.cpp )

#add_library( GUI      OBJECT GUI.cpp Draw2D.cpp Draw3D.cpp )

file(GLOB SDL2OGL_HEADERS_ *.h *.cpp )
add_custom_target(SDL2OGL_headers SOURCES ${SDL2OGL_HEADERS_} )
