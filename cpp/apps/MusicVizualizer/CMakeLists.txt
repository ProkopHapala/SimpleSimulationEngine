

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/engines
        ${COMMON_SRCS}/SDL2OGL
        ${OpenGL_INCLUDE_DIR}
        ${GLU_INCLUDE_DIR}
        ${SDL2_INCLUDE_DIR}
)

find_library(SDL_MIXER_LIBRARY
  NAMES SDL2_mixer
  HINTS
    ENV SDLMIXERDIR
    ENV SDLDIR
  PATH_SUFFIXES lib
)

#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
find_package( SDL2_mixer REQUIRED )

include_directories( /usr/local/include/SDL2 )
#link_directories   ( /usr/local/lib )



set( NAME "MusicVisualizer" )

find_package(GLEW REQUIRED)
if (GLEW_FOUND)
    include_directories(${GLEW_INCLUDE_DIRS})
    link_libraries(${GLEW_LIBRARIES})
endif()
add_executable( ${NAME}_main
	${NAME}_main.cpp 
	#$<TARGET_OBJECTS:SDL2OGL>
	$<TARGET_OBJECTS:SDL2OGL3>
)
target_link_libraries( ${NAME}_main${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${SDL2_MIXER_LIBRARY} )

add_custom_command(
    TARGET ${NAME}_main
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources   ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
