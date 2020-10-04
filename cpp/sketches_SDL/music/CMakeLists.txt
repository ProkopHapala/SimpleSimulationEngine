
project(sketches_SDL_music)

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

add_executable       ( test_playMP3 test_playMP3.cpp $<TARGET_OBJECTS:SDL2OGL>            )
target_link_libraries( test_playMP3 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${SDL2_MIXER_LIBRARY} )

add_executable       ( test_visualize test_visualize.cpp $<TARGET_OBJECTS:SDL2OGL>            )
target_link_libraries( test_visualize ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${SDL2_MIXER_LIBRARY} )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_playMP3 test_visualize
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
