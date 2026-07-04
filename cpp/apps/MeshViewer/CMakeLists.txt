
include_directories( ${INCLUDE_DIRS} )

add_executable       ( meshViewer meshViewer_main.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( meshViewer ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

file(GLOB meshViewer_HEADERS_ *.h )
add_custom_target( meshViewer_HEADERS SOURCES ${meshViewer_HEADERS_} )

add_executable       ( test_svg_export test_svg_export.cpp )
target_link_libraries( test_svg_export )

add_custom_command(
    TARGET ${PROJECT_NAME} meshViewer
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)
