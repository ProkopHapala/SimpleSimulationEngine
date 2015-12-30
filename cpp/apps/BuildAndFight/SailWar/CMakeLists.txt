
find_package(OpenGL REQUIRED)
find_package(GLU REQUIRED)
find_package(SDL2 REQUIRED)

include_directories(
        ${COMMON_SRCS}
        ${OpenGL_INCLUDE_DIR}
        ${GLU_INCLUDE_DIR}
        ${SDL2_INCLUDE_DIR}
)

add_executable( SailWar main.cpp )
target_link_libraries( SailWar ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

file( COPY ${CMAKE_CURRENT_SOURCE_DIR}/data DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )