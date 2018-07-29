
project(sketches_SDL_Shooter)

add_executable       ( test_ShotHit test_ShotHit.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_ShotHit ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

#add_executable       ( test_BroadPhase test_BroadPhase.cpp $<TARGET_OBJECTS:SDL2OGL> )
#target_link_libraries( test_BroadPhase ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable       ( test_SweepAndPrune test_SweepAndPrune.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_SweepAndPrune ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable       ( test_BoxAndSweep test_BoxAndSweep.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_BoxAndSweep ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_executable       ( test_GridHash test_GridHash.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_GridHash ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} Shooter )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_ShotHit test_SweepAndPrune
    #COMMAND ln -sf ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)
