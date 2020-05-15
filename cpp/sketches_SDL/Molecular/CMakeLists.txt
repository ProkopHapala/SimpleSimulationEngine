
project(sketches_SDL_Molecular)

#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )

add_executable       ( test_sp3space test_sp3space.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_sp3space ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_spRotations test_spRotations.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:Lingebra>  )
target_link_libraries( test_spRotations ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_CLCFSF test_CLCFSF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> $<TARGET_OBJECTS:Lingebra>  )
target_link_libraries( test_CLCFSF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_ConfDynamics test_ConfDynamics.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt> )
target_link_libraries( test_ConfDynamics ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_Multipoles test_Multipoles.cpp $<TARGET_OBJECTS:SDL2OGL> )
target_link_libraries( test_Multipoles ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_SoftMolecularDynamics test_SoftMolecularDynamics.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_SoftMolecularDynamics ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_eFF test_eFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_eFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_MMFFmini test_MMFFmini.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_MMFFmini ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_FARFF test_FARFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_FARFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RARFF test_RARFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RARFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RARFF2 test_RARFF2.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RARFF2 ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RspFF test_RspFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RspFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_FTRFF test_FTRFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_FTRFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_RRFF test_RRFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_RRFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_EOFF test_EOFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_EOFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_ESFF test_ESFF.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:DynamicOpt>  )
target_link_libraries( test_ESFF ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_executable       ( test_BondAdaptedMesh test_BondAdaptedMesh.cpp $<TARGET_OBJECTS:SDL2OGL>  )
target_link_libraries( test_BondAdaptedMesh ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} )

add_custom_command(
    TARGET ${PROJECT_NAME} test_SoftMolecularDynamics test_eFF
    #COMMAND ln -sf ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    #COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_CURRENT_BINARY_DIR}/data
)
