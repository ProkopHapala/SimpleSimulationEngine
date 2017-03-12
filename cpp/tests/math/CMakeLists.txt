
#ADD_CUSTOM_TARGET( shaders ALL  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/shaders ${CMAKE_CURRENT_BINARY_DIR}/shaders )

#add_executable       ( test_Diagonalization test_Diagonalization.cpp )
#target_link_libraries( test_Diagonalization )

add_library(Test_Lingebra SHARED Test_Lingebra.cpp $<TARGET_OBJECTS:Lingebra> )

add_custom_command(
     TARGET ${PROJECT_NAME} Test_Lingebra
     COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/py ${CMAKE_CURRENT_BINARY_DIR}/py
#    COMMAND ${CMAKE_COMMAND} -E create_symlink ${MY_SRC_DIR}/common_resources ${CMAKE_CURRENT_BINARY_DIR}/common_resources
)



