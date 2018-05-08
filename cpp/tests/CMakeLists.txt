
# in future there will be termina-based tests which are visualized by external python scripts

project(tests)

add_subdirectory(math)
add_subdirectory(lua)
add_subdirectory(parsing)

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} )