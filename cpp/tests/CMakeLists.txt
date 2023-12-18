
# in future there will be termina-based tests which are visualized by external python scripts

project(tests)

add_subdirectory(math)
add_subdirectory(parsing)

if( WITH_LUA )
add_subdirectory(lua)
endif() # WITH_LUA

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} )