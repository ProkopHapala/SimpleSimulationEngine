
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/cpp )

add_library(KosmoSuite SHARED KosmoSuite.cpp )

add_subdirectory (cpp)
