
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${CMAKE_SOURCE_DIR}/apps/MolecularEditor
)

add_library(Molecular SHARED Molecular.cpp $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )

