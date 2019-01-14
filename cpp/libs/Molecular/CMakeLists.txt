
include_directories( ${COMMON_SRCS}/dynamics/Molecular )
add_library(Molecular SHARED Molecular.cpp $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )

include_directories( ${COMMON_SRCS}/dynamics/Molecular )
add_library(RigidMol SHARED RigidMol.cpp $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )

