
#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/molecular )
add_library(Molecular   SHARED Molecular.cpp  $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(RigidMol    SHARED RigidMol.cpp   $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(ReactiveFF  SHARED ReactiveFF.cpp )

add_library(CLCFGO_lib  SHARED CLCFGO_lib.cpp )
add_library(eFF_lib     SHARED eFF_lib.cpp )


