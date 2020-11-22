
#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/CombatModels )
#add_library(Molecular   SHARED Molecular.cpp  $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(CombatModels  SHARED CombatModels.cpp )


