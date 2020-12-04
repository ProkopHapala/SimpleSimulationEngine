
#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/math    )
include_directories( ${COMMON_SRCS}/CombatModels )
include_directories( ${MY_SRC_DIR}/apps/OrbitalWar    )
#add_library(Molecular   SHARED Molecular.cpp  $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(CombatModels   SHARED CombatModels.cpp   )
add_library(SpaceCombatLib SHARED SpaceCombatLib.cpp )


