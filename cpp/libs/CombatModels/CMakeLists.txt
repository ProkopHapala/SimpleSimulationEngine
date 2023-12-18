
#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/math         )
include_directories( ${COMMON_SRCS}/CombatModels )
include_directories( ${COMMON_SRCS}/Orbital      )

#include_directories( ${MY_SRC_DIR}/apps/OrbitalWar   )
include_directories( ${MY_SRC_DIR}/apps/LandTactics   )
set( LTdir ${MY_SRC_DIR}/apps/LandTactics )
#add_library(Molecular   SHARED Molecular.cpp  $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(CombatModels   SHARED CombatModels.cpp   )
add_library(LandCombatLib  SHARED LandCombatLib.cpp ${LTdir}/LTUnit.cpp ${LTdir}/LTUnitType.cpp ${LTdir}/LTcommon.h $<TARGET_OBJECTS:Body2D> )
add_library(SpaceCombatLib SHARED SpaceCombatLib.cpp )


