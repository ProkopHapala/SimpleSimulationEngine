
#include_directories( ${COMMON_SRCS}/dynamics/Molecular )
include_directories( ${COMMON_SRCS}/math         )
include_directories( ${COMMON_SRCS}/CombatModels )
include_directories( ${COMMON_SRCS}/Orbital      )

#include_directories( ${MY_SRC_DIR}/apps/OrbitalWar   )

set( LTdir ${MY_SRC_DIR}/apps/LandTactics )
set( LCdir ${MY_SRC_DIR}/apps/LandCraft )
include_directories( ${LTdir} )
include_directories( ${LCdir} )

#add_library(Molecular   SHARED Molecular.cpp  $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(CombatModels   SHARED CombatModels.cpp   )
add_library(SpaceCombatLib SHARED SpaceCombatLib.cpp )
add_library(LandCombatLib  SHARED LandCombatLib.cpp ${LTdir}/LTUnit.cpp ${LTdir}/LTUnitType.cpp ${LTdir}/LTcommon.h $<TARGET_OBJECTS:Body2D> )
add_library(LandCraftLib   SHARED LandCraftLib.cpp  ${LCdir}/LandCraftWorld.cpp  $<TARGET_OBJECTS:TerrainGrid2D> $<TARGET_OBJECTS:Noise>  )
