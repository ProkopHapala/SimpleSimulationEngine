
project(Apps)

find_package( OpenGL REQUIRED )
find_package( GLU    REQUIRED )
find_package( SDL2   REQUIRED )

add_subdirectory (SailWar)
add_subdirectory (AeroCombat)
add_subdirectory (FormationTactics)
add_subdirectory (MinimalTactics)
add_subdirectory (LandTactics)
add_subdirectory (LandCraft)
add_subdirectory (NonInertial)
add_subdirectory (MultiFight3D)
add_subdirectory (BlockHouseTactics)
add_subdirectory (CastleBuilder)
add_subdirectory (OrbitalWar)
add_subdirectory (CAD)
add_subdirectory (ShapePainter)
add_subdirectory (Tanks)
add_subdirectory (NavalBattle)
add_subdirectory (SwordPlay)
add_subdirectory (MolecularEditor)
#add_subdirectory (MolecularEditor2)
add_subdirectory (DemoCrat)

if( WITH_MUSIC)
add_subdirectory (MusicVizualizer)
endif() # WITH_MUSIC

if( WITH_NET)
add_subdirectory (SailWar_Multi)
endif() # WITH_NET


file (GLOB_RECURSE _COMMON_SRCS_     ${COMMON_SRCS}/*.h     ${COMMON_SRCS}/*.cpp     ${COMMON_SRCS}/*.c     )
file (GLOB_RECURSE _COMMON_SDL_SRCS_ ${COMMON_SDL_SRCS}/*.h ${COMMON_SDL_SRCS}/*.cpp ${COMMON_SDL_SRCS}/*.c )

add_custom_target( ${PROJECT_NAME}_src_dep SOURCES ${_COMMON_SRCS_} ${_COMMON_SDL_SRCS_} )
