
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
)

add_library( Body       OBJECT Body.cpp       )
add_library( Body2D     OBJECT Body2D.cpp     )
add_library( AeroSurf   OBJECT AeroSurf.cpp   )
add_library( AeroSurf2D OBJECT AeroSurf2D.cpp )
add_library( AeroCraft  OBJECT AeroCraft.cpp  )

add_library( DynamicOpt OBJECT DynamicOpt.cpp )
add_library( SoftBody   OBJECT SoftBody.cpp   )
add_library( TrussDynamics_d OBJECT TrussDynamics_d.cpp )
add_library( TrussDynamics_f OBJECT TrussDynamics_f.cpp )

add_library( Truss      OBJECT Truss.cpp      )
add_library( ShockWaves OBJECT Shock1D.cpp    )

add_library( CommodityNetwork OBJECT CommodityNetwork.cpp )

#add_library( MoleculeWorld2D OBJECT MoleculeWorld2D.cpp )
#add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )

#add_subdirectory (Molecular)

file(GLOB dynamics_HEADERS_ *.h )
add_custom_target(dynamics_headers SOURCES ${dynamics_HEADERS_} )
