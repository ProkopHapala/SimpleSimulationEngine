
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/dynamics/Molecular
)

add_library( MoleculeWorld2D OBJECT MoleculeWorld2D.cpp )
add_library( MolecularEngine OBJECT MolecularWorld.cpp MoleculeType.cpp AtomTypes.cpp  )

file(GLOB Molecular_HEADERS_ *.h )
add_custom_target(Molecular_headers SOURCES ${Molecular_HEADERS_} )
