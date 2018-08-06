
include_directories( ${COMMON_SRCS}/dynamics/Molecular )

#add_library(Molecular SHARED Molecular.cpp $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )
add_library(Flight SHARED libFlight.cpp $<TARGET_OBJECTS:Body> $<TARGET_OBJECTS:AeroSurf> $<TARGET_OBJECTS:AeroCraft> )


