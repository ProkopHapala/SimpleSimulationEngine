
add_library(Molecular SHARED Molecular.cpp $<TARGET_OBJECTS:MolecularEngine> $<TARGET_OBJECTS:DynamicOpt> )

