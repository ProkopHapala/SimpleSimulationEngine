
# Includes from common modules that SpaceWorld and its components depend on
include_directories( ${COMMON_SRCS}/Orbital )
include_directories( ${COMMON_SRCS}/dynamics )
include_directories( ${COMMON_SRCS}/math )

# Define the shared library
# It is built from our new C-wrapper (SpaceTacticsLib.cpp)
# and the core C++ logic it wraps (SpaceWorld.cpp).
# We also link against common object libraries that provide underlying functionality.
add_library(SpaceTacticsLib SHARED SpaceTacticsLib.cpp)

# Pull in required object libraries (e.g., PointBody implementation from Body.cpp)
target_link_libraries(SpaceTacticsLib PRIVATE Body)
