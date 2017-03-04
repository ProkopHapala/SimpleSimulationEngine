
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
)

add_library(Shock SHARED Shock.cpp $<TARGET_OBJECTS:ShockWaves> )

