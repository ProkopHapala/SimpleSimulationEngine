
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/SDL2OGL
)

add_library( Body       OBJECT Body.cpp       )
add_library( Body2D     OBJECT Body2D.cpp     )
add_library( AeroSurf   OBJECT AeroSurf.cpp   )
add_library( AeroSurf2D OBJECT AeroSurf2D.cpp )
add_library( Collisions OBJECT Collisions.cpp )
add_library( SoftBody   OBJECT SoftBody.cpp   )


file(GLOB dynamics_HEADERS_ *.h )
add_custom_target(dynamics_headers SOURCES ${dynamics_HEADERS_} )
