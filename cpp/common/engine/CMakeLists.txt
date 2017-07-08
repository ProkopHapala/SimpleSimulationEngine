
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
)

add_library( Collisions     OBJECT Collisions.cpp )
add_library( Terrain25D     OBJECT Terrain25D.cpp )
add_library( Object3D       OBJECT Object3D.cpp )
add_library( Object25D      OBJECT Object25D.cpp )
add_library( Projectile3D   OBJECT Projectile3D.cpp )
add_library( Gun3D          OBJECT Gun3D.cpp )
add_library( Turret         OBJECT Turret.cpp )
add_library( Ship2D         OBJECT Ship2D.cpp )
add_library( Warrior3D      OBJECT Warrior3D.cpp )
add_library( Warrior25D     OBJECT Warrior25D.cpp )

file(GLOB engine_HEADERS_ *.h )
add_custom_target(engine_headers SOURCES ${engine_HEADERS_} )
