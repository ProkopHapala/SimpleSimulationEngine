
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

add_library( Shooter SHARED Shooter.cpp
${COMMON_SRCS}/engine/Ship2D.cpp
${COMMON_SRCS}/engine/Warrior3D.cpp
${COMMON_SRCS}/engine/Warrior25D.cpp
${COMMON_SRCS}/engine/Turret.cpp
${COMMON_SRCS}/engine/Gun3D.cpp
${COMMON_SRCS}/engine/Projectile3D.cpp
${COMMON_SRCS}/engine/Object3D.cpp
${COMMON_SRCS}/engine/Terrain25D.cpp
${COMMON_SRCS}/engine/Collisions.cpp
${COMMON_SRCS}/dynamics/Body.cpp
${COMMON_SRCS}/dynamics/Body2D.cpp
${COMMON_SRCS}/dynamics/AeroSurf.cpp
${COMMON_SRCS}/dynamics/AeroSurf2D.cpp
${COMMON_SRCS}/dynamics/AeroCraft.cpp
${COMMON_SRCS}/dynamics/DynamicOpt.cpp
${COMMON_SRCS}/dynamics/SoftBody.cpp
# would generate this error : http://stackoverflow.com/questions/19768267/relocation-r-x86-64-32s-against-linking-error
#$<TARGET_OBJECTS:Body>
#$<TARGET_OBJECTS:AeroSurf>
#$<TARGET_OBJECTS:Collisions>
#$<TARGET_OBJECTS:DynamicOpt>
#$<TARGET_OBJECTS:SoftBody>
)

#file(GLOB engine_HEADERS_ *.h )
#add_custom_target(engine_headers SOURCES ${engine_HEADERS_} )


