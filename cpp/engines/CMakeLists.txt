
include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/dynamics
)

add_library( Collisions OBJECT Collisions.cpp )

#add_library( Shooter OBJECT Shooter.cpp Warrior3D.cpp Projectile3D.cpp Object3D.cpp Terrain25D.h 
add_library( Shooter SHARED Ship2D.cpp Shooter.cpp Warrior3D.cpp Warrior25D.cpp Turret.cpp Gun3D.cpp Projectile3D.cpp Object3D.cpp Terrain25D.h Collisions.cpp
${COMMON_SRCS}/dynamics/Body.cpp
${COMMON_SRCS}/dynamics/Body2D.cpp
${COMMON_SRCS}/dynamics/AeroSurf.cpp
${COMMON_SRCS}/dynamics/AeroSurf2D.cpp
#${COMMON_SRCS}/dynamics/Collisions.cpp
${COMMON_SRCS}/dynamics/DynamicOpt.cpp
${COMMON_SRCS}/dynamics/SoftBody.cpp
# would generate this error : http://stackoverflow.com/questions/19768267/relocation-r-x86-64-32s-against-linking-error
#$<TARGET_OBJECTS:Body>
#$<TARGET_OBJECTS:AeroSurf>
#$<TARGET_OBJECTS:Collisions>
#$<TARGET_OBJECTS:DynamicOpt>
#$<TARGET_OBJECTS:SoftBody> 
)

#add_library( Engines OBJECT Shooter )

file(GLOB engines_HEADERS_ *.h )
add_custom_target(engines_HEADERS SOURCES ${engines_HEADERS_} )
