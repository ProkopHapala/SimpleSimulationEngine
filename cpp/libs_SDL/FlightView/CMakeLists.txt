
include_directories(
        ${MY_SRC_DIR}/apps/AeroCombat
        ${MY_SRC_DIR}/libs/libFlight
)

add_library( FlightView SHARED FlightView.cpp $<TARGET_OBJECTS:SDL2OGL> $<TARGET_OBJECTS:Body> $<TARGET_OBJECTS:AeroSurf> $<TARGET_OBJECTS:AeroCraft> )
#add_library( FlightView SHARED FlightView.cpp $<TARGET_OBJECTS:SDL2OGL> ${OpenGL_LIBRARIES} ${GLU_LIBRARY} ${SDL2_LIBRARY} )
#add_library(Flight SHARED libFlight.cpp $<TARGET_OBJECTS:Body> $<TARGET_OBJECTS:AeroSurf> $<TARGET_OBJECTS:AeroCraft> )
