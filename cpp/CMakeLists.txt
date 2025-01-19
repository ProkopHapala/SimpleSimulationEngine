cmake_minimum_required ( VERSION 2.8 )
project ( SimpleSimulationEngine )

option(RELEASE     "build optimized release (-O3 -Ofast) rather than debug? (-g -Og)"               OFF)
option(WITH_SDL    "build with GUI and 3D graphics (SDL+OpenGL)"                                    OFF)
option(WITH_OMP    "use OpenMP paralelization"                                                      OFF)
option(WITH_OPENCL "build with OpenCL GPU acceleration? (can accelare >100x some forcefields and grids  )" OFF)
option(WITH_LUA    "build with LUA scripting                   "                                    OFF)
option(WITH_ASAN   "use runtime memory sanitizer (asan) to trace segfaults and leaks?"              OFF)
option(WITH_NET    "build with GUI and 3D graphics (SDL+OpenGL)"                                    OFF)
option(WITH_MUSIC  "build with GUI and 3D graphics (SDL+OpenGL)"                                    OFF)
message("OPTIONS: -DRELEASE=${RELEASE} -DWITH_SDL=${WITH_SDL} -DWITH_OMP=${WITH_OMP} -DWITH_OPENCL=${WITH_OPENCL} -DWITH_LUA=${WITH_LUA} -DWITH_NET=${WITH_NET} -DWITH_ASAN=${WITH_ASAN}" )

# ==============================
# ==== FLAGS
# ==============================

if(RELEASE)
    set(CMAKE_BUILD_TYPE Release)
else()
    set(CMAKE_BUILD_TYPE Debug)
endif()

set(default_build_type Release)
if( UNIX )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=gnu++0x" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++20" )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -fPIC -fno-strict-aliasing" )
    #SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++2a" )
endif()

message ( "CMAKE_BUILD_TYPE "   ${CMAKE_BUILD_TYPE} )
if (  CMAKE_BUILD_TYPE MATCHES Debug )
    if ( WITH_ASAN )
        SET( AXULIARY_COMPILE_FLAGS "-g -Og -fsanitize=address -fno-omit-frame-pointer")
        #SET( AXULIARY_COMPILE_FLAGS "-g -Og -fsanitize=address -fno-omit-frame-pointer -static-libasan") # static-asan does not work see https://github.com/google/sanitizers/issues/796
    else()
        #SET( AXULIARY_COMPILE_FLAGS "-Wall -g -Og")
        SET( AXULIARY_COMPILE_FLAGS "-g -Og")
    endif()
    #set(CMAKE_CXX_CPPCHECK "cppcheck")
    message("DEBUG BUILD FFLAGS : " )
else()
    #SET( AXULIARY_COMPILE_FLAGS "-w -O3 -fPIC ")
    #SET( AXULIARY_COMPILE_FLAGS "-w -O2 -march=native")
    #SET( AXULIARY_COMPILE_FLAGS "-w -O2 -march=native")
    SET( AXULIARY_COMPILE_FLAGS "-w -Ofast -march=native")
    #SET( AXULIARY_COMPILE_FLAGS "-w -Ofast -march=native -fopt-info-vec-all")
endif()
if ( WITH_OMP )
    find_package(OpenMP)
    #if(OpenMP_CXX_FOUND)
    #    target_link_libraries(MyTarget PUBLIC OpenMP::OpenMP_CXX)
    #endif()
    SET( AXULIARY_COMPILE_FLAGS "${AXULIARY_COMPILE_FLAGS} -fopenmp" )
endif()

SET( IGNORE_WARRNING_FLAGS "-Wno-int-to-pointer-cast -Wno-switch -Wno-reorder -Wno-deprecated -Wno-maybe-uninitialized -Wno-char-subscripts -Wno-write-strings -Wno-format -Wno-parentheses -Wno-unused-but-set-variable  -Wno-narrowing -Wno-unused-result  -Wno-sign-compare -Wno-strict-aliasing -Wno-unused-variable -Wno-unused-value -Wno-comment -Wno-misleading-indentation " )
SET( WARRNING_TO_ERROR     "-Werror=return-type -Werror=init-self -Werror=uninitialized -Werror=implicit-fallthrough -Werror=tautological-compare -Werror=delete-non-virtual-dtor")
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${AXULIARY_COMPILE_FLAGS} ${IGNORE_WARRNING_FLAGS} ${WARRNING_TO_ERROR}" )
message("BUILD FFLAGS : " )
message ( ${CMAKE_CXX_FLAGS} )

# ==============================
# ==== DIRS
# ==============================

SET( MY_SRC_DIR    "${CMAKE_SOURCE_DIR}" )
SET( MY_BIN_SHARE  "${MY_SRC_DIR}/Build"       )

SET( COMMON_SRCS      "${MY_SRC_DIR}/common"     )
SET( COMMON_SDL_SRCS  "${MY_SRC_DIR}/common_SDL" )
SET( COMMON_RESOURCES "${MY_SRC_DIR}/common_resources"       )
SET( COMMON_SHADERS   "${MY_SRC_DIR}/common_resources/shaders" )

include_directories(
        ${COMMON_SRCS}
        ${COMMON_SRCS}/algorithms
        ${COMMON_SRCS}/dataStructures
        ${COMMON_SRCS}/math
        ${COMMON_SRCS}/geometry
        ${COMMON_SRCS}/dynamics
        ${COMMON_SRCS}/maps
        #${COMMON_SRCS}/optimization
        ${COMMON_SRCS}/utils
        ${COMMON_SRCS}/CombatModels
        ${COMMON_SRCS}/engine
)

set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake_utils )

#file (GLOB_RECURSE _ALL_SRCS_ *.h *.cpp *.c )
file (GLOB_RECURSE _COMMON_SRCS_     ${COMMON_SRCS}/*.h         ${COMMON_SRCS}/*.cpp       ${COMMON_SRCS}/*.c     )
file (GLOB_RECURSE _COMMON_SDL_SRCS_ ${COMMON_SDL_SRCS}/*.h     ${COMMON_SDL_SRCS}/*.cpp   ${COMMON_SDL_SRCS}/*.c )
file (GLOB_RECURSE _COMMON_SHADERS_  ${COMMON_SHADERS}/*.glslv  ${COMMON_SHADERS}/*.glslf                         )

# ==============================
# ==== BUILD TARGETS
#===============================

add_subdirectory( ${MY_SRC_DIR}/common   ${MY_BIN_SHARE}/common    )
add_subdirectory( ${MY_SRC_DIR}/libs     )
add_subdirectory( ${MY_SRC_DIR}/tests    )   # currently empty

# ====== OpenCL

if ( WITH_OPENCL )
    include_directories(${COMMON_SRCS}/OCL)
    add_subdirectory( ${MY_SRC_DIR}/libs_OCL )
endif()

# ====== SDL@OpenGL

if ( WITH_SDL )
    # https://cmake.org/cmake/help/latest/policy/CMP0072.html
    #set(OpenGL_GL_PREFERENCE LEGACY)
    set(OpenGL_GL_PREFERENCE GLVND)
    include_directories(
        ${COMMON_SDL_SRCS}
        ${COMMON_SDL_SRCS}/network
        ${COMMON_SDL_SRCS}/SDL2
        ${COMMON_SDL_SRCS}/SDL2OGL
        ${COMMON_SDL_SRCS}/SDL2OGL3
    )
    file (GLOB_RECURSE _COMMON_SDL_SRCS_ ${COMMON_SDL_SRCS}/*.h ${COMMON_SDL_SRCS}/*.cpp ${COMMON_SDL_SRCS}/*.c )
    add_subdirectory( ${MY_SRC_DIR}/common_SDL   ${MY_BIN_SHARE}/common_SDL )
    add_subdirectory( ${MY_SRC_DIR}/libs_SDL     )
    add_subdirectory( ${MY_SRC_DIR}/sketches_SDL )
    add_subdirectory( ${MY_SRC_DIR}/apps         )
    if ( WITH_OPENCL AND WITH_SDL )
        add_subdirectory( ${MY_SRC_DIR}/apps_OCL     )
        add_subdirectory( ${MY_SRC_DIR}/sketches_OCL )
    endif()

endif()

# ====== OpenCL + SDL@OpenGL


