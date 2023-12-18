find_package( OpenCL REQUIRED )
#find_package( clFFT  REQUIRED )

include_directories( ${COMMON_SRCS}/OpenCL ${COMMON_SRCS}/Orbital  )

add_library(SpaceCombatOCL_lib SHARED SpaceCombatOCL_lib.cpp )
target_link_libraries(SpaceCombatOCL_lib ${OPENCL_LIBRARY})
#target_link_libraries(SpaceCombatOCL_lib ${OPENCL_LIBRARY} ${CLFFT_LIBRARIES})