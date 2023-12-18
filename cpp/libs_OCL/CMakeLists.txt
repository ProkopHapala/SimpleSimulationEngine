find_package( OpenCL REQUIRED )
#find_package( clFFT  REQUIRED )

include_directories( ${COMMON_SRCS}/OpenCL    )

add_library(SpaceCraftingOCL_lib SHARED OCL_SpaceCraftingOCL_lib.cpp )
target_link_libraries(SpaceCraftingOCL_lib ${OPENCL_LIBRARY})
#target_link_libraries(SpaceCraftingOCL_lib ${OPENCL_LIBRARY} ${CLFFT_LIBRARIES})