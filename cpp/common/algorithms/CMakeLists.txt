
set( NAME "algorithms" )

add_executable( ${NAME} 
	main.cpp 
	${COMMON_SRCS}/testUtils.h
	${CMAKE_CURRENT_SOURCE_DIR}/arrayAlgs.h
)

target_link_libraries( ${NAME} )

