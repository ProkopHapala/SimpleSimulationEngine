
set( NAME "algorithms" )

add_executable( ${NAME} 
	main.cpp 
	${COMMON_SRCS}/testUtils.h
	arrayAlgs.h
)

target_link_libraries( ${NAME} )

