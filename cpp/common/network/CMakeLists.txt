find_package( OpenGL   REQUIRED )
find_package( GLU      REQUIRED )
find_package( SDL2     REQUIRED )

#set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/local/lib)
find_package( SDL2_net REQUIRED )

include_directories( /usr/local/include/SDL2 )
#link_directories   ( /usr/local/lib )

add_executable       ( client client.cpp  )
target_link_libraries( client ${SDL2_NET_LIBRARY} )

add_executable       ( server server.cpp )
target_link_libraries( server ${SDL2_NET_LIBRARY} )

add_executable       ( test_UDPClient test_UDPClient.cpp  )
target_link_libraries( test_UDPClient ${SDL2_LIBRARY} ${SDL2_NET_LIBRARY} )

add_executable       ( test_UDPServer test_UDPServer.cpp )
target_link_libraries( test_UDPServer ${SDL2_LIBRARY} ${SDL2_NET_LIBRARY} )

add_executable       ( test_UDPNode_client test_UDPNode_client.cpp $<TARGET_OBJECTS:SDL2OGL>            )
target_link_libraries( test_UDPNode_client ${OpenGL_LIBRARY} ${GLU_LIBRARY} ${SDL2_LIBRARY} ${SDL2_NET_LIBRARY} )

add_executable       ( test_UDPNode_server test_UDPNode_server.cpp )
target_link_libraries( test_UDPNode_server ${SDL2_LIBRARY} ${SDL2_NET_LIBRARY}    )


