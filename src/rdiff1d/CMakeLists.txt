PROJECT(rdiff1d)

ADD_EXECUTABLE(rdiff1d main.cpp config.lua
	collection.hpp collection.cpp
)

TARGET_LINK_YOCTO(rdiff1d chemical spade thr lua++ lua mk)

INSTALL( TARGETS rdiff1d RUNTIME DESTINATION bin )

