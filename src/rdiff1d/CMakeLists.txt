PROJECT(rdiff1d)

ADD_EXECUTABLE(rdiff1d main.cpp config.lua geometry.lua
	collection.hpp collection.cpp
	chemsys.hpp    chemsys.cpp
	initializer.hpp initializer.cpp
	
	params.hpp params.cpp
	cell.hpp   cell.cpp
)

TARGET_LINK_YOCTO(rdiff1d chemical spade thr lua++ lua mk)

INSTALL( TARGETS rdiff1d RUNTIME DESTINATION bin )

