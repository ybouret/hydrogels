PROJECT(rdiff1d)

ADD_EXECUTABLE(rdiff1d main.cpp 
	config.lua boum.lua simple.lua gelatine.lua setup.lua
	collection.hpp collection.cpp
	chemsys.hpp    chemsys.cpp
	initializer.hpp initializer.cpp
	
	params.hpp params.cpp
	cell.hpp   cell.cpp cell-front.cpp
)

TARGET_LINK_YOCTO(rdiff1d chem fame lua++ lua mk)

INSTALL( TARGETS rdiff1d RUNTIME DESTINATION bin )

