ADD_EXECUTABLE(diff2d main.cpp skin.lua
	types.hpp	types.cpp
	library.hpp	library.cpp
	chemsys.hpp	chemsys.cpp
	param.hpp	param.cpp
	workspace.hpp	workspace.cpp
	initializer.hpp initializer.cpp
)
TARGET_LINK_YOCTO(diff2d aqueous spade thr lua lua++ mk fs)

INSTALL( TARGETS diff2d RUNTIME DESTINATION bin )

