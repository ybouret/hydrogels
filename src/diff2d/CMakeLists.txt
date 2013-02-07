ADD_EXECUTABLE(diff2d main.cpp 
	types.hpp	types.cpp
)
TARGET_LINK_YOCTO(diff2d aqueous spade thr lua lua++ mk fs)

INSTALL( TARGETS diff2d RUNTIME DESTINATION bin )

