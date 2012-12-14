PROJECT(benchmark)

ADD_EXECUTABLE(benchmark main.cpp) 

TARGET_LINK_YOCTO(benchmark spade lua lua++ mk fs)

INSTALL( TARGETS benchmark RUNTIME DESTINATION bin )

