PROJECT(benchmark)

ADD_EXECUTABLE(benchmark main.cpp config.lua) 

TARGET_LINK_YOCTO(benchmark spade mk lua++ lua)

INSTALL( TARGETS benchmark RUNTIME DESTINATION bin )

