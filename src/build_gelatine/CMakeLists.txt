PROJECT(build_gelatine)

ADD_EXECUTABLE(build_gelatine main.cpp)
TARGET_LINK_YOCTO(build_gelatine mk)

INSTALL( TARGETS build_gelatine RUNTIME DESTINATION bin )

