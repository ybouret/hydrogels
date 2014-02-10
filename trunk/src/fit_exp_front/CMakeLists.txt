PROJECT(fitExpFront)


INCLUDE(${YROOT}/share/yocto/fltk.cmake)

ADD_EXECUTABLE(fitExpFront main.cpp)
TARGET_LINK_YOCTO(fitExpFront mk)

INSTALL( TARGETS fitExpFront RUNTIME DESTINATION bin )

