PROJECT(fitGelFront)


INCLUDE(${YROOT}/share/yocto/fltk.cmake)

ADD_EXECUTABLE(fitGelFront main.cpp)
TARGET_LINK_YOCTO(fitGelFront mk)

INSTALL( TARGETS fitGelFront RUNTIME DESTINATION bin )

