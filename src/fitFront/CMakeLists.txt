PROJECT(fitFront)


INCLUDE(${YROOT}/share/yocto/fltk.cmake)

ADD_EXECUTABLE(fitFront main.cpp)
TARGET_LINK_YOCTO(fitFront mk)
