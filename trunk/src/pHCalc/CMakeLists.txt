PROJECT(pHCalc)

INCLUDE(${YROOT}/share/yocto/fltk.cmake)

IF(FLTK_FOUND)
	FLUID_UIC(calc)
	INCLUDE_DIRECTORIES(.)
	ADD_EXECUTABLE(pHCalc main.cpp main.hpp ${calc_SOURCES})
	TARGET_LINK_YOCTO(pHCalc fltk aqueous mk)
	TARGET_LINK_FLTK(pHCalc)
	
ENDIF()
