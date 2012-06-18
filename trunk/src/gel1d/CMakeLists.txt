PROJECT(gel1d)

INCLUDE(${YROOT}/share/yocto/fltk.cmake)

IF(FLTK_FOUND)
	ADD_DEFINITIONS( -DHAS_FLTK=1 )
	FLUID_UIC(ui)
ENDIF()

ADD_EXECUTABLE(gel1d main.cpp config.lua front.lua ${ui_SOURCES})
TARGET_LINK_YOCTO(gel1d aqueous lua lua++ swamp mk)

IF(FLTK_FOUND)
	TARGET_LINK_LIBRARIES(gel1d y-fltk)
	TARGET_LINK_FLTK(gel1d)
ENDIF()
