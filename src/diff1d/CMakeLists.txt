PROJECT(diff1d)

SET(FLTK_FOUND FALSE)

IF(FALSE)
INCLUDE(${YROOT}/share/yocto/fltk.cmake)

IF(FLTK_FOUND)
	ADD_DEFINITIONS( -DHAS_FLTK=1 )
	FLUID_UIC(ui)
ENDIF()

ENDIF()

#ADD_EXECUTABLE(diff1d main.cpp config.lua front.lua ${ui_SOURCES})
ADD_EXECUTABLE(diff1d main.cpp config.lua
types.hpp types.cpp
library.cpp library.hpp
chemsys.cpp chemsys.hpp
initializer.cpp initializer.hpp
param.hpp param.cpp
workspace.hpp workspace.cpp
worker.hpp    worker.cpp 
cell.hpp      cell.cpp cell-utils.cpp
)
TARGET_LINK_YOCTO(diff1d aqueous spade thr lua lua++ mk)

IF(FLTK_FOUND)
	TARGET_LINK_LIBRARIES(gel1d y-fltk)
	TARGET_LINK_FLTK(gel1d)
ENDIF()

INSTALL( TARGETS diff1d RUNTIME DESTINATION bin )

