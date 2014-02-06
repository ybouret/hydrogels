PROJECT(fit_gelatine)

FIND_PROGRAM(RCPACK NAMES "${YROOT}/bin/yocto-rcpack")

MESSAGE(STATUS "rcpack=${RCPACK}" )

ADD_EXECUTABLE(fit_gelatine main.cpp)
TARGET_LINK_YOCTO(fit_gelatine mk)

GET_TARGET_PROPERTY(the_exe fit_gelatine LOCATION )
SET(RESOURCES ${CMAKE_CURRENT_SOURCE_DIR}/rc)

ADD_CUSTOM_COMMAND(
	TARGET      fit_gelatine
	POST_BUILD
	COMMAND ${RCPACK} ${the_exe} ${RESOURCES} 
    COMMENT "Packing resources"
    )

INSTALL( TARGETS fit_gelatine RUNTIME DESTINATION bin )

