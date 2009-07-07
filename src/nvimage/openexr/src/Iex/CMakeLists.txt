PROJECT (Iex)

SET(IEX_SRCS	
	Iex.h
	IexBaseExc.cpp
	IexBaseExc.h
	IexErrnoExc.h
	IexMacros.h
	IexMathExc.h
	IexThrowErrnoExc.cpp
	IexThrowErrnoExc.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF(ILMBASE_SHARED)	
	ADD_LIBRARY(Iex SHARED ${IEX_SRCS})
ELSE(ILMBASE_SHARED)
	ADD_LIBRARY(Iex ${IEX_SRCS})
ENDIF(ILMBASE_SHARED)

#TARGET_LINK_LIBRARIES(Iex ${LIBS})

INSTALL(TARGETS Iex
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib/static)

INSTALL(FILES 
	Iex.h
	IexBaseExc.h
	IexErrnoExc.h
	IexMacros.h
	IexMathExc.h
	IexThrowErrnoExc.h
	DESTINATION include/OpenEXR)