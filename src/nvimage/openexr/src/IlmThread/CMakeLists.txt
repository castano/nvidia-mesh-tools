PROJECT (IlmThread)

SET(ILMTHREAD_SRCS
	IlmThread.cpp
	IlmThread.h
	IlmThreadMutex.cpp
	IlmThreadMutex.h
	IlmThreadPool.cpp
	IlmThreadPool.h
	IlmThreadSemaphore.cpp
	IlmThreadSemaphore.h)

IF(WIN32)
	SET(ILMTHREAD_SRCS ${ILMTHREAD_SRCS}
		IlmThreadMutexWin32.cpp
		IlmThreadSemaphoreWin32.cpp
		IlmThreadWin32.cpp)
ELSE(WIN32)
	SET(ILMTHREAD_SRCS ${ILMTHREAD_SRCS}
		IlmThreadMutexPosix.cpp
		IlmThreadPosix.cpp
		IlmThreadSemaphorePosix.cpp
		IlmThreadSemaphorePosixCompat.cpp)
ENDIF(WIN32)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../Iex)

IF(ILMBASE_SHARED)
	ADD_LIBRARY(IlmThread SHARED ${ILMTHREAD_SRCS})
ELSE(ILMBASE_SHARED)
	ADD_LIBRARY(IlmThread ${ILMTHREAD_SRCS})
ENDIF(ILMBASE_SHARED)

#TARGET_LINK_LIBRARIES(IlmThread ${LIBS})

INSTALL(TARGETS IlmThread
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib/static)

INSTALL(FILES 
	IlmThread.h
	IlmThreadMutex.h
	IlmThreadPool.h
	IlmThreadSemaphore.h
	DESTINATION include/OpenEXR)
