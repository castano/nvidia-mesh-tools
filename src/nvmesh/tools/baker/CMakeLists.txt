
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${NV_SOURCE_DIR}/extern)

IF(OPENEXR_FOUND)
#	SET(LIBS ${LIBS} ${OPENEXR_LIBRARIES})
	INCLUDE_DIRECTORIES(${OPENEXR_INCLUDE_PATHS})
ENDIF(OPENEXR_FOUND)

SET(SRCS
	Baker.cpp 
	Samplers.h Samplers.cpp 
	GeometryImage.h 
	DetailedMeshPass.h DetailedMeshPass.cpp
	BaseMeshPass.h BaseMeshPass.cpp
	BaseSurface.h
	AccBaseSurface.h AccBaseSurface.cpp
	TriBaseSurface.h TriBaseSurface.cpp
	TiledTask.h TiledTask.cpp
	CmdOptions.h)

ADD_EXECUTABLE(nvbaker ${SRCS})

TARGET_LINK_LIBRARIES(nvbaker ${LIBS} nvcore nvmath nvmesh nvimage)

INSTALL(TARGETS nvbaker DESTINATION bin)

