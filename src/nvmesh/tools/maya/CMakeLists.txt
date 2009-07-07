
IF(MAYA_FOUND)
	SOURCE_GROUP("MEL Scripts" "\\.mel")

	INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR} ${MAYA_INCLUDE_PATH})

	ADD_DEFINITIONS(-DNT_PLUGIN)

	# Shared nvmaya library
	ADD_LIBRARY(nvmaya MayaUtils.h MayaUtils.cpp MayaMeshBuilder.h MayaMeshBuilder.cpp)
	TARGET_LINK_LIBRARIES(nvmaya ${LIBS} ${MAYA_LIBRARIES} nvcore nvmath nvimage nvmesh)


	# MayaAnalyzer plugin
	ADD_LIBRARY(MayaAnalyzer SHARED MayaAnalyzer.cpp)
	TARGET_LINK_LIBRARIES(MayaAnalyzer ${LIBS} ${MAYA_LIBRARIES} nvcore nvmath nvimage nvmesh nvmaya)

	INSTALL(TARGETS MayaAnalyzer DESTINATION bin/maya)

	SET_TARGET_PROPERTIES(MayaAnalyzer
		PROPERTIES
		PREFIX     "Nv"
		SUFFIX     "${MAYA_EXTENSION}")


	# MayaExporter plugin
	ADD_LIBRARY(MayaExporter SHARED MayaExportPlugin.cpp MayaExporter.cpp MayaExporter.h NvMayaExportScript.mel)
	TARGET_LINK_LIBRARIES(MayaExporter ${LIBS} ${MAYA_LIBRARIES} nvcore nvmath nvimage nvmesh nvmaya)

	INSTALL(TARGETS MayaExporter DESTINATION bin/maya)

	SET_TARGET_PROPERTIES(MayaExporter
		PROPERTIES
		PREFIX     "Nv"
		SUFFIX     "${MAYA_EXTENSION}")

	
	# MayaCommandLineExporter tool
	ADD_EXECUTABLE(MayaCommandLineExporter MayaCommandLineExporter.cpp MayaExporter.cpp MayaExporter.h)
	TARGET_LINK_LIBRARIES(MayaCommandLineExporter ${LIBS} ${MAYA_LIBRARIES} nvcore nvmath nvimage nvmesh nvmaya)
	
	INSTALL(TARGETS MayaCommandLineExporter DESTINATION bin)

ENDIF(MAYA_FOUND)

