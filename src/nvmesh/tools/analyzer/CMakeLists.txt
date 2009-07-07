
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(nvmeshanalyzer MeshAnalyzer.cpp MeshAnalysis.h MeshAnalysis.cpp)
TARGET_LINK_LIBRARIES(nvmeshanalyzer ${LIBS} nvcore nvmath nvimage nvmesh)

INSTALL(TARGETS nvmeshanalyzer DESTINATION bin)
