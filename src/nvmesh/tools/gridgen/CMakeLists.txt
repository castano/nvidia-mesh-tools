
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(nvgridgen GridGen.cpp)
TARGET_LINK_LIBRARIES(nvgridgen ${LIBS} nvcore nvmath nvimage nvmesh)

INSTALL(TARGETS nvgridgen DESTINATION bin)
