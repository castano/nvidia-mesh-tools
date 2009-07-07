
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(nvbzrexport BzrExport.cpp)
TARGET_LINK_LIBRARIES(nvbzrexport ${LIBS} nvcore nvmath nvimage nvmesh)

INSTALL(TARGETS nvbzrexport DESTINATION bin)
