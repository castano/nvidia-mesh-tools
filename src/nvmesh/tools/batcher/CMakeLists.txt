
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(nvmeshbatcher MeshBatcher.cpp)
TARGET_LINK_LIBRARIES(nvmeshbatcher ${LIBS} nvcore nvmath nvimage nvmesh)

INSTALL(TARGETS nvmeshbatcher DESTINATION bin)
