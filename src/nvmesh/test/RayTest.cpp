// This code is in the public domain -- castanyo@yahoo.es

#include <stdio.h>

#include <nvcore/StdStream.h>
#include <nvcore/Ptr.h>
#include <nvcore/Timer.h>

#include <nvmesh/TriMesh.h>
#include <nvmesh/import/MeshImportOBJ.h>
#include <nvmesh/kdtree/MeshKDTree.h>
#include <nvmesh/kdtree/KDTree.h>

using namespace nv;


struct MyMessageHandler : public MessageHandler
{
	void log(const char * str, va_list arg)
	{
		va_list tmp;
		va_copy(tmp, arg);
		vprintf(str, arg);

#if _DEBUG && NV_OS_WIN32
//		static StringBuilder buffer(1024);
//		buffer.format(str, arg);
//		OutputDebugStringA(buffer.str());
#endif

		va_end(tmp);
	}
};


int main(int argc, const char * argv[])
{
	MyMessageHandler messageHandler;
	debug::setMessageHandler(&messageHandler);

	if (argc != 2) {
		printf("Usage: nvraytest filename\n");
		return 0;
	}

	const char * fileName = argv[1];

	Timer timer;

	timer.start();

	AutoPtr<TriMesh> mesh;

	{
		StdInputStream stream(fileName);
		if (stream.isError()) {
			printf("Error, cannot open '%s'\n", fileName);
			return 0;
		}
		
		printf("Importing '%s'\n", fileName);
		MeshImportOBJ importer;
		importer.import(&stream);
		
		printf("\nBuilding triangle mesh:\n");
		mesh = importer.builder().buildTriMesh();
	}

	printf("Triangle count = %u\n", mesh->faceCount());
	printf("Vertex count = %u\n", mesh->vertexCount());

	printf("%d sec\n", timer.elapsed() / 1000);

	timer.start();

	AutoPtr<KDTree> kdTree( buildKDTree(mesh.ptr()) );

	printf("%d sec\n", timer.elapsed() / 1000);

	nvDebug("DONE!\n");
	
	return 0;
}
