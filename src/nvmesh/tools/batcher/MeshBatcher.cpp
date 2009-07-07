// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include <stdio.h>
#include <float.h> // FLT_MAX

#include <nvcore/StdStream.h>
#include <nvcore/Ptr.h>
#include <nvcore/RefCounted.h>
#include <nvcore/Timer.h>

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/subdiv/AccMeshBuilder.h>
#include <nvmesh/subdiv/AccMesh.h>
#include <nvmesh/subdiv/AccPatch.h>
#include <nvmesh/subdiv/RemapFaces.h>
#include <nvmesh/import/MeshImport.h>
#include <nvmesh/render/MeshOptimizer.h>

using namespace	nv;


enum InputMode
{
	InputMode_Triangle,
	InputMode_WingedTriangle,
	InputMode_Quad,
	InputMode_Patch,
};

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

// Primitives per batch in Fermi.
static int getFermiPrimitivePerBatch(uint size)
{
	nvCheck (size <= 32);

	if (size <= 3) return 32;
	if (size == 4) return 24;
	if (size <= 6) return 16;
	if (size <= 15) return 5;
	if (size <= 24) return 4;
	if (size <= 32) return 1;
	return 0;
}


int main(int argc, const char *	argv[])
{
	MyMessageHandler messageHandler;
	debug::setMessageHandler(&messageHandler);

	const char * fileName =	NULL;
	int	inputMode =	InputMode_Triangle;
	VertexCache::Mode cacheMode = VertexCache::Mode_Batch;
	int cacheSize = 32;
	MeshOptimizer::Method optimizationMethod = MeshOptimizer::Method_FermiBatcher;
	
	for	(int i = 1;	i <	argc; i++)
	{
		if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--inputMode") == 0)
		{
			inputMode =	atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--cacheMode") == 0)
		{
			cacheMode = (VertexCache::Mode)atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-n") == 0)
		{
			cacheSize = atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--optimize") == 0)
		{
			optimizationMethod = (MeshOptimizer::Method)atoi(argv[i+1]);
			i++;
		}
		else if	(argv[i][0]	!= '-')
		{
			fileName = argv[i];
		}
	}

	if (fileName ==	NULL)
	{
		printf("Usage: nvmeshbatcher [options] filename\n");
		printf("options are:\n");
		printf("\t-i, --inputMode: Specify input mode\n");
		printf("\t\t0 =	triangle [default]\n");
		printf("\t\t1 =	winged triangle\n");
		printf("\t\t2 =	quad\n");
		printf("\t\t3 =	patches\n");
		printf("\t-o, --optimize: Optimize input\n");
		printf("\t\t0 = Forsyth\n");
		printf("\t\t1 = K-Cache-Reorder\n");
		printf("\t\t2 = Tipsy\n");
		printf("\t\t3 = Castano [default]\n");
		printf("\t\t4 = NVTriStrip\n");
		printf("\t-c, --cacheMode: Cache mode\n");
		printf("\t\t0 =	batch [default]\n");
		printf("\t\t1 =	FIFO\n");
		printf("\t\t2 =	LRU\n");
		printf("\t-n: Cache size [default=32]\n");
		return 0;
	}

	AutoPtr<MeshImport>	importer(MeshImport::importer(fileName));

	if (importer ==	NULL) {
		printf("Error, unkown file type	'%s'\n", fileName);
		return 0;
	}

	StdInputStream stream(fileName);
	if (stream.isError()) {
		printf("Error, cannot open '%s'\n",	fileName);
		return 0;
	}

	importer->import(&stream);

	AutoPtr<HalfEdge::Mesh> mesh(importer->builder().buildHalfEdgeMesh());

	const uint vertexCount = mesh->vertexCount();


	Array<uint>	indexArray;
	uint primitiveSize = 0;

	if (inputMode == InputMode_Triangle)
	{
		//mesh = HalfEdge::triangulate(mesh);

		const uint faceCount = mesh->faceCount();
		const uint vertexCount = mesh->vertexCount();

		// Output basic	info:
		printf("Input mesh:\n");
		printf("    Face count = %u\n",	faceCount);
		printf("    Vertex count = %u\n", vertexCount);
		printf("    Optimal ACMR = %.3f\n", float(vertexCount) / faceCount);

		indexArray.resize(3	* faceCount);

		uint f = 0;
		for	(HalfEdge::Mesh::ConstFaceIterator it(mesh->faces()); !it.isDone(); it.advance(), f++)
		{
			const HalfEdge::Face * face = it.current();

			uint i = 0;
			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++)
			{
				const HalfEdge::Vertex * vertex = it.current()->vertex();
				indexArray[3*f + i + 0] = vertex->id();
			}
			nvCheck(i == 3);
		}

		primitiveSize = 3;
	}
	else if (inputMode == InputMode_WingedTriangle)
	{
		//mesh = HalfEdge::triangulate(mesh); // @@ We assume the input is triangulated already.

		const uint faceCount = mesh->faceCount();

		// Output basic	info:
		printf("Input mesh:\n");
		printf("    Face count = %u\n", faceCount);
		printf("    Vertex count = %u\n", vertexCount);
		printf("    Optimal ACMR = %.3f\n", float(vertexCount) / faceCount);

		indexArray.resize(6 * faceCount);

		uint f = 0;
		for (HalfEdge::Mesh::ConstFaceIterator it(mesh->faces()); !it.isDone(); it.advance(), f++)
		{
			const HalfEdge::Face * face = it.current();

			uint i = 0;
			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i+=2)
			{
				const HalfEdge::Edge * edge = it.current();

				indexArray[6*f + i + 0] = edge->vertex()->id();
				if (edge->isBoundary())
				{
					indexArray[6*f + i + 1] = indexArray[6*f + i + 0];
				}
				else
				{
					indexArray[6*f + i + 1] = edge->pair()->prev()->from()->id();
					//indexArray[6*f + i + 1] = edge->pair()->next()->to()->id();
				}
			}
			nvCheck(i == 6);
		}

		primitiveSize = 6;
	}
	else if (inputMode == InputMode_Quad)
	{
		// @@ Eliminate or triangulate polygons.

		const uint faceCount = mesh->faceCount();
		const uint vertexCount = mesh->vertexCount();

		uint triFaceCount = 0;
		uint quadFaceCount = 0;
		for	(HalfEdge::Mesh::ConstFaceIterator it(mesh->faces()); !it.isDone(); it.advance())
		{
			const HalfEdge::Face * face = it.current();

			const uint edgeCount = face->edgeCount();
			if (edgeCount == 3) triFaceCount++;
			else if (edgeCount == 4) quadFaceCount++;
		}

		// Output basic	info:
		printf("Input mesh:\n");
		printf("    Quad face count = %u\n", quadFaceCount);
		printf("    Vertex count = %u\n", vertexCount);
		printf("    Optimal ACMR = %.3f\n", float(vertexCount) / quadFaceCount);

		indexArray.resize(4 * quadFaceCount);

		uint f = 0;
		for	(HalfEdge::Mesh::ConstFaceIterator it(mesh->faces()); !it.isDone(); it.advance())
		{
			const HalfEdge::Face * face = it.current();

			const uint edgeCount = face->edgeCount();
			if (edgeCount != 4) continue; // @@ Only quads for now.

			uint i = 0;
			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++)
			{
				const HalfEdge::Edge * edge = it.current();

				indexArray[4*f + i] = edge->vertex()->id();
			}
			nvCheck(i == 4);

			f++;
		}

		primitiveSize = 4;
	}
	else if (inputMode == InputMode_Patch)
	{
		RemapFaces::minimizeTopologyCount(mesh.ptr());

		AccMeshBuilder builder(mesh.ptr());

		AccMesh * accMesh = builder.buildAccMesh(/*BuildFlags_DoNotSortPatches*/);

		const uint patchCount = accMesh->patchCount();
		const uint vertexCount = mesh->vertexCount();

		// Output basic	info:
		printf("Input mesh:\n");
		printf("    Patch count = %u\n", patchCount);
		printf("    Vertex count = %u\n", vertexCount);
		printf("    Optimal ACMR = %.3f\n", float(vertexCount) / patchCount);


		for	(uint p = 0; p < patchCount; p++)
		{
			const FaceOneRing * faceOneRing = accMesh->patchAt(p).faceOneRing;

			uint patchVertexCount = faceOneRing->vertexCount();
			if (patchVertexCount > 32)
			{
				printf("Warning: patch %u has %u vertices.\n", p, patchVertexCount);
				patchVertexCount = 32;
			}

			primitiveSize = max(primitiveSize, patchVertexCount);
		}

		indexArray.resize(primitiveSize * patchCount);

		for	(uint p = 0; p < patchCount; p++)
		{
			const FaceOneRing * faceOneRing = accMesh->patchAt(p).faceOneRing;
			const uint patchVertexCount = faceOneRing->vertexCount();

			uint v = 0;
			for (v = 0; v < patchVertexCount; v++)
			{
				indexArray[primitiveSize*p + v] = faceOneRing->vertexAt(v)->id();
			}
			for (; v < primitiveSize; v++)
			{
				indexArray[primitiveSize*p + v] = indexArray[primitiveSize*p + patchVertexCount - 1];
			}
		}
	}
	else
	{
		printf("Error: Input mode not supported.\n");
		return 1;
	}

	const uint primitivePerBatch = getFermiPrimitivePerBatch(primitiveSize);

	const uint primitiveCount	= indexArray.size() / primitiveSize;
	VertexCache::Stats stats = MeshOptimizer::processIndices(indexArray, primitiveSize, cacheSize, cacheMode, primitivePerBatch);

	printf("    Primitive size = %u\n", primitiveSize);
	printf("    Primitive per batch = %u\n", primitivePerBatch);

	printf("Input result:\n");
	printf("    Transform count = %u\n", stats.transformCount);
	printf("    ACMR = %.3f\n", float(stats.transformCount) / primitiveCount);
	printf("    Transform per vertex = %.3f\n", float(stats.transformCount) / vertexCount);
	printf("    Thread per vertex = %.3f\n", float(stats.batchCount * cacheSize) / vertexCount);
	printf("    Batch count = %d\n", stats.batchCount);
	printf("    Full batch count = %d\n", stats.fullBatchCount);


	MeshOptimizer batcher(vertexCount, indexArray, primitiveSize, cacheSize, cacheMode);

	batcher.optimize(optimizationMethod, &indexArray, NULL);

	stats = MeshOptimizer::processIndices(indexArray, primitiveSize, cacheSize, cacheMode, primitivePerBatch);

	printf("Optimized result:\n");
	printf("    Transform count = %u\n", stats.transformCount);
	printf("    ACMR = %.3f\n", float(stats.transformCount) / primitiveCount);
	printf("    Transform per vertex = %.3f\n", float(stats.transformCount) / vertexCount);
	printf("    Thread per vertex = %.3f\n", float(stats.batchCount * cacheSize) / vertexCount);
	printf("    Batch count = %d\n", stats.batchCount);
	printf("    Full batch count = %d\n", stats.fullBatchCount);

	return 0;
}
