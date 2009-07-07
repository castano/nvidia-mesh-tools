// This code is in the public domain -- castanyo@yahoo.es

#include <stdio.h>

#include <nvcore/StdStream.h>
#include <nvcore/Ptr.h>

#include <nvmesh/TriMesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>
#include <nvmesh/import/MeshImportOBJ.h>
#include <nvmesh/MeshTopology.h>
#include <nvmesh/subdiv/Subdivide.h>


using namespace nv;


void analyze(HalfEdge::Mesh * mesh)
{
	printf("\nBuilding topology:\n");
	MeshTopology topology(mesh);
	
	printf("is closed = %s\n", topology.isClosed() ? "true" : "false");
	printf("is connected = %s\n", topology.isConnected() ? "true" : "false");
	
	if (!topology.isConnected()) {
		printf("connected components = %d\n", topology.connectedCount());
	}
	
	if (topology.isClosed() && topology.isConnected()) {
		printf("genus = %d\n", topology.genus());
		printf("euler number = %d\n", topology.euler());
	}
	
	if (!topology.isClosed()) {
		printf("holes = %d\n", topology.holeCount());
	}


	HashMap<uint, uint> valenceMap;
	HashMap<uint, uint> primitiveMap;
	HashMap<uint, uint> combinedMap;
	
	uint vertexCount = mesh->vertexCount();
	for(uint v = 0; v < vertexCount; v++) {
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		
		if (vertex->isBoundary()) continue;

		uint valence = vertex->valence();
		uint num;
		if (valenceMap.get(valence, &num)) {
			valenceMap.set(valence, num + 1);
		}
		else {
			valenceMap.add(valence, 1);
		}
	}
	
	uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++) {
		HalfEdge::Face * face = mesh->faceAt(f);
		
		if (face->isBoundary()) continue;

		uint count = face->edgeCount();
		uint num;
		if (primitiveMap.get(count, &num)) {
			primitiveMap.set(count, num + 1);
		}
		else {
			primitiveMap.add(count, 1);
		}
		
		// combinations of primitive/valencies.
		for(HalfEdge::Face::EdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			uint valence = it.current()->from()->valence();
			
			uint key = valence << 8 | count;
			uint num;
			if (combinedMap.get(key, &num)) {
				combinedMap.set(key, num + 1);
			}
			else {
				combinedMap.add(key, 1);
			}
		}
	}
	
	printf("valencies = %d\n", valenceMap.count());
	printf("primitives = %d\n", primitiveMap.count());
	printf("split combinations = %d\n", combinedMap.count());
	combinedMap.clear();
	
	// Analyze the number of combinations without the CC split.
	for(uint f = 0; f < faceCount; f++) {
		HalfEdge::Face * face = mesh->faceAt(f);

		if (face->isBoundary()) continue;

		uint min = 256;
		HalfEdge::Edge * start = NULL;

		for(HalfEdge::Face::EdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			HalfEdge::Vertex * vertex = it.current()->from();
			uint valence = vertex->valence();

			if (min > valence) {
				min = valence;
				start = it.current();
			}
		}

		nvCheck(start != NULL && min < 256);

		uint key = 0;
		for(HalfEdge::Face::EdgeIterator it(face->edges(start)); !it.isDone(); it.advance())
		{
			HalfEdge::Vertex * vertex = it.current()->from();
			uint valence = vertex->valence();

			key = key << 8 | valence;
		}

		uint num;
		if (combinedMap.get(key, &num)) {
			combinedMap.set(key, num + 1);
		}
		else {
			combinedMap.add(key, 1);
		}
	}

	printf("nosplit combinations = %d\n", combinedMap.count());
}

void split(HalfEdge::Mesh * mesh)
{
	printf("\nSplitting mesh:\n");
	AutoPtr<HalfEdge::Mesh> splitMesh(Subdivide::doCatmullClarkSplit(mesh));
	
	printf("Face count = %u\n", splitMesh->faceCount());
	printf("Vertex count = %u\n", splitMesh->vertexCount());
	printf("Edge count = %u\n", splitMesh->edgeCount());

	analyze(splitMesh.ptr());
}


void analyzeIrregularFaces(HalfEdge::Mesh * mesh)
{
	uint boundaryFaceCount = 0;
	uint regularFaceCount = 0;

	uint oneBranch = 0;
	uint twoBranch = 0;
	uint threeBranch = 0;

	uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++) {
		HalfEdge::Face * face = mesh->faceAt(f);
		
		if (face->isBoundary()) {
			boundaryFaceCount++;
		}

		bool regular = true;

		uint valenceTable[12] = {
			0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0
		};

		for(HalfEdge::Face::EdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			HalfEdge::Vertex * vertex = it.current()->from();
			uint valence = vertex->valence();

			if (valence != 4) {
				regular = false;
			}

			if (valence < 12) {
				valenceTable[valence]++;
			}
		}

		int branchCount = 0;
		for(int i = 0; i < 12; i++) {
			if (valenceTable[i] != 0) branchCount++;
		}

		if (branchCount == 2) {
			oneBranch++;
		}
		else if (branchCount == 3) {
			twoBranch++;
		}
		else if (branchCount == 4) {
			threeBranch++;
		}

		if (regular)
		{
			regularFaceCount++;
		}
	}

	printf("boundary face count = %d\n", boundaryFaceCount);
	printf("regular face count = %d\n", regularFaceCount);
	printf("irregular face count = %d (%d/%d/%d)\n", faceCount-regularFaceCount, oneBranch, twoBranch, threeBranch);

}

void analyzeTopologyCombinations(HalfEdge::Mesh * mesh)
{
	uint boundaryFaceCount = 0;
	uint regularFaceCount = 0;

	HashMap<uint, uint> map;

	uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++) {
		HalfEdge::Face * face = mesh->faceAt(f);
		
		if (face->isBoundary()) {
			boundaryFaceCount++;
		}

		uint valenceTable[4] = {0, 0, 0, 0};

		int i = 0;
		for(HalfEdge::Face::EdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			HalfEdge::Vertex * vertex = it.current()->from();
			uint valence = vertex->valence();

			if (vertex->isBoundary()) valence = (~valence + 1) & 0xFF;

			valenceTable[i] = valence;
			i++;
		}

		if (valenceTable[0] == 4 && valenceTable[1] == 4 && valenceTable[2] == 4 && valenceTable[3] == 4)
		{
			regularFaceCount++;
		}

		if (valenceTable[0] > valenceTable[1]) swap(valenceTable[0], valenceTable[1]);
		if (valenceTable[0] > valenceTable[2]) swap(valenceTable[0], valenceTable[2]);
		if (valenceTable[0] > valenceTable[3]) swap(valenceTable[0], valenceTable[3]);
		if (valenceTable[1] > valenceTable[2]) swap(valenceTable[1], valenceTable[2]);
		if (valenceTable[1] > valenceTable[3]) swap(valenceTable[1], valenceTable[3]);
		if (valenceTable[2] > valenceTable[3]) swap(valenceTable[2], valenceTable[3]);

		uint key = (valenceTable[0] << 24) | (valenceTable[1] << 16) | (valenceTable[2] << 8) | (valenceTable[3] << 0);
		
		uint value = 0;
		map.get(key, &value); 
		map.set(key, value + 1);
	}

	printf("face count = %d\n", faceCount);
	printf("boundary face count = %d\n", boundaryFaceCount);
	printf("regular face count = %d\n", regularFaceCount);
	printf("topology combination count = %d\n", map.count());

	foreach(i, map)
	{
		char valence0 = (map[i].key >> 24) & 0xFF;
		char valence1 = (map[i].key >> 16) & 0xFF;
		char valence2 = (map[i].key >> 8) & 0xFF;
		char valence3 = (map[i].key >> 0) & 0xFF;

		printf("%X : %d %d %d %d : %u\n", map[i].key, (int)valence0, (int)valence1, (int)valence2, (int)valence3, map[i].value);
	}
}

void outputStencils(HalfEdge::Mesh * mesh)
{

//	printf("corner 0:\n", );
//	printf("\tvalence = %d\n", );

}



void findCompressionRatio(HalfEdge::Mesh * mesh)
{
	// pos, nor, uv -> (3 + 3 + 2) * sizeof(float) = 8 * 4 = 32 bytes
	const int vertexSize = 32;

	// pos, nor, uv, tan -> (3 + 3 + 2 + 4) * sizeof(float) = 12 * 4 = 48 bytes
	//const int vertexSize = 48;
	
	const int baseMeshSize = mesh->vertexCount() * vertexSize;
	
	// pos = 4x4 * sizeof(float3) = 16 * 12 = 192
	// tan = 3x4 * 2 * sizeof(float3) = 12 * 2 * 12 = 288
	// tex = 2x2 * sizeof(float2) = 4 * 8 = 32
	// total = 192 + 288 + 32
	const int controlFaceSize = 512;
	
	// 16 bit displacements.
	const int displacementSize = 2;
	
	const int refinementLevels[] = {4, 8, 16, 32};
	
	printf("\n\n");
	
	for(int i = 0; i < 4; i++) {
		const int refinementLevel = refinementLevels[i];
		const int refinementLevelSq = refinementLevel*refinementLevel;

		printf("tessellation level %d\n", refinementLevel);
		
		const int precomputedSize = vertexSize * refinementLevelSq * mesh->faceCount();
		printf("   precomputed =  %f\n", float(precomputedSize)/1024);
		
		const int displacedFaceSize = refinementLevelSq * displacementSize;
		
		const int bezierMeshSize = (controlFaceSize + displacedFaceSize) * mesh->faceCount();
		printf("   bezier =  %f\n", float(bezierMeshSize)/1024);
		
		const int fullMeshSize = baseMeshSize + displacedFaceSize * mesh->faceCount();
		printf("   full =  %f\n", float(fullMeshSize)/1024);
	}
	
}


int main(int argc, const char * argv[])
{
	if (argc != 2) {
		printf("Usage: nvmeshtest filename\n");
		return 0;
	}

	const char * fileName = argv[1];

	StdInputStream stream(fileName);
	if (stream.isError()) {
		printf("Error, cannot open '%s'\n", fileName);
		return 0;
	}
	
	printf("Importing '%s'\n", fileName);
	MeshImportOBJ importer;
	importer.import(&stream);
	
/*
	printf("\nBuilding triangle mesh:\n");
	AutoPtr<TriMesh> triMesh(importer.builder().buildTriMesh());

	printf("Triangle count = %u\n", triMesh->faceCount());
	printf("Vertex count = %u\n", triMesh->vertexCount());
*/
	// @@ Build tri adjacency.

	printf("\nBuilding half edge mesh:\n");
	AutoPtr<HalfEdge::Mesh> halfEdgeMesh(importer.builder().buildHalfEdgeMesh());

	// @@ Compute holes.

	printf("Face count = %u\n", halfEdgeMesh->faceCount());
	printf("Vertex count = %u\n", halfEdgeMesh->vertexCount());
	printf("Edge count = %u\n", halfEdgeMesh->edgeCount());
/*
	analyze(halfEdgeMesh.ptr());

	split(halfEdgeMesh.ptr());
*/
	analyzeIrregularFaces(halfEdgeMesh.ptr());

	analyzeTopologyCombinations(halfEdgeMesh.ptr());

	//findCompressionRatio(halfEdgeMesh.ptr());

	return 0;
}
