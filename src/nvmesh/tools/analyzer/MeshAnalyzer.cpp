// This code is in the public domain -- castanyo@yahoo.es

#include <stdio.h>

#include <nvcore/StdStream.h>
#include <nvcore/Ptr.h>

#include <nvmesh/TriMesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>
#include <nvmesh/import/MeshImport.h>
#include <nvmesh/MeshTopology.h>
#include <nvmesh/subdiv/Subdivide.h>
#include <nvmesh/subdiv/RemapFaces.h>
#include <nvmesh/param/ParameterizationQuality.h>
#include <nvmesh/param/Atlas.h>
#include <nvmesh/param/LeastSquaresConformalMap.h>
#include <nvmesh/param/SingleFaceMap.h>


using namespace nv;


static void analyzeTopology(const HalfEdge::Mesh * mesh)
{
	printf("\nTopology Info:\n");
	MeshTopology topology(mesh);
	
	printf("  is closed = %s\n", topology.isClosed() ? "true" : "false");
	if (!topology.isClosed()) {
		printf("  holes = %d\n", topology.holeCount());
	}

	printf("  is connected = %s\n", topology.isConnected() ? "true" : "false");
	if (!topology.isConnected()) {
		printf("  connected components = %d\n", topology.connectedCount());
	}
	
	if (topology.isClosed() && topology.isConnected()) {
		printf("  genus = %d\n", topology.genus());
		printf("  euler number = %d\n", topology.euler());
	}
}


static void analyzePatchConfigurations(const HalfEdge::Mesh * mesh)
{
	uint boundaryFaceCount = 0;
	uint regularFaceCount = 0;
	uint semiRegularFaceCount = 0;
	uint regularBoundaryFaceCount = 0;

	uint triangleFaceCount = 0;
	uint quadFaceCount = 0;
	uint polygonFaceCount = 0;

	HashMap<uint, uint> map;

	const uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = mesh->faceAt(f);
		
		if (face->isBoundary()) {
			boundaryFaceCount++;
		}
		
		const uint edgeCount = face->edgeCount();

		if (edgeCount == 3)
		{
			triangleFaceCount++;
		}
		else if (edgeCount == 4)
		{
			quadFaceCount++;
		}
		else
		{
			// Skip polygons.
			polygonFaceCount++;
			continue;
		}

		uint key = RemapFaces::topologyId(face);
		
		if (key == 0x04040404)
		{
			regularFaceCount++;
		}
		else if ((key & 0x00FFFFFF) == 0x00040404 || (key & 0xFFFFFF00) == 0x04040400)
		{
			semiRegularFaceCount++;
		}

		if (key == 0x0404FDFD)
		{
			regularBoundaryFaceCount++;
		}
		
		uint value = 0;
		map.get(key, &value); 
		map.set(key, value + 1);
	}

	printf("\nPatch Configuration Info:\n");
	printf("  patch count = %u\n", faceCount);
	printf("    triangle count = %u (%d%%)\n", triangleFaceCount, (100 * triangleFaceCount / faceCount));
	printf("    quad count = %u (%d%%)\n", quadFaceCount, (100 * quadFaceCount / faceCount));
	printf("    polygon count = %u (%d%%)\n", polygonFaceCount, (100 * polygonFaceCount / faceCount));
	printf("  regular patch count = %u (%d%%)\n", regularFaceCount, (100 * regularFaceCount / faceCount));
	printf("  semi-regular patch count = %u (%d%%)\n", semiRegularFaceCount, (100 * semiRegularFaceCount / faceCount));
	printf("  boundary patch count = %u (%d%%)\n", boundaryFaceCount, (100 * boundaryFaceCount / faceCount));
	printf("  regular boundary patch count = %u (%d%%)\n", regularBoundaryFaceCount, (100 * regularBoundaryFaceCount / faceCount));
	printf("  patch configuration count = %u\n", map.count());

	if (boundaryFaceCount != faceCount)
	{
		printf("\n  Interior Patches:\n");

		foreach(i, map)
		{
			char valence0 = (map[i].key >> 24) & 0xFF;
			char valence1 = (map[i].key >> 16) & 0xFF;
			char valence2 = (map[i].key >> 8) & 0xFF;
			char valence3 = (map[i].key >> 0) & 0xFF;

			if (valence0 >= 0 && valence1 >= 0 && valence2 >= 0 && valence3 >= 0) 
			{
				printf("   - %d %d %d %d : %u\n", (int)valence0, (int)valence1, (int)valence2, (int)valence3, map[i].value);
			}
		}
	}

	if (boundaryFaceCount)
	{
		printf("\n  Boundary Patches:\n");

		foreach(i, map)
		{
			char valence0 = (map[i].key >> 24) & 0xFF;
			char valence1 = (map[i].key >> 16) & 0xFF;
			char valence2 = (map[i].key >> 8) & 0xFF;
			char valence3 = (map[i].key >> 0) & 0xFF;

			if (valence0 < 0 || valence1 < 0 || valence2 < 0 || valence3 < 0) 
			{
				printf("   - %d %d %d %d : %u\n", (int)valence0, (int)valence1, (int)valence2, (int)valence3, map[i].value);
			}
		}
	}
}


static void analyzeExtraordinaryVertices(const HalfEdge::Mesh * mesh)
{
	HashMap<uint, uint> map;

	uint boundaryVertexCount = 0;
	uint irregularVertexCount = 0;
	uint irregularBoundaryVertexCount = 0;

	const uint vertexCount = mesh->vertexCount();
	for(uint v = 0; v < vertexCount; v++) {
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);

		uint valence = vertex->valence();

		if (!vertex->isFirstColocal()) 
		{
			// Process colocal vertices only once.
			continue;
		}

		if (vertex->isBoundary())
		{
			boundaryVertexCount++;

			if (valence != 3) {
				irregularBoundaryVertexCount++;
			}

			// Flag boundary vertices.
			valence = (~valence + 1);
		}
		else
		{
			if (valence != 4) {
				irregularVertexCount++;
			}
		}


		uint value = 0;
		map.get(valence, &value); 
		map.set(valence, value + 1);
	}

	printf("\nExtraordinary Vertex Info:\n");
	printf("  extraordinary vertex count: %u\n", irregularVertexCount);
	printf("  boundary vertex count: %u\n", boundaryVertexCount);
	printf("  boundary extraordinary vertex count: %u\n", irregularBoundaryVertexCount);

	if (boundaryVertexCount != vertexCount)
	{
		printf("\n  Interior Vertices:\n");

		foreach(i, map)
		{
			int valence = map[i].key;
			if (valence >= 0) {
				printf("   + %d : %u\n", valence, map[i].value);
			}
		}
	}

	if (boundaryVertexCount)
	{
		printf("\n  Boundary Vertices:\n");

		foreach(i, map)
		{
			int valence = map[i].key;
			if (valence < 0) {
				printf("   + %d : %u\n", -valence, map[i].value);
			}
		}
	}

}


static void analyzeParameterizationQuality(const HalfEdge::Mesh * mesh)
{
	// Analyze existing parameterization.
	ParameterizationQuality parameterizationQuality(mesh);

	printf("\nParameterization Quality:\n");
	printf("  Valid parameterization: %s\n", parameterizationQuality.isValid() ? "yes" : "no");
	printf("  RMS stretch metric: %f\n", parameterizationQuality.rmsStretchMetric());
	printf("  MAX stretch metric: %f\n", parameterizationQuality.maxStretchMetric());
	printf("  RMS conformal metric: %f\n", parameterizationQuality.rmsConformalMetric());
	printf("  RMS authalic metric: %f\n", parameterizationQuality.maxAuthalicMetric());

	Atlas atlas(mesh);
	atlas.extractCharts();

	const uint chartCount = atlas.chartCount();
	for (uint i = 0; i < chartCount; i++)
	{
		Chart * chart = atlas.chartAt(i);
		nvDebugCheck(chart != NULL);

		// @@ If has multiple holes -> Close holes with "virtual faces".

		if (chart->isDisk())
		{
			if (chart->faceCount() == 1)
			{
				computeSingleFaceMap(chart->mesh());
			}
			else
			{
				computeLeastSquaresConformalMap(chart->mesh());
			}
		}
	}

	// atlas.scaleCharts();
	// atlas.packCharts();

	HalfEdge::Mesh * remappedMesh = atlas.mergeCharts();

	ParameterizationQuality parameterizationQuality2(remappedMesh);
	
	printf("\nParameterization Quality:\n");
	printf("  Valid parameterization: %s\n", parameterizationQuality2.isValid() ? "yes" : "no");
	printf("  RMS stretch metric: %f\n", parameterizationQuality2.rmsStretchMetric());
	printf("  MAX stretch metric: %f\n", parameterizationQuality2.maxStretchMetric());
	printf("  RMS conformal metric: %f\n", parameterizationQuality2.rmsConformalMetric());
	printf("  RMS authalic metric: %f\n", parameterizationQuality2.maxAuthalicMetric());

}



int main(int argc, const char * argv[])
{
	bool split = false;
	const char * fileName = NULL;

	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "--split") == 0)
		{
			split = true;
		}
		else if (argv[i][0] != '-')
		{
			fileName = argv[i];
		}
	}

	if (fileName == NULL) {
		printf("Usage: nvmeshanalyzer [--split] filename\n");
		return 0;
	}

	AutoPtr<MeshImport> importer(MeshImport::importer(fileName));

	if (importer == NULL) {
		printf("Error, unkown file type '%s'\n", fileName);
		return 0;
	}

	StdInputStream stream(fileName);
	if (stream.isError()) {
		printf("Error, cannot open '%s'\n", fileName);
		return 0;
	}

	importer->import(&stream);

	AutoPtr<HalfEdge::Mesh> halfEdgeMesh(importer->builder().buildHalfEdgeMesh());

	RemapFaces::minimizeTopologyCount(halfEdgeMesh.ptr());

	// For some reason, enabling this line, makes gcc produce incorrect code.
	/*if (split)
	{
		halfEdgeMesh = Subdivide::doCatmullClarkSplit(halfEdgeMesh.ptr());
	}*/

	// Output basic info:
	printf("\nBasic Mesh Info:\n");
	printf("  Face count = %u\n", halfEdgeMesh->faceCount());
	printf("  Vertex count = %u\n", halfEdgeMesh->colocalVertexCount());
	printf("  Edge count = %u\n", halfEdgeMesh->edgeCount());

	// Output topology info:
	analyzeTopology(halfEdgeMesh.ptr());

	// Output patch configuration statistics:
	analyzePatchConfigurations(halfEdgeMesh.ptr());
	
	// Output extraordinary vertex statistics:
	analyzeExtraordinaryVertices(halfEdgeMesh.ptr());

	// Output parameterization info:
	analyzeParameterizationQuality(halfEdgeMesh.ptr());

	return 0;
}
