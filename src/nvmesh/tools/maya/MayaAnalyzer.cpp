
#include <nvcore/Ptr.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/subdiv/RemapFaces.h>

// Required by maya!
#include <iostream>

#include <maya/MFnPlugin.h> 
#include <maya/MPxCommand.h> 
#include <maya/MGlobal.h>
#include <maya/MString.h>
#include <maya/MDagPath.h>
#include <maya/MFnDagNode.h>
#include <maya/MSelectionList.h>
#include <maya/MItSelectionList.h>
#include <maya/MIOStream.h>
#include <maya/MColor.h>
#include <maya/MFnMesh.h>
#include <maya/MArgList.h>

#include "MayaUtils.h"

using namespace nv;



static void analyzePatchConfigurations(const HalfEdge::Mesh * mesh, const MDagPath & dagPath, bool colorize)
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

	if (colorize)
	{
		MStatus status;
		MFnMesh meshFn(dagPath, &status);

		// Compute max value:
		uint maxValue = 0;
		foreach(i, map)
		{
			maxValue = max(maxValue, map[i].value);
		}
		if (maxValue < faceCount / 3) maxValue = faceCount / 3;


		// Set face colors:
		for(uint f = 0; f < faceCount; f++)
		{
			const HalfEdge::Face * face = mesh->faceAt(f);
			uint key = RemapFaces::topologyId(face);

			uint value = 0;
			map.get(key, &value); 

			float hue = 120.0f * float(value) / float(maxValue);	// hue is in [0-360) range
			MColor color = MColor(MColor::kHSV, hue, 1.0f, 1.0f);

			meshFn.setFaceColor(color, f);
		}
	}
}



class nvAnalyze : public MPxCommand
{ 
public: 
    MStatus doIt( const MArgList & args ); 
    static void * creator() { return new nvAnalyze; }  
}; 

MStatus nvAnalyze::doIt( const MArgList & args ) 
{
	bool colorize = false;

	// Parse the arguments.
	for (uint i = 0; i < args.length(); i++)
	{
		if (MString("-c") == args.asString(i) || MString("-color") == args.asString(i)) 
		{
			colorize = true;
		}
	}

	MDagPath node;
	MObject component;
	MFnDagNode nodeFn;
	

	MSelectionList list; 
	MGlobal::getActiveSelectionList(list); 

	MItSelectionList iter(list, MFn::kMesh);
    if (iter.isDone())
	{
		printf("Error: Nothing is selected.\n");
        return MS::kFailure;
    }

	for (; !iter.isDone(); iter.next()) 
	{ 
		iter.getDagPath(node, component);
        nodeFn.setObject(node);
		
		AutoPtr<HalfEdge::Mesh> mesh(MayaUtils::buildHalfEdgeMesh(node));

		RemapFaces::minimizeTopologyCount(mesh.ptr());

		printf("Analyzing '%s':\n", nodeFn.name().asChar());

		analyzePatchConfigurations(mesh.ptr(), node, colorize);
	}

	return MS::kSuccess; 
}


DLL_EXPORT MStatus initializePlugin(MObject obj)
{
	MFnPlugin plugin(obj, "NVIDIA Corp.", "1.0", "Any"); 
	plugin.registerCommand("nvAnalyze", nvAnalyze::creator);
	return MS::kSuccess; 
}

DLL_EXPORT MStatus uninitializePlugin(MObject obj)
{
	MFnPlugin plugin( obj );
	plugin.deregisterCommand("nvAnalyze");
	return MS::kSuccess;
}
