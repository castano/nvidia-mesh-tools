// This code is in the public domain -- castanyo@yahoo.es

#include "MeshAdjacency.h"

using namespace nv;

void TriMeshAdjacency::buildAdjacency(TriMesh * mesh)
{
	nvCheck(mesh != NULL);
	
	// Check out Eric Lengyel's code.
	
}
	
	
// Proximal methods.
void TriMeshAdjacency::computeProximals()
{
}

void TriMeshAdjacency::addProximal(uint r, uint v)
{
}

bool TriMeshAdjacency::hasProximal(uint v) const
{
	return false;
}

uint TriMeshAdjacency::countProximals(uint v) const
{
	return 1;
}

uint TriMeshAdjacency::firstProximal(uint v) const
{
	return v;
}

