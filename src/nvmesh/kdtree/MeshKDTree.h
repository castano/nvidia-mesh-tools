// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#ifndef NV_MESH_MESHKDTREE_H
#define NV_MESH_MESHKDTREE_H

#include <nvmesh/nvmesh.h>

namespace nv
{
	class KDTree;
	class TriMesh;
	class QuadTriMesh;

	KDTree * buildKDTree(const TriMesh * mesh);
	KDTree * buildKDTree(const QuadTriMesh * mesh);

} // nv namespace

#endif // NV_MESH_MESHKDTREE_H
