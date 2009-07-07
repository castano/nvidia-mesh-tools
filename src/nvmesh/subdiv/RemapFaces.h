// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_REMAPFACES_H
#define NV_MESH_REMAPFACES_H

#include <nvmesh/nvmesh.h>

namespace nv
{
	namespace HalfEdge 
	{ 
		class Mesh; 
		class Face;
	}

	namespace RemapFaces
	{
		// Remap faces to avoid parametrization discontinuities.
		void consistentPatchOrientation(HalfEdge::Mesh * mesh);
		
		uint topologyId(const HalfEdge::Face * face);
		void minimizeTopologyCount(HalfEdge::Mesh * mesh);
	}

} // nv namespace

#endif // NV_MESH_REMAPFACES_H
