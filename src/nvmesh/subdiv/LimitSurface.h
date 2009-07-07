// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_LIMITSURFACE_H
#define NV_MESH_LIMITSURFACE_H

#include <nvmath/Vector.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
//#include <nvmesh/subdiv/Stencil.h>

namespace nv
{
	namespace CatmullClark
	{
		
		void projectToLimitSurface(HalfEdgeMesh * mesh);
		
		Vector3 limitPosition(const HalfEdgeMesh::Vertex * vertex);
		Vector3 limitNormal(const HalfEdgeMesh::Vertex * vertex);
		Vector3 limitTangent(const HalfEdgeMesh::Edge * edge);
		
		//Stencil * limitPositionStencil(const HalfEdgeMesh::Vertex * vertex);
		//Stencil * limitTangentStencil(const HalfEdgeMesh::Vertex * edge);
		
	} // CatmullClark namespace

	namespace Loop
	{
		
		void projectToLimitSurface(HalfEdgeMesh * mesh);
		
		Vector3 limitPosition(const HalfEdgeMesh::Vertex * vertex);
		Vector3 limitNormal(const HalfEdgeMesh::Vertex * vertex);
		Vector3 limitTangent(const HalfEdgeMesh::Edge * edge);
		
		//Stencil * limitPositionStencil(const HalfEdgeMesh::Vertex * vertex);
		//Stencil * limitTangentStencil(const HalfEdgeMesh::Vertex * edge);
		
	} // Loop namespace

} // nv namespace

#endif // NV_MESH_LIMITSURFACE_H
