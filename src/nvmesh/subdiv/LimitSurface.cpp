// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_LIMITSURFACE_H
#define NV_MESH_LIMITSURFACE_H

#include <nvmesh/subdiv/LimitSurface.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdgeEdge.h>

using namespace nv;


void CatmullClark::projectToLimitSurface(HalfEdgeMesh * mesh)
{
}

Vector3 CatmullClark::limitPosition(const HalfEdgeMesh::Vertex * vertex)
{
	return Vector3(zero);
}

Vector3 CatmullClark::limitNormal(const HalfEdgeMesh::Vertex * vertex)
{
	return Vector3(zero);
}

Vector3 CatmullClark::limitTangent(const HalfEdgeMesh::Edge * edge)
{
	return Vector3(zero);
}


void Loop::projectToLimitSurface(HalfEdgeMesh * mesh)
{
}

Vector3 Loop::limitPosition(const HalfEdgeMesh::Vertex * vertex)
{
	return Vector3(zero);
}

Vector3 Loop::limitNormal(const HalfEdgeMesh::Vertex * vertex)
{
	return Vector3(zero);
}

Vector3 Loop::limitTangent(const HalfEdgeMesh::Edge * edge)
{
	return Vector3(zero);
}

#endif // NV_MESH_SUBDIVIDE_H
