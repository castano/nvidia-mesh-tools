// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include "SingleFaceMap.h"

using namespace nv;


Vector3 faceNormal(const HalfEdge::Face * face)
{
	// @@ Compute best fit plane?

	// For now we just use the normal of the first triangle.
	Vector3 p0 = face->edge()->from()->pos();
	Vector3 p1 = face->edge()->to()->pos();
	Vector3 p2 = face->edge()->next()->to()->pos();

	return normalize(cross(p1 - p0, p2 - p0));
}


void nv::computeSingleFaceMap(HalfEdge::Mesh * mesh)
{
	nvDebugCheck(mesh != NULL);
	nvDebugCheck(mesh->faceCount() == 1);

	HalfEdge::Face * face = mesh->faceAt(0);
	nvCheck(face != NULL);

	Vector3 p0 = face->edge()->from()->pos();
	Vector3 p1 = face->edge()->to()->pos();

	Vector3 X = normalize(p1 - p0);
	Vector3 Z = faceNormal(face);
	Vector3 Y = normalize(cross(Z, X));

	uint i = 0;
	for (HalfEdge::Face::EdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++)
	{
		HalfEdge::Vertex * vertex = it.vertex();
		nvCheck(vertex != NULL);

		if (i == 0)
		{
			vertex->setTex(Vector2(0, 0));
		}
		else
		{
			Vector3 pn = vertex->pos();

			float xn = dot((pn - p0), X);
			float yn = dot((pn - p0), Y);

			vertex->setTex(Vector2(xn, yn));
		}
	}
}

