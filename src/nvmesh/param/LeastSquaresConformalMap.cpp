// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmath/Sparse.h>
#include <nvmath/Solver.h>

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>

#include "LeastSquaresConformalMap.h"
#include "ParameterizationQuality.h"

#include <stdio.h>

using namespace nv;

namespace
{

	static uint countMeshTriangles(const HalfEdge::Mesh * mesh)
	{
		const uint faceCount = mesh->faceCount();

		uint triangleCount = 0;

		for (uint f = 0; f < faceCount; f++)
		{
			const HalfEdge::Face * face = mesh->faceAt(f);

			triangleCount += face->edgeCount() - 2;
		}

		return triangleCount;
	}

	static const HalfEdge::Vertex * findBoundaryVertex(const HalfEdge::Mesh * mesh)
	{
		const uint vertexCount = mesh->vertexCount();

		for (uint v = 0; v < vertexCount; v++)
		{
			const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
			if (vertex->isBoundary()) return vertex;
		}

		return NULL;
	}


	// Test all pairs!	
	void findDiameterVertices(HalfEdge::Mesh * mesh, HalfEdge::Vertex ** a, HalfEdge::Vertex ** b)
	{
		nvDebugCheck(mesh != NULL);
		nvDebugCheck(a != NULL);
		nvDebugCheck(b != NULL);
		
		const HalfEdge::Vertex * vertex = findBoundaryVertex(mesh);
		nvDebugCheck(vertex == NULL);

		float maxLength = 0.0f;
		
		*a = *b = NULL;

		const HalfEdge::Edge * firstEdge = vertex->edge();
		nvDebugCheck(firstEdge->face() == NULL);

		const HalfEdge::Edge * edge0 = firstEdge;
		const HalfEdge::Edge * edge1 = firstEdge->next();
		/*
		while (edge1 != firstEdge)
		{
			edge0 = edge0->next();

			if (edge0

		};
		*/
	}


	// Test all pairs and check distance in textrue space.
	void findTextureDiameterVertices(HalfEdge::Mesh * mesh, HalfEdge::Vertex ** a, HalfEdge::Vertex ** b)
	{
		nvDebugCheck(mesh != NULL);
		nvDebugCheck(a != NULL);
		nvDebugCheck(b != NULL);

		const uint vertexCount = mesh->vertexCount();

		float maxLength = 0.0f;

		for (uint v0 = 1; v0 < vertexCount; v0++)
		{
			HalfEdge::Vertex * vertex0 = mesh->vertexAt(v0);
			nvDebugCheck(vertex0 != NULL);
			
			if (!vertex0->isBoundary()) continue;
			
			for (uint v1 = 0; v1 < v0; v1++)
			{
				HalfEdge::Vertex * vertex1 = mesh->vertexAt(v1);
				nvDebugCheck(vertex1 != NULL);

				if (!vertex1->isBoundary()) continue;
				
				float len = length(vertex0->pos() - vertex1->pos());

				if (len > maxLength)
				{
					maxLength = len;
					
					*a = vertex0;
					*b = vertex1;
				}
			}
		}
		
		nvDebugCheck(*a != NULL && *b != NULL);
	}

	// Fast sweep in 3 directions
	void findApproximateDiameterVertices(HalfEdge::Mesh * mesh, HalfEdge::Vertex ** a, HalfEdge::Vertex ** b)
	{
		nvDebugCheck(mesh != NULL);
		nvDebugCheck(a != NULL);
		nvDebugCheck(b != NULL);
		
		const uint vertexCount = mesh->vertexCount();
		
		HalfEdge::Vertex * minVertex[3];
		HalfEdge::Vertex * maxVertex[3];

#pragma message(NV_FILE_LINE "FIXME: Do not initialize min/max vertices with vertex that may not be on a boundary.")
		minVertex[0] = minVertex[1] = minVertex[2] = mesh->vertexAt(0);
		maxVertex[0] = maxVertex[1] = maxVertex[2] = mesh->vertexAt(0);

		for (uint v = 1; v < vertexCount; v++)
		{
			HalfEdge::Vertex * vertex = mesh->vertexAt(v);
			nvDebugCheck(vertex != NULL);

			if (!vertex->isBoundary())
			{
				// Skip interior vertices.
				continue;
			}

			if (vertex->pos().x() < minVertex[0]->pos().x()) minVertex[0] = vertex;
			else if (vertex->pos().x() > maxVertex[0]->pos().x()) maxVertex[0] = vertex;

			if (vertex->pos().y() < minVertex[1]->pos().y()) minVertex[1] = vertex;
			else if (vertex->pos().y() > maxVertex[1]->pos().y()) maxVertex[1] = vertex;

			if (vertex->pos().z() < minVertex[2]->pos().z()) minVertex[2] = vertex;
			else if (vertex->pos().z() > maxVertex[2]->pos().z()) maxVertex[2] = vertex;
		}

		float lengths[3];
		for (int i = 0; i < 3; i++)
		{
			lengths[i] = length(minVertex[i]->pos() - maxVertex[i]->pos());
		}

		if (lengths[0] > lengths[1] && lengths[0] > lengths[2])
		{
			*a = minVertex[0];
			*b = maxVertex[0];
		}
		else if (lengths[1] > lengths[2])
		{
			*a = minVertex[1];
			*b = maxVertex[1];
		}
		else
		{
			*a = minVertex[2];
			*b = maxVertex[2];
		}
	}

	// These two functions are from Bruno Levy and were under GPL!! He would surely grant permission, but I should ask anyway.
	
    // Computes the coordinates of the vertices of a triangle
    // in a local 2D orthonormal basis of the triangle's plane.
    static void project_triangle(Vector3::Arg p0, Vector3::Arg p1, Vector3::Arg p2, Vector2 * z0, Vector2 * z1, Vector2 * z2)
	{
		Vector3 X = normalize(p1 - p0);
		Vector3 Z = normalize(cross(X, (p2 - p0)));
		Vector3 Y = normalize(cross(Z, X));
		
		float x0 = 0.0f;
		float y0 = 0.0f;
		float x1 = length(p1 - p0);
		float y1 = 0.0f;
		float x2 = dot((p2 - p0), X);
		float y2 = dot((p2 - p0), Y);
		
		*z0 = Vector2(x0, y0);
		*z1 = Vector2(x1, y1);
		*z2 = Vector2(x2, y2);
    }

    // LSCM equation, geometric form :
    // (Z1 - Z0)(U2 - U0) = (Z2 - Z0)(U1 - U0)
    // Where Uk = uk + i.vk is the complex number 
    //                       corresponding to (u,v) coords
    //       Zk = xk + i.yk is the complex number 
    //                       corresponding to local (x,y) coords
    // cool: no divide with this expression,
    //  makes it more numerically stable in
    //  the presence of degenerate triangles.
    
    void setup_conformal_map_relations(SparseMatrix & A, int row, const HalfEdge::Vertex * v0, const HalfEdge::Vertex * v1, const HalfEdge::Vertex * v2)
	{
		int id0 = v0->id();
		int id1 = v1->id();
		int id2 = v2->id();
		
		Vector3 p0 = v0->pos();
		Vector3 p1 = v1->pos();
		Vector3 p2 = v2->pos();
		
		Vector2 z0, z1, z2 ;
		project_triangle(p0, p1, p2, &z0, &z1, &z2);
		
		Vector2 z01 = z1 - z0;
		Vector2 z02 = z2 - z0;
		
		float a = z01.x();
		float b = z01.y();
		float c = z02.x();
		float d = z02.y();
		nvCheck(b == 0.0f);

		// Note  : 2*id + 0 --> u
		//         2*id + 1 --> v
		int u0_id = 2 * id0 + 0;
		int v0_id = 2 * id0 + 1;
		int u1_id = 2 * id1 + 0;
		int v1_id = 2 * id1 + 1;
		int u2_id = 2 * id2 + 0;
		int v2_id = 2 * id2 + 1;

		// Note : b = 0

		// Real part
		A.setCoefficient(u0_id, 2 * row + 0, -a+c);
		A.setCoefficient(v0_id, 2 * row + 0,  b-d);
		A.setCoefficient(u1_id, 2 * row + 0,   -c);
		A.setCoefficient(v1_id, 2 * row + 0,    d);
		A.setCoefficient(u2_id, 2 * row + 0,    a);

		// Imaginary part
		A.setCoefficient(u0_id, 2 * row + 1, -b+d);
		A.setCoefficient(v0_id, 2 * row + 1, -a+c);
		A.setCoefficient(u1_id, 2 * row + 1,   -d);
		A.setCoefficient(v1_id, 2 * row + 1,   -c);
		A.setCoefficient(v2_id, 2 * row + 1,    a);
	}
	
} // namespace


void nv::computeLeastSquaresConformalMap(HalfEdge::Mesh * mesh)
{
	nvDebugCheck(mesh != NULL);

	/*{
		ParameterizationQuality parameterizationQuality(mesh);

		printf("BEFORE:\n");
		printf("  Valid parameterization: %s\n", parameterizationQuality.isValid() ? "yes" : "no");
		printf("  RMS stretch metric: %f\n", parameterizationQuality.rmsStretchMetric());
		printf("  MAX stretch metric: %f\n", parameterizationQuality.maxStretchMetric());
		printf("  RMS conformal metric: %f\n", parameterizationQuality.rmsConformalMetric());
		printf("  RMS authalic metric: %f\n", parameterizationQuality.maxAuthalicMetric());
	}*/

	// For this to work properly, mesh should not have colocals that have the same 
	// attributes, unless you want the vertices to actually have different texcoords.

	const uint vertexCount = mesh->vertexCount();
	const uint D = 2 * (vertexCount /*- 2*/);
	const uint N = 2 * countMeshTriangles(mesh);

	// N is the number of equations (one per triangle)
	// D is the number of variables (one per vertex; there are 2 pinned vertices).
	nvDebugCheck(N >= D - 4);

	SparseMatrix A(D, N);
	FullVector b(N);
	FullVector x(D);

	// Fill b:
	b.fill(0.0f);
	
	// Fill x:
	HalfEdge::Vertex * v0;
	HalfEdge::Vertex * v1;
	findApproximateDiameterVertices(mesh, &v0, &v1);
	//findTextureDiameterVertices(mesh, &v0, &v1);

	for (uint v = 0; v < vertexCount; v++)
	{
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		// Initial solution.
		x[2 * v + 0] = vertex->tex().x();
		x[2 * v + 1] = vertex->tex().y();
	}

	// Fill A:
	const uint faceCount = mesh->faceCount();
	for (uint f = 0, t = 0; f < faceCount; f++)
	{
		HalfEdge::Face * face = mesh->faceAt(f);
		nvDebugCheck(face != NULL);
		
		const HalfEdge::Vertex * vertex0 = NULL;

		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			const HalfEdge::Edge * edge = it.current();
			nvCheck(edge != NULL);
			
			if (vertex0 == NULL)
			{
				vertex0 = edge->vertex();
			}
			else if (edge->to() != vertex0)
			{
				const HalfEdge::Vertex * vertex1 = edge->from();
				const HalfEdge::Vertex * vertex2 = edge->to();
				
				setup_conformal_map_relations(A, t, vertex0, vertex1, vertex2);
				
				t++;
			}
		}
	}

	const uint lockedParameters[] =
	{
		2 * v0->id() + 0, 
		2 * v0->id() + 1,
		2 * v1->id() + 0,
		2 * v1->id() + 1
	};

	// Solve
	LeastSquaresSolver(A, b, x, lockedParameters, 4);
	
	// Map x back to texcoords:
	for (uint v = 0; v < vertexCount; v++)
	{
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		Vector2 tc(x[2 * v + 0], x[2 * v + 1]);
		vertex->setTex(tc);
	}

	/*{
		ParameterizationQuality parameterizationQuality(mesh);

		printf("AFTER:\n");
		printf("  Valid parameterization: %s\n", parameterizationQuality.isValid() ? "yes" : "no");
		printf("  RMS stretch metric: %f\n", parameterizationQuality.rmsStretchMetric());
		printf("  MAX stretch metric: %f\n", parameterizationQuality.maxStretchMetric());
		printf("  RMS conformal metric: %f\n", parameterizationQuality.rmsConformalMetric());
		printf("  RMS authalic metric: %f\n", parameterizationQuality.maxAuthalicMetric());
	}*/
}
