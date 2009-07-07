// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include "AccMeshBuilder.h"

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>

#include <nvmesh/subdiv/AccMesh.h>
#include <nvmesh/subdiv/AccPatch.h>
#include <nvmesh/subdiv/RemapFaces.h>

#include <nvmesh/geometry/MeshNormals.h>
#include <nvmesh/geometry/TangentSpace.h>

#include <nvmath/Basis.h>

#include <nvcore/Ptr.h>
#include <nvcore/Radix.h>
#include <nvcore/Containers.h>

using namespace nv;

namespace
{
	static float alphaWeight(int i, int n)
	{
		const float cp = cosf(PI/n);
		const float denom = (n * sqrtf(4 + cp*cp));

		if (n == 4)
		{
			// Special path to ensure symmetry.
			if (i == 0) return (1.0f / n + cp / denom);
			if (i == 2) return -(1.0f / n + cp / denom);
			return 0.0f;
		}
		else
		{
			return (1.0f / n + cp / denom) * cosf((2.0f * PI * i) / n); 
		}
	}

	static float betaWeight(int i, int n)
	{
		const float cp = cosf(PI/n);
		const float denom = (n * sqrtf(4 + cp*cp));

		if (n == 4)
		{
			// Special path to ensure symmetry.
			if (i == 0 || i == 3) return (1.0f / denom) * sqrtf(2.0f) / 2.0f;
			/*if (i == 1 || i == 2)*/ return -(1.0f / denom) * sqrtf(2.0f) / 2.0f;
		}
		else
		{
			return (1.0f / denom) * cosf((2.0f * PI * i + PI) / n);
		}
	}

	static float pseudoValence(const HalfEdge::Vertex * vertex)
	{
		float valence = float(vertex->valence());

		if (vertex->isBoundary())
		{
			// We treat boundary vertices as being half a closed mesh. Corners are special case.
			// n = 4 for corners and n = 2*(n-1) for boundaries.
			if (valence == 2) return 4;
			return (valence - 1) * 2;
		}

		return valence;
	}

} // namespace 

AccMeshBuilder::AccMeshBuilder(const HalfEdge::Mesh * mesh, BoundaryMode bm/*= BoundaryMode_Spline*/) :
	m_mesh(mesh),
	m_boundaryMode(bm)
{
	nvCheck(mesh != NULL);
}

AccMeshBuilder::~AccMeshBuilder()
{
}

AccMesh * AccMeshBuilder::buildAccMesh(uint buildFlags/*= 0*/) const
{
	// Note: best results are obtained if you first call `RemapFaces::minimizeTopologyCount(mesh)`.

	const uint vertexCount = m_mesh->vertexCount();
	const uint faceCount = m_mesh->faceCount();

	// Compute global tangent basis.
	Array<Basis> basisArray;
	computeGlobalTangentSpace(basisArray);

	AutoPtr<AccMesh> accMesh(new AccMesh());
	
	/*for (HalfEdge::Mesh::ConstVertexIterator it(m_mesh->vertices()); !it.isDone(); it.advance())
	{
		const HalfEdge::Vertex * vertex = it.current();
		accMesh->addVertex(vertex);
	}*/

	for (uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = m_mesh->faceAt(f);
		nvDebugCheck(face != NULL);

        const bool isBoundary = FaceOneRing::isBoundary(face);

		// Skip boundaries.
		if (m_boundaryMode == BoundaryMode_None && isBoundary)
		{
			continue;
		}

		const HalfEdge::Edge * firstEdge = NULL;
		FaceTopology * faceTopology = accMesh->addFaceTopology(face, &firstEdge);
		nvCheck (firstEdge != NULL);

		AutoPtr<FaceOneRing> faceOneRing(new FaceOneRing(face, firstEdge, faceTopology));

		bool addQuadPatch = (buildFlags & BuildFlags_GenerateQuadPatches) != 0 && FaceOneRing::isQuad(face);
		bool addRegularPatch = (buildFlags & BuildFlags_GenerateRegularPatches) != 0 && FaceOneRing::isRegular(face);
		bool addTriPatch = (buildFlags & BuildFlags_GenerateTriangularPatches) != 0 && FaceOneRing::isTriangle(face);


        if (addQuadPatch)
        {
			uint p;
            if (addRegularPatch)
            {
				p = accMesh->addRegularPatch(faceOneRing.ptr());

                if (buildFlags & BuildFlags_BuildBezierPatches)
                {
                    AccPatch * patch = buildBezierAccPatch(faceOneRing.ptr());
                    if (patch != NULL) accMesh->setRegularPatch(p, patch);
                }
                if (buildFlags & BuildFlags_BuildGregoryPatches)
                {
                    AccPatch * patch = buildGregoryAccPatch(faceOneRing.ptr());
                    if (patch != NULL) accMesh->setRegularPatch(p, patch);
                }
                if (buildFlags & BuildFlags_BuildPmPatches)
                {
                    AccPatch * patch = buildPmRegularAccPatch(faceOneRing.ptr());
                    if (patch != NULL) accMesh->setRegularPatch(p, patch);
                }
				if (buildFlags & BuildFlags_BuildTexCoordPatches)
				{
					AccPatch * patch = buildTexCoordPatch(faceOneRing.ptr());
					if (patch != NULL) accMesh->setRegularPatch(p, patch);
				}
				if (buildFlags & BuildFlags_BuildTangentSpacePatches)
				{
					//AccPatch * patch = buildTangentSpacePatch(faceOneRing.ptr(), basisArray);
					//if (patch != NULL) accMesh->setRegularPatch(p, patch);
				}
			}
            else
            {
                p = accMesh->addQuadPatch(faceOneRing.ptr());

                if (buildFlags & BuildFlags_BuildBezierPatches)
                {
                    AccPatch * patch = buildBezierAccPatch(faceOneRing.ptr());
                    if (patch != NULL) accMesh->setQuadPatch(p, patch);
                }
                if (buildFlags & BuildFlags_BuildGregoryPatches)
                {
                    AccPatch * patch = buildGregoryAccPatch(faceOneRing.ptr());
                    if (patch != NULL) accMesh->setQuadPatch(p, patch);
                }
                if (buildFlags & BuildFlags_BuildPmPatches)
                {
                    AccPatch * patch = buildPmQuadAccPatch(faceOneRing.ptr());
                    if (patch != NULL) accMesh->setQuadPatch(p, patch);
                }
                if (buildFlags & BuildFlags_BuildTexCoordPatches)
                {
                    AccPatch * patch = buildTexCoordPatch(faceOneRing.ptr());
                    if (patch != NULL) accMesh->setQuadPatch(p, patch);
                }
				if (buildFlags & BuildFlags_BuildTangentSpacePatches)
				{
					//AccPatch * patch = buildTangentSpacePatch(faceOneRing.ptr(), basisArray);
					//if (patch != NULL) accMesh->setQuadPatch(p, patch);
				}
            }

			faceOneRing.release();
        }
        else if (addTriPatch)
        {
            const uint p = accMesh->addTriPatch(faceOneRing.ptr());

            if (buildFlags & BuildFlags_BuildGregoryPatches)
            {
                AccPatch * patch = buildGregoryAccPatch(faceOneRing.ptr());
                if (patch != NULL) accMesh->setTriPatch(p, patch);
            }
            if (buildFlags & BuildFlags_BuildPmPatches)
            {
                AccPatch * patch = buildPmTriangleAccPatch(faceOneRing.ptr());
                if (patch != NULL) accMesh->setTriPatch(p, patch);
            }
			if ((buildFlags & BuildFlags_BuildPolarPatches) && FaceOneRing::isPolarTriangle(face))
			{
                AccPatch * patch = buildPolarPatch(faceOneRing.ptr());
                if (patch != NULL) accMesh->setTriPatch(p, patch);
            }
            if (buildFlags & BuildFlags_BuildTexCoordPatches)
            {
                AccPatch * patch = buildTexCoordPatch(faceOneRing.ptr());
                if (patch != NULL) accMesh->setTriPatch(p, patch);
            }
			if (buildFlags & BuildFlags_BuildTangentSpacePatches)
			{
				//AccPatch * patch = buildTangentSpacePatch(faceOneRing.ptr(), basisArray);
				//if (patch != NULL) accMesh->setTriPatch(p, patch);
			}

			faceOneRing.release();
        }
	}

	printf("Vertex Topology Count = %d\n", accMesh->vertexTopologyCount());
	printf("Face Topology Count = %d\n", accMesh->faceTopologyCount());

	// @@ Create patch groups?

	return accMesh.release();
}



//-------------------------------
//  BezierAccPatch 
//-------------------------------

BezierAccPatch * AccMeshBuilder::buildBezierAccPatch(const FaceOneRing * face) const
{
	nvDebugCheck(face->edgeCount() == 4);

	BezierAccPatch * patch = new BezierAccPatch(face);

	computePositionStencil(patch);
	computeTangentStencil(patch);

	patch->evaluateControlPoints();

	FaceTopology * faceTopology = patch->faceOneRing()->faceTopology();
	if (faceTopology->bezierAccStencil == NULL)
	{
		faceTopology->bezierAccStencil = patch->stencil();
	}
	else
	{
		// Make sure stencils are the same. If this fails this is probably due to self-overlapping stencils.
		nvCheck(*faceTopology->bezierAccStencil == *patch->stencil());

		// @@ Delete new stencils, reference the same.
	}

	return patch;
}

// Positions
void AccMeshBuilder::computePositionStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	// Position control points are in this order:
	//  0 -- 1 -- 2 -- 3
	//  |    |    |    |
	//  4 -- 5 -- 6 -- 7
	//  |    |    |    |
	//  8 -- 9 --10 --11
	//  |    |    |    |
	// 12 --13 --14 --15

	// 0, 3, 15, 12
	computeCornerPositionStencil(patch);

	// 1, 2, 7, 11, 14, 13, 8, 4 
	computeEdgePositionStencil(patch);

	// 5, 6, 10, 9
	computeInteriorPositionStencil(patch);
}

void AccMeshBuilder::computeCornerPositionStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);

	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	nvDebugCheck(firstEdge != NULL);

	const uint cornerIndices[4] = {0, 3, 15, 12};

	BezierAccStencil * stencil = patch->stencil();

	uint v = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Vertex * vertex = it.current()->from();
		const uint cornerIndex = cornerIndices[v];

		if (vertex->isBoundary())
		{
			nvDebugCheck(m_boundaryMode == BoundaryMode_Spline);

			// Compute limit vertex position.
			stencil->positionStencil(cornerIndex, vertex) = 2.0f / 3.0f;

			const HalfEdge::Edge * edge0 = vertex->edge();
			nvCheck(edge0->face() == NULL);
			nvCheck(edge0->to() != vertex);

			stencil->positionStencil(cornerIndex, edge0->to()) = 1.0f / 6.0f;

			const HalfEdge::Edge * edge1 = vertex->edge()->prev();
			nvCheck(edge1->face() == NULL);
			nvCheck(edge1->from() != vertex);

			stencil->positionStencil(cornerIndex, edge1->from()) = 1.0f / 6.0f;

			nvDebugCheck(stencil->positionStencil(cornerIndex).isNormalized());
		}
		else
		{
			const float valence = float(vertex->valence());

			stencil->positionStencil(cornerIndex, vertex) = valence * valence;

			// Traverse one ring consistently for higher accuracy.
			const HalfEdge::Edge * edge0 = it.current();

			for(HalfEdge::Vertex::ConstEdgeIterator eit(edge0); !eit.isDone(); eit.advance())
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(vertex->pos() == edge->from()->pos());

				stencil->positionStencil(cornerIndex, edge->to()) = 4.0f;

				if (edge->next()->next()->next() == edge)
				{
					// Distribute weight to all vertices.
					stencil->positionStencil(cornerIndex, vertex) += 1.0f / 3.0f;
					stencil->positionStencil(cornerIndex, edge->to()) += 1.0f / 3.0f;
					stencil->positionStencil(cornerIndex, edge->next()->to()) += 1.0f / 3.0f;
				}
				else
				{
					stencil->positionStencil(cornerIndex, edge->next()->to()) = 1.0f;
				}
			}

			// Normalize stencil.
			stencil->positionStencil(cornerIndex).normalize();
		}
	}
}

// Compute interior control point stencils.
static void addInteriorPositionStencil(BezierAccPatch * patch, uint idx, const HalfEdge::Face * face, const HalfEdge::Edge * edge)
{
	nvDebugCheck(patch != NULL);
	nvDebugCheck(idx < 16);
	nvDebugCheck(face != NULL);
	nvDebugCheck(edge != NULL);

	const uint edgeCount = face->edgeCount();
	const float valence = pseudoValence(edge->from());

	float weights[4];
	if (edgeCount == 4)
	{
		weights[0] = 3 * valence;
		weights[1] = 6;
		weights[2] = 3;
		weights[3] = 6;
	}
	else
	{
		nvCheck(edgeCount == 3);
		weights[0] = 3 * valence + 1;
		weights[1] = 7;
		weights[2] = 7;
	}

	BezierAccStencil * stencil = patch->stencil();

	uint i = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges(edge)); !it.isDone(); it.advance(), i++)
	{
		const HalfEdge::Vertex * vertex = it.current()->from();
		stencil->positionStencil(idx, vertex) += weights[i];
	}
}

void AccMeshBuilder::computeEdgePositionStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);

	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	nvDebugCheck(firstEdge != NULL);

	BezierAccStencil * stencil = patch->stencil();

	// Compute edge control points.
	const uint edgeIndices[8] = {1, 2, 7, 11, 14, 13, 8, 4};

	uint v = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Edge * edge0 = it.current();
		const HalfEdge::Edge * edge1 = it.current()->pair();

		const HalfEdge::Face * face0 = edge0->face();
		const HalfEdge::Face * face1 = edge1->face();

		const uint index0 = edgeIndices[2 * v + 0];
		const uint index1 = edgeIndices[2 * v + 1];

		nvDebugCheck(face0 == face);
		if (face1 != NULL)
		{
			// Interior boundary control points are the average of the neighbor interior control points.
			addInteriorPositionStencil(patch, index0, face0, edge0);
			addInteriorPositionStencil(patch, index0, face1, edge1->next());

			addInteriorPositionStencil(patch, index1, face0, edge0->next());
			addInteriorPositionStencil(patch, index1, face1, edge1);

			stencil->positionStencil(index0).normalize();
			stencil->positionStencil(index1).normalize();
		}
		else
		{
			nvDebugCheck(m_boundaryMode == BoundaryMode_Spline);

			stencil->positionStencil(index0, edge0->vertex()) = 2.0f / 3.0f;
			stencil->positionStencil(index0, edge1->vertex()) = 1.0f / 3.0f;

			stencil->positionStencil(index1, edge0->vertex()) = 1.0f / 3.0f;
			stencil->positionStencil(index1, edge1->vertex()) = 2.0f / 3.0f;

			nvDebugCheck(stencil->positionStencil(index0).isNormalized());
			nvDebugCheck(stencil->positionStencil(index1).isNormalized());
		}
	}
}

void AccMeshBuilder::computeInteriorPositionStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);

	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	nvDebugCheck(firstEdge != NULL);

	BezierAccStencil * stencil = patch->stencil();

	// Compute interior control points.
	const uint interiorIndices[4] = {5, 6, 10, 9};

	uint v = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Edge * edge = it.current();
		const uint interiorIndex = interiorIndices[v];

		addInteriorPositionStencil(patch, interiorIndex, face, edge);

		stencil->positionStencil(interiorIndex).normalize();
	}
}


// Tangents
void AccMeshBuilder::computeTangentStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	// Tangent control points are in this order:

	//  +---+---+---+
	//  0   1   2   3
	//  +---+---+---+
	//  4   5   6   7
	//  +---+---+---+
	//  8   9  10  11
	//  +---+---+---+

	// 0, 3, 8, 11
	computeCornerTangentStencil(patch);

	// 4, 5, 6, 7
	computeInteriorTangentStencil(patch);

	// 1, 2, 9, 10
	computeEdgeTangentStencil(patch);

	// Symmetric boundary tangent control points are in this order:
	
	//  +-15-+---+-4-+
	//  0    1   2   3
	//  +-14-+---+-5-+
	//  |    |   |   7
	//  +-13-+---+-6-+
	//  11   10  9   8
	//  +-12-+---+-7-+

	//computeSymmetricTangentStencil(patch);
}


void AccMeshBuilder::computeCornerTangentStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);

	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	nvDebugCheck(firstEdge != NULL);

	BezierAccStencil * stencil = patch->stencil();

	static const uint tangentIndices[4] = {0, 8, 11, 3};
	static const uint bitangentIndices[4] = {0, 3, 11, 8};

	static const float tangentFactors[4] = {1, -1, -1, 1};
	static const float bitangentFactors[4] = {1, 1, -1, -1};

	uint v = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Edge * faceEdge = it.current();
		const HalfEdge::Vertex * vertex = faceEdge->from();

		const uint valence = vertex->valence();

		if (vertex->isBoundary())
		{
			const uint vertexCount = patch->vertexCount();
			StencilMask	r0(vertexCount);
			StencilMask	r1(vertexCount);

			computeBoundaryTangentStencils(patch, vertex, r0, r1);

			// Compute rotation factors.
			const HalfEdge::Edge * tanEdge = faceEdge;
			const HalfEdge::Edge * bitEdge = faceEdge->prev()->pair();

			int tanIdx = -1;
			int bitIdx = -1;

			int i = 0;
			for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance(), i++)
			{
				if (eit.current() == tanEdge) tanIdx = i;
				if (eit.current() == bitEdge) bitIdx = i;
			}
			nvDebugCheck(tanIdx != -1);
			nvDebugCheck(bitIdx != -1);

			const int k = valence - 1;
		    const float omega = PI / k;
		    float ctan = cosf(tanIdx * omega);
		    float stan = sinf(tanIdx * omega);
		    float cbit = cosf(bitIdx * omega);
		    float sbit = sinf(bitIdx * omega);

		    if (equal(ctan, 0.0f)) ctan = 0;
		    if (equal(stan, 0.0f)) stan = 0;
		    if (equal(cbit, 0.0f)) cbit = 0;
		    if (equal(sbit, 0.0f)) sbit = 0;

            if (k == 1)
            {
                ctan = -1;
                cbit = 1;
                stan = 0.0001f;
                sbit = 0.0001f;
            }

			for (uint i = 0; i < vertexCount; i++)
			{
				stencil->tangentStencil(tangentIndices[v])[i] = r0[i] * stan + r1[i] * ctan;
				stencil->bitangentStencil(bitangentIndices[v])[i] = r0[i] * sbit + r1[i] * cbit;
			}

			if (v == 1 || v == 3) swap(stencil->tangentStencil(tangentIndices[v]), stencil->bitangentStencil(bitangentIndices[v]));

			stencil->tangentStencil(tangentIndices[v]) *= tangentFactors[v];
			stencil->bitangentStencil(bitangentIndices[v]) *= bitangentFactors[v];
		}
		else
		{
			// @@ Compute corner tangents as in theboundaries.

			Vector3 tan(zero);
			Vector3 bit(zero);

			// We do this to ensure consistent orientation or the tangents.
			const HalfEdge::Edge * start = it.current();
			if (v == 1 || v == 3) start = it.current()->prev()->pair();

			int i = 0;
			for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges(start)); !eit.isDone(); eit.advance(), i++)
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(edge != NULL);
				nvDebugCheck(edge->from()->pos() == vertex->pos());

				stencil->tangentStencil(tangentIndices[v], edge->to()) += alphaWeight(i, valence);

				float beta = betaWeight(i, valence);
				if (edge->pair()->prev()->prev()->prev() == edge->pair())
				{
					stencil->tangentStencil(tangentIndices[v], vertex) += beta / 3.0f;
					stencil->tangentStencil(tangentIndices[v], edge->pair()->from()) += beta / 3.0f;
					stencil->tangentStencil(tangentIndices[v], edge->pair()->prev()->from()) += beta / 3.0f;
				}
				else
				{
					stencil->tangentStencil(tangentIndices[v], edge->pair()->prev()->from()) += beta;
				}
			}

			start = it.current();
			if (v == 0 || v == 2) start = it.current()->prev()->pair();

			i = 0;
			for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges(start)); !eit.isDone(); eit.advance(), i++)
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(edge != NULL);
				nvDebugCheck(edge->from()->pos() == vertex->pos());

				stencil->bitangentStencil(bitangentIndices[v], edge->to()) += alphaWeight(i, valence);

				float beta = betaWeight(i, valence);
				if (edge->pair()->prev()->prev()->prev() == edge->pair())
				{
					stencil->bitangentStencil(bitangentIndices[v], vertex) += beta / 3.0f;
					stencil->bitangentStencil(bitangentIndices[v], edge->pair()->from()) += beta / 3.0f;
					stencil->bitangentStencil(bitangentIndices[v], edge->pair()->prev()->from()) += beta / 3.0f;
				}
				else
				{
					stencil->bitangentStencil(bitangentIndices[v], edge->pair()->prev()->from()) += beta;
				}
			}

			stencil->tangentStencil(tangentIndices[v]) *= tangentFactors[v];
			stencil->bitangentStencil(bitangentIndices[v]) *= bitangentFactors[v];
		}
	}

}

void AccMeshBuilder::computeEdgeTangentStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);

	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	nvDebugCheck(firstEdge != NULL);

	BezierAccStencil * stencil = patch->stencil();

	static const uint v10Index[4] = {1, 9, 10, 2};
	static const uint v20Index[4] = {2, 10, 9, 1};
	static const float vFactor[4] = {1, -1, -1, 1};

	static const uint u00Index[4] = {0, 3, 11, 8};
	static const uint u10Index[4] = {4, 7, 7, 4};
	static const uint u20Index[4] = {8, 11, 3, 0};
	static const float uFactor[4] = {1, 1, -1, -1};

	static const uint b10Index[4] = {1, 7, 14, 8};
	static const uint b11Index[4] = {5, 6, 10, 9};

	static const uint b20Index[4] = {2, 11, 13, 4};
	static const uint b21Index[4] = {6, 10, 9, 5};

	// Edge tangents
	uint v = 0;
	for(HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Edge * edge = it.current();
		nvDebugCheck(edge != NULL);

		const HalfEdge::Vertex * v0 = edge->from();
		nvDebugCheck(v0 != NULL);

		const float valence0 = pseudoValence(v0);
		const float c0 = cosf(2 * PI / valence0);

		const HalfEdge::Vertex * v1 = edge->to();
		nvDebugCheck(v1 != NULL);

		const float valence1 = pseudoValence(v1);
		const float c1 = cosf(2 * PI / valence1);

		// Special case for regular interior edges, to make sure tangents are exact.
		// This allows us to ignore boundary conditions with regular patches.
		/*if (!edge->isBoundary() && valence0 == 4 && valence1 == 4)
		{
			nvCheck(equal(c0, 0.0f));
			nvCheck(equal(c1, 0.0f));

			const HalfEdge::Edge * edge0 = edge;
			const HalfEdge::Edge * edge1 = edge->pair();
			nvDebugCheck(edge0 != NULL && edge1 != NULL);

			const HalfEdge::Face * face0 = edge0->face();
			const HalfEdge::Face * face1 = edge1->face();
			nvDebugCheck(face0 != NULL && face1 != NULL);

			::interiorControlPoint(face0, edge0, patch, v10Index[v], (v & 1), 0, 1.5f * vFactor[v]);
			::interiorControlPoint(face1, edge1->next(), patch, v10Index[v], (v & 1), 1, -1.5f * vFactor[v]);

			::interiorControlPoint(face0, edge0->next(), patch, v20Index[v], (v & 1), 0, 1.5f * vFactor[v]);
			::interiorControlPoint(face1, edge1, patch, v20Index[v], (v & 1), 1, -1.5f * vFactor[v]);
		}
		else*/
		{
			const uint vertexCount = patch->vertexCount();
			for (uint i = 0; i < vertexCount; i++)
			{
				float b10 = stencil->positionStencil(b10Index[v])[i];
				float b11 = stencil->positionStencil(b11Index[v])[i];
				float b20 = stencil->positionStencil(b20Index[v])[i];
				float b21 = stencil->positionStencil(b21Index[v])[i];

				float u10, u00, u20;
				if (v & 1) {
					u00 = stencil->bitangentStencil(u00Index[v])[i] * uFactor[v];
					u10 = stencil->bitangentStencil(u10Index[v])[i] * uFactor[v];
					u20 = stencil->bitangentStencil(u20Index[v])[i] * uFactor[v];
				}
				else {
					u00 = stencil->tangentStencil(u00Index[v])[i] * uFactor[v];
					u10 = stencil->tangentStencil(u10Index[v])[i] * uFactor[v];
					u20 = stencil->tangentStencil(u20Index[v])[i] * uFactor[v];
				}

				float v10 = (2 * c0 * u10 - 1 * c1 * u00) / 3 + 3 * (b11 - b10);
				float v20 = (1 * c0 * u20 - 2 * c1 * u10) / 3 + 3 * (b21 - b20);

				if (v & 1) {
					stencil->tangentStencil(v10Index[v])[i] = v10 * vFactor[v];
					stencil->tangentStencil(v20Index[v])[i] = v20 * vFactor[v];
				}
				else {
					stencil->bitangentStencil(v10Index[v])[i] = v10 * vFactor[v];
					stencil->bitangentStencil(v20Index[v])[i] = v20 * vFactor[v];
				}
			}
		}
	}
}

void AccMeshBuilder::computeInteriorTangentStencil(BezierAccPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	BezierAccStencil * stencil = patch->stencil();

	const uint vertexCount = patch->vertexCount();
	for (uint i = 0; i < vertexCount; i++)
	{
		stencil->tangentStencil(4)[i] = 3 * (stencil->positionStencil(2)[i] - stencil->positionStencil(1)[i]);
		stencil->tangentStencil(5)[i] = 3 * (stencil->positionStencil(6)[i] - stencil->positionStencil(5)[i]);
		stencil->tangentStencil(6)[i] = 3 * (stencil->positionStencil(10)[i] - stencil->positionStencil(9)[i]);
		stencil->tangentStencil(7)[i] = 3 * (stencil->positionStencil(14)[i] - stencil->positionStencil(13)[i]);

		stencil->bitangentStencil(4)[i] = 3 * (stencil->positionStencil(8)[i] - stencil->positionStencil(4)[i]);
		stencil->bitangentStencil(5)[i] = 3 * (stencil->positionStencil(9)[i] - stencil->positionStencil(5)[i]);
		stencil->bitangentStencil(6)[i] = 3 * (stencil->positionStencil(10)[i] - stencil->positionStencil(6)[i]);
		stencil->bitangentStencil(7)[i] = 3 * (stencil->positionStencil(11)[i] - stencil->positionStencil(7)[i]);
	}
}
/*
void AccMeshBuilder::computeSymmetricTangentStencil(BezierAccPatch * patch) const
{
	// Interior symmetric traversal tangents are easy:
	const uint vertexCount = patch->vertexCount();
	for (uint i = 0; i < vertexCount; i++)
	{
		// w10 = 3 * (b11 - b10)
		// w20 = 3 * (b21 - b20)
		
		patch->symmetricTangentStencil(1)[i] = 3 * (patch->positionStencil(5)[i] - patch->positionStencil(4)[i]);
		patch->symmetricTangentStencil(2)[i] = 3 * (patch->positionStencil(9)[i] - patch->positionStencil(8)[i]);
		patch->symmetricTangentStencil(5)[i] = 3 * (patch->positionStencil(7)[i] - patch->positionStencil(6)[i]);
		patch->symmetricTangentStencil(6)[i] = 3 * (patch->positionStencil(11)[i] - patch->positionStencil(10)[i]);

		patch->symmetricBitangentStencil(1)[i] = 3 * (patch->positionStencil(5)[i] - patch->positionStencil(1)[i]);
		patch->symmetricBitangentStencil(2)[i] = 3 * (patch->positionStencil(6)[i] - patch->positionStencil(2)[i]);
		patch->symmetricBitangentStencil(5)[i] = 3 * (patch->positionStencil(13)[i] - patch->positionStencil(9)[i]);
		patch->symmetricBitangentStencil(6)[i] = 3 * (patch->positionStencil(14)[i] - patch->positionStencil(10)[i]);
	}
	
	// w00 can be computed as:
	// w00 = v00 - c0 * u00
	
	// or as:
	// w00 = sum(alpha~i * pi + beta~i * qi)
	
	// no matter how, we have to make sure the result is the same for adjacent patches.

}
*/


//-------------------------------
//  GregoryAccPatch 
//-------------------------------

GregoryAccPatch * AccMeshBuilder::buildGregoryAccPatch(const FaceOneRing * face) const
{
	const uint edgeCount = face->edgeCount();
	nvDebugCheck(edgeCount == 3 || edgeCount == 4);

	GregoryAccPatch * patch = new GregoryAccPatch(face);

	computeCornerStencil(patch);

	computeEdgeStencil(patch);

	computeInteriorStencil(patch);

	patch->evaluateControlPoints();

	FaceTopology * faceTopology = patch->faceOneRing()->faceTopology();
	if (faceTopology->gregoryAccStencil == NULL)
	{
		faceTopology->gregoryAccStencil = patch->stencil();
	}
	else
	{
		// Make sure stencils are the same. If this fails this is probably due to self-overlapping stencils.
		nvCheck(*faceTopology->gregoryAccStencil == *patch->stencil());

		// @@ Delete new stencils, reference the same.
	}

	return patch;
}

void AccMeshBuilder::computeCornerStencil(GregoryAccPatch * patch) const
{
	//                                         /--- quads ---/   /--- tris ---/
	static const uint cornerIndices[7]      = { 8, 11, 19, 16,     6,  9, 12};

	const uint primitiveOffset = patch->isQuad() ? 0 : 4;

	const HalfEdge::Face * face = patch->face();
	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	GregoryAccStencil * stencil = patch->stencil();

	// Compute corner control points.
	uint v = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Vertex * vertex = it.current()->from();
		const uint valence = vertex->valence();
		const uint cornerIndex = cornerIndices[primitiveOffset+v];

		if (vertex->isBoundary())
		{
			nvDebugCheck(m_boundaryMode == BoundaryMode_Spline);

			// Compute limit vertex position.
			stencil->positionStencil(cornerIndex, vertex) = 2.0f / 3.0f;

			const HalfEdge::Edge * edge0 = vertex->edge();
			nvCheck(edge0->face() == NULL);
			nvCheck(edge0->to() != vertex);

			stencil->positionStencil(cornerIndex, edge0->to()) = 1.0f / 6.0f;

			const HalfEdge::Edge * edge1 = vertex->edge()->prev();
			nvCheck(edge1->face() == NULL);
			nvCheck(edge1->from() != vertex);

			stencil->positionStencil(cornerIndex, edge1->from()) = 1.0f / 6.0f;

			nvDebugCheck(stencil->positionStencil(cornerIndex).isNormalized());
		}
		else
		{
			stencil->positionStencil(cornerIndex, vertex) = float(3 * valence * valence);

			for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edge()); !eit.isDone(); eit.advance())
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(vertex->pos() == edge->from()->pos());

				stencil->positionStencil(cornerIndex, edge->to()) = 12.0f;

				if (FaceOneRing::isTriangle(edge->face()))
				{
					// Distribute weight to all vertices.
					stencil->positionStencil(cornerIndex, vertex) += 1.0f;
					stencil->positionStencil(cornerIndex, edge->to()) += 1.0f;
					stencil->positionStencil(cornerIndex, edge->next()->to()) += 1.0f;
				}
				else
				{
					stencil->positionStencil(cornerIndex, edge->next()->to()) = 3.0f;
				}
			}

			// Normalize stencil.
			stencil->positionStencil(cornerIndex).normalize();
		}
	}
}


void AccMeshBuilder::computeEdgeStencil(GregoryAccPatch * patch) const
{
	//                                         /--- quads ---/   /--- tris ---/
	static const uint cornerIndices[7]      = { 8, 11, 19, 16,     6,  9, 12};
	static const uint edge1Indices[7]       = { 9, 13, 18, 14,     7, 10, 13};
	static const uint edge2Indices[7]       = {12, 10, 15, 17,    14,  8, 11};

	const uint primitiveOffset = patch->isQuad() ? 0 : 4;

	const float tangentScales[14] = {
		0.0f, 0.0f, 0.0f, 0.667791f, 1.0f, 1.11268f, 1.1284f, 1.10289f, 1.06062f, 1.01262f, 0.963949f, 0.916926f, 0.872541f, 0.831134f
	};

	const HalfEdge::Face * face = patch->face();
	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	GregoryAccStencil * stencil = patch->stencil();

	// Compute corner / edge control points.
	uint v = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Vertex * vertex = it.current()->from();
		const uint valence = vertex->valence();
		const uint cornerIndex = cornerIndices[primitiveOffset+v];

		int i1 = 0, i2 = 0, j = 0;
		for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edge()); !eit.isDone(); eit.advance(), j++)
		{
			const HalfEdge::Edge * edge = eit.current();

			// find index of "our" edge for edge control points
			if (edge == it.current()) {
				i1 = j;
			}
			if (edge == it.current()->prev()->pair()) {
				i2 = j;
			}
		}

		if (vertex->isBoundary())
		{
			const uint vertexCount = patch->vertexCount();
			StencilMask	r0(vertexCount);
			StencilMask	r1(vertexCount);

			computeBoundaryTangentStencils(patch, vertex, r0, r1);

			const int k = valence - 1;
			const float omega = PI / k;

			const uint edgeIndex1 = edge1Indices[primitiveOffset + v];
			const uint edgeIndex2 = edge2Indices[primitiveOffset + v];

			if (it.current()->isBoundary())
			{
				nvDebugCheck(m_boundaryMode == BoundaryMode_Spline);

				nvDebugCheck(it.current()->from() == vertex);

			    stencil->positionStencil(edgeIndex1, vertex) = 2.0f / 3.0f;
			    stencil->positionStencil(edgeIndex1, it.current()->to()) = 1.0f / 3.0f;

                nvDebugCheck(stencil->positionStencil(edgeIndex1).isNormalized());

                if (valence == 2)
                {
			        for (uint i = 0; i < vertexCount; i++)
			        {
				        stencil->positionStencil(edgeIndex1)[i] += r0[i] * 0.0001f;
			        }                    
                }
			}
			else
			{
				stencil->positionStencil(edgeIndex1) = stencil->positionStencil(cornerIndex);

				// Compute index of it.current() around vertex.
				int idx = 0;
				for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance(), idx++)
				{
					if (eit.current() == it.current()) break;
				}
				nvDebugCheck(idx != valence);

				const float c = cosf(idx * omega);
				const float s = sinf(idx * omega);

				for (uint i = 0; i < vertexCount; i++)
				{
					stencil->positionStencil(edgeIndex1)[i] += (r0[i] * s + r1[i] * c) / 3.0f;
				}
			}

			if (it.current()->prev()->isBoundary())
			{
				nvDebugCheck(m_boundaryMode == BoundaryMode_Spline);

				nvDebugCheck(it.current()->prev()->pair()->from() == vertex);

				stencil->positionStencil(edgeIndex2, vertex) = 2.0f / 3.0f;
				stencil->positionStencil(edgeIndex2, it.current()->prev()->pair()->to()) = 1.0f / 3.0f;

				nvDebugCheck(stencil->positionStencil(edgeIndex2).isNormalized());

                if (valence == 2)
                {
			        for (uint i = 0; i < vertexCount; i++)
			        {
				        stencil->positionStencil(edgeIndex2)[i] += r0[i] * 0.0001f;
			        }                    
                }
			}
			else
			{
				
				stencil->positionStencil(edgeIndex2) = stencil->positionStencil(cornerIndex);

				// Compute index of it.current() around vertex.
				int idx = 0;
				for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance(), idx++)
				{
					if (eit.current() == it.current()->prev()->pair()) break;
				}
				nvDebugCheck(idx != valence);

				const float c = cosf(idx * omega);
				const float s = sinf(idx * omega);

				for (uint i = 0; i < vertexCount; i++)
				{
					stencil->positionStencil(edgeIndex2)[i] += (r0[i] * s + r1[i] * c) / 3;
				}
			}
		}
		else
		{
			const float costerm = cosf(PI / valence);
			const float sqrtterm = sqrtf(4.0f + costerm*costerm);

			//const float tangentScale = 1.0f;
			const float tangentScale = tangentScales[min(valence, 13U)];

			const float alpha = (1.0f +  costerm / sqrtterm) / (3.0f * valence) * tangentScale;
			const float beta  = 1.0f / (3.0f * valence * sqrtterm) * tangentScale;


			const uint edgeIndex1 = edge1Indices[primitiveOffset + v];
			const uint edgeIndex2 = edge2Indices[primitiveOffset + v];

			stencil->positionStencil(edgeIndex1) = stencil->positionStencil(cornerIndex);
			stencil->positionStencil(edgeIndex2) = stencil->positionStencil(cornerIndex);

			int j = 0;
			for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance(), j++)
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(vertex->pos() == edge->from()->pos());

				const float costerm1_a = cosf(PI * 2 * (j-i1) / valence);
				const float costerm1_b = cosf(PI * (2 * (j-i1)-1) / valence); // -1 instead of +1 b/c of edge->next()->to()

				const float costerm2_a = cosf(PI * 2 * (j-i2) / valence);
				const float costerm2_b = cosf(PI * (2 * (j-i2)-1) / valence); // -1 instead of +1 b/c of edge->next()->to()


				stencil->positionStencil(edgeIndex1, edge->to()) += alpha * costerm1_a;
				stencil->positionStencil(edgeIndex2, edge->to()) += alpha * costerm2_a;

				if (FaceOneRing::isTriangle(edge->face()))
				{
					// Distribute weight to all vertices.
					stencil->positionStencil(edgeIndex1, vertex) += beta * costerm1_b / 3.0f;				// @@ This probably does not provide watertight results!! (1/3 + 1/3 + 1/3 != 1)
					stencil->positionStencil(edgeIndex1, edge->to()) += beta * costerm1_b / 3.0f;
					stencil->positionStencil(edgeIndex1, edge->next()->to()) += beta * costerm1_b / 3.0f;

					stencil->positionStencil(edgeIndex2, vertex) += beta * costerm2_b / 3.0f;				// @@ This probably does not provide watertight results!! (1/3 + 1/3 + 1/3 != 1)
					stencil->positionStencil(edgeIndex2, edge->to()) += beta * costerm2_b / 3.0f;
					stencil->positionStencil(edgeIndex2, edge->next()->to()) += beta * costerm2_b / 3.0f;
				}
				else
				{
					stencil->positionStencil(edgeIndex1, edge->next()->to()) += beta * costerm1_b;

					stencil->positionStencil(edgeIndex2, edge->next()->to()) += beta * costerm2_b;
				}
			}
		}
	}
}


void AccMeshBuilder::computeInteriorStencil(GregoryAccPatch * patch) const
{
	//                                      /--- quads ---/   /--- tris ---/
	static const uint corner1Indices[7]   = { 8, 11, 19, 16,     6,  9, 12};
	static const uint corner2Indices[7]   = {11, 19, 16,  8,     9, 12,  6};
	static const uint edge1Indices[7]     = { 9, 13, 18, 14,     7, 10, 13};
	static const uint edge2Indices[7]     = {10, 15, 17, 12,     8, 11, 14};
	static const uint interior1Indices[7] = { 1,  3,  6,  4,     1,  3,  5};
	static const uint interior2Indices[7] = { 2,  7,  5,  0,     2,  4,  0};

	const uint primitiveOffset = patch->isQuad() ? 0 : 4;
	const uint vertexCount = patch->isQuad() ? 4 : 3;

	const HalfEdge::Face * face = patch->face();
	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	GregoryAccStencil * stencil = patch->stencil();

	// interior control points
	uint v = 0;
	for(HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Edge * edge = it.current();

		if (edge->isBoundary())
		{
			const float valence1 = pseudoValence(edge->from());
			const float valence2 = pseudoValence(edge->to());

			float weights1[4];
			float weights2[4];
			if (patch->isQuad())
			{
				weights1[0] = 3 * valence1;
				weights1[1] = 6;
				weights1[2] = 3;
				weights1[3] = 6;

				weights2[0] = 6;
				weights2[1] = 3 * valence2;
				weights2[2] = 6;
				weights2[3] = 3;
			}
			else
			{
				nvDebugCheck(patch->isTriangle());
				weights1[0] = 3 * valence1 + 1;
				weights1[1] = 7;
				weights1[2] = 7;

				weights2[0] = 7;
				weights2[1] = 3 * valence2 + 1;
				weights2[2] = 7;
			}

			uint idx1 = interior1Indices[primitiveOffset+v];
			uint idx2 = interior2Indices[primitiveOffset+v];

			uint i = 0;
			for (HalfEdge::Face::ConstEdgeIterator it(face->edges(edge)); !it.isDone(); it.advance(), i++)
			{
				const HalfEdge::Vertex * vertex = it.current()->from();
				stencil->positionStencil(idx1, vertex) += weights1[i];
				stencil->positionStencil(idx2, vertex) += weights2[i];
			}

			stencil->positionStencil(idx1).normalize();
			stencil->positionStencil(idx2).normalize();
		}
		else
		{
			const HalfEdge::Vertex * e0 = edge->from();
			const float costerm0 = cosf(2.0f * PI / pseudoValence(e0));

			const HalfEdge::Vertex * f0 = edge->to();
			const float costerm1 = cosf(2.0f * PI / pseudoValence(f0));

			//  p0 +------+ q0
			//     |      |
			//  f0 +======+ e0 <=== current edge
			//     |      |
			//  p1 +------+ q1

			const HalfEdge::Vertex * q0 = edge->next()->to();
			const HalfEdge::Vertex * p0 = edge->prev()->from();

			const HalfEdge::Vertex * p1 = edge->pair()->next()->to();
			const HalfEdge::Vertex * q1 = edge->pair()->prev()->from();


			StencilMask	x(patch->vertexCount());
			StencilMask	y(patch->vertexCount());

			for (uint i = 0; i < patch->vertexCount(); i++)
			{
				x[i] =
					(costerm1 * stencil->positionStencil(corner1Indices[primitiveOffset+v])[i] -
					(2*costerm0 + costerm1) * stencil->positionStencil(edge1Indices[primitiveOffset+v])[i] +
					2*costerm0 * stencil->positionStencil(edge2Indices[primitiveOffset+v])[i]) / 3.0f;
			}

			// y = (2*( midedgeA1 - midedgeB1) + 4*(centroidA - centroidB))/18.0f;
			y[patch->vertexIndex(p0)] = 1;
			y[patch->vertexIndex(p1)] = -1;

			// Add centroidA
			if (patch->isTriangle())
			{
				y[patch->vertexIndex(p0)] += 4.0f / 3.0f;
				y[patch->vertexIndex(e0)] += 4.0f / 3.0f;
				y[patch->vertexIndex(f0)] += 4.0f / 3.0f;
			}
			else
			{
				y[patch->vertexIndex(p0)] += 1;
				y[patch->vertexIndex(q0)] += 1;
				y[patch->vertexIndex(e0)] += 1;
				y[patch->vertexIndex(f0)] += 1;
			}

			// Sub centroidB
			if (FaceOneRing::isTriangle(edge->pair()->face()))
			{
				y[patch->vertexIndex(p1)] -= 4.0f / 3.0f;
				y[patch->vertexIndex(e0)] -= 4.0f / 3.0f;
				y[patch->vertexIndex(f0)] -= 4.0f / 3.0f;

			}
			else
			{
				y[patch->vertexIndex(p1)] -= 1;
				y[patch->vertexIndex(q1)] -= 1;
				y[patch->vertexIndex(e0)] -= 1;
				y[patch->vertexIndex(f0)] -= 1;
			}

			y /= 18.0f;

			if (patch->isTriangle())
			{
				x *= 3.0f / 4.0f;
				y *= 3.0f / 4.0f;
			}

			// This change makes the triangle boundaries smoother, but distorts the quads next to them.
			/*if (patch->isTriangle() || FaceOneRing::isTriangle(edge->pair()->face()))
			{
				y *= 4.0f / 3.0f;
			}*/

			stencil->positionStencil(interior1Indices[primitiveOffset+v]) = stencil->positionStencil(edge1Indices[primitiveOffset+v]);
			stencil->positionStencil(interior1Indices[primitiveOffset+v]) += x;
			stencil->positionStencil(interior1Indices[primitiveOffset+v]) += y;



			for (uint i = 0; i < patch->vertexCount(); i++)
			{
				x[i] =
					(costerm0 * stencil->positionStencil(corner2Indices[primitiveOffset+v])[i] -
					(2*costerm1 + costerm0) * stencil->positionStencil(edge2Indices[primitiveOffset+v])[i] + 
					2*costerm1 * stencil->positionStencil(edge1Indices[primitiveOffset+v])[i]) / 3.0f;
			}

			// y = (2*( midedgeA2 - midedgeB2) + 4*(centroidA - centroidB))/18.0f;
			y = 0.0f;

			// (2*( midedgeA2 - midedgeB2)
			y[patch->vertexIndex(q0)] = 1;
			y[patch->vertexIndex(q1)] = -1;

			// Add centroidA
			if (patch->isTriangle())
			{
				y[patch->vertexIndex(p0)] += 4.0f / 3.0f;
				y[patch->vertexIndex(e0)] += 4.0f / 3.0f;
				y[patch->vertexIndex(f0)] += 4.0f / 3.0f;
			}
			else
			{
				y[patch->vertexIndex(p0)] += 1;
				y[patch->vertexIndex(q0)] += 1;
				y[patch->vertexIndex(e0)] += 1;
				y[patch->vertexIndex(f0)] += 1;
			}

			// Sub centroidB
			if (FaceOneRing::isTriangle(edge->pair()->face()))
			{
				y[patch->vertexIndex(p1)] -= 4.0f / 3.0f;
				y[patch->vertexIndex(e0)] -= 4.0f / 3.0f;
				y[patch->vertexIndex(f0)] -= 4.0f / 3.0f;

			}
			else
			{
				y[patch->vertexIndex(p1)] -= 1;
				y[patch->vertexIndex(q1)] -= 1;
				y[patch->vertexIndex(e0)] -= 1;
				y[patch->vertexIndex(f0)] -= 1;
			}

			y /= 18.0f;

			if (patch->isTriangle())
			{
				x *= 3.0f / 4.0f;
				y *= 3.0f / 4.0f;
			}

			// This change makes the triangle boundaries smoother, but distorts the quads next to them.
			/*if (patch->isTriangle() || FaceOneRing::isTriangle(edge->pair()->face()))
			{
				y *= 4.0f / 3.0f;
			}*/

			stencil->positionStencil(interior2Indices[primitiveOffset+v]) = stencil->positionStencil(edge2Indices[primitiveOffset+v]);
			stencil->positionStencil(interior2Indices[primitiveOffset+v]) += x;
			stencil->positionStencil(interior2Indices[primitiveOffset+v]) += y;
		}
	}
}




//-------------------------------
//  PolarPatch 
//-------------------------------

PolarPatch * AccMeshBuilder::buildPolarPatch(const FaceOneRing * face) const
{
	nvDebugCheck(face->edgeCount() == 4);

	// ...

	return NULL;
}


//-------------------------------
//  PmQuadPatch 
//-------------------------------

PmQuadAccPatch * AccMeshBuilder::buildPmQuadAccPatch(const FaceOneRing * face) const
{
	const uint edgeCount = face->edgeCount();
	nvDebugCheck(edgeCount == 3 || edgeCount == 4);

	PmQuadAccPatch * patch = new PmQuadAccPatch(face);

	computeCoefficients(patch);

	return patch;
}


void AccMeshBuilder::computeCoefficients(PmQuadAccPatch * patch) const
{
	float sigma[7]={0.41f, 0.5f, 0.55f, 0.5797f, 0.5985f, 0.6111f, 0.62f};
    const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);
	Vector3 vertices[4], cornerPos[4];
	uint valence[4], index[8];
	const HalfEdge::Edge * edges[4];
	uint v = 0;

    const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), ++v)
	{
		// vertices in ccw order (use right-hand coordinate system)
		vertices[v] = it.current()->from()->pos();
		valence[v] = it.current()->from()->valence();
		edges[v] = it.current();
	}

	uint sides, maxValence;
    sides=4;
	maxValence= max(max(valence[0], valence[1]), max(valence[2],valence[3]));
	nvCheck(maxValence - 3 < 7);

	Array<Vector3> tangentPos; tangentPos.resize(sides * maxValence);
	Array<Vector3> edgePos; edgePos.resize(sides * maxValence);
	Array<Vector3> directNbr; directNbr.resize(maxValence);
	Array<Vector3> diagNbr; diagNbr.resize(maxValence);

	// vertex-based computations
	for (v = 0; v < sides; ++v)
	{
	    // compute corner points
		cornerPos[v]=Vector3(zero);
		const HalfEdge::Vertex * vertex = edges[v]->from();
		uint i=0;
		for(HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance())
		{
			const HalfEdge::Edge * edge = eit.current();
			nvDebugCheck(vertices[v] == edge->from()->pos());
		    // directNbr --> diagNbr in cw order
			directNbr[i]=edge->to()->pos();
			diagNbr[i]=edge->pair()->prev()->from()->pos();
			if (diagNbr[i]==edge->pair()->next()->to()->pos())
               diagNbr[i]= (edge->to()->pos() + edge->pair()->next()->to()->pos())*0.5;

			cornerPos[v] += float(valence[v])*vertices[v] + 4.0f*directNbr[i] + diagNbr[i];

			i++;
		}
		cornerPos[v] = cornerPos[v]/float(valence[v]*(valence[v]+5));

		// compute edge points
		for (uint i=0; i<valence[v]; ++i)
		{
			int next=(i+1)%valence[v];
		    int pre=(i+valence[v]-1)%valence[v];
		    edgePos[v*maxValence+i] = (8*vertices[v] + 2*(directNbr[pre] + directNbr[next]) 
			         + 4*directNbr[i] + (diagNbr[i] + diagNbr[pre]))/18;
			if (directNbr[i] == vertices[(v+1)%sides]) 
				index[2*v]=i;
			if (directNbr[i] == vertices[(v+sides-1)%sides]) 
				index[2*v+1]=i;
		}

 
		// compute tangent points
		for (uint i=0; i<valence[v]; ++i)
		{
			tangentPos[v*maxValence+i]=Vector3(zero);
			for (uint j=0; j<valence[v]; j++) {
				tangentPos[v*maxValence+i]   += cosf((2*PI/valence[v])*((int)i-(int)j))*edgePos[v*maxValence+j];
			}
			tangentPos[v*maxValence+i]=cornerPos[v]+(1.0f/(sigma[valence[v]-3]*valence[v]))*tangentPos[v*maxValence+i];
		}

	}


	Vector3 facePos[4], b211[4], b121[4], b112[4];
	//compute face points
	for (uint v=0; v<sides; ++v) 
	{
		uint next=(v+1)%sides;
		uint prev=(v+sides-1)%sides;
		uint op=(v+2)%sides;
		facePos[v]=(4*vertices[v] + 2*(vertices[next] + vertices[prev]) + vertices[op])/9;

	}
    //facet-based computations: (2)compute b211 and b121
	Vector3 center = Vector3(zero);
	for (uint v=0; v<sides; ++v) 
	{
		uint next=(v+1)%sides;
		uint currentIndex = index[2*v];
		uint preIndex = index[2*next+1];
		
		Vector3 U0=tangentPos[v*maxValence+currentIndex]-cornerPos[v];
		Vector3 U1=tangentPos[next*maxValence+preIndex]-tangentPos[v*maxValence+currentIndex];  
		Vector3 U2=cornerPos[next]-tangentPos[next*maxValence+preIndex];  
		float mu=1-cosf(2*PI/sides);
		Vector3 pert=3*(facePos[v]- edgePos[v*maxValence+currentIndex])/(4*mu*(sinf(2*PI/valence[v])+sinf(2*PI/valence[next])));
		b211[v] = (cornerPos[v]+3*tangentPos[v*maxValence+currentIndex])/4 + (1+cosf(2*PI/valence[v]))/4*U1
			        +(1-cosf(2*PI/valence[next]))/8*U0+ pert;	
		
		float s1=(1+cosf(2*PI/valence[v]))/(4*mu);
		float s2=(2*mu-(1+cosf(2*PI/valence[next])))/(8*mu);

		b211[v] = (cornerPos[v]+3*tangentPos[v*maxValence+currentIndex])/4 + s1*U1
			        +s2*U0+ pert;
		// reverse orientation
		pert=3*(facePos[next]- edgePos[next*maxValence+preIndex])/(4*mu*(sinf(2*PI/valence[v])+sinf(2*PI/valence[next])));

		s1=(1+cosf(2*PI/valence[next]))/(4*mu);
		s2=(2*mu-(1+cosf(2*PI/valence[v])))/(8*mu);
		b121[v]= (cornerPos[next]+3*tangentPos[next*maxValence+preIndex])/4 - s1*U1
			       -s2*U2 + pert;

		// center: bi-cubic patch evaluated at (0.5, 0.5)
  		   
		center += cornerPos[v] +9*facePos[v] + 3*(edgePos[v*maxValence+currentIndex]+edgePos[next*maxValence+preIndex]);
	}
	center = center/float(sides*16);

	//facet-based computations: (3)compute b112
	// assign 24 coefficients per patch, 6 per sector as follows
	//         5
	//       3  4
	//    0--1--2--

	b112[0]=(3.0f*b211[0] - 3.0f*b121[1] + b211[1] - b121[2] - b211[2] + b121[3] - 3.0f*b211[3] + 3.0f*b121[0])/16.0f + center;
	b112[1]=(3.0f*b211[1] - 3.0f*b121[2] - 3.0f*b211[0] + 3.0f*b121[1] + b211[2] - b121[3] - b211[3] + b121[0])/16.0f + center;
	b112[2]=(-3.0f*b211[1] + 3.0f*b121[2] - b211[0] + b121[1] + 3.0f*b211[2] - 3.0f*b121[3] + b211[3] - b121[0])/16.0f + center;
	b112[3]=(-b211[1] + b121[2] + b211[0] - b121[1] - 3.0f*b211[2] + 3.0f*b121[3] + 3.0f*b211[3] - 3.0f*b121[0])/16.0f + center; 
    for (v=0; v<sides; ++v) {
		uint next=(v+1)%sides;
		uint currentIndex = index[2*v];
		uint preIndex = index[2*next+1];
		patch->m_positionArray[6*v]=cornerPos[v];
		patch->m_positionArray[6*v+1]=tangentPos[v*maxValence+currentIndex];
		patch->m_positionArray[6*v+2]=tangentPos[next*maxValence+preIndex];
		patch->m_positionArray[6*v+3]=b211[v];
		patch->m_positionArray[6*v+4]=b121[v];
		patch->m_positionArray[6*v+5]=b112[v];
    }
}


//-------------------------------
//  PmTriangleAccPatch 
//-------------------------------

PmTriangleAccPatch * AccMeshBuilder::buildPmTriangleAccPatch(const FaceOneRing * face) const
{
	const uint edgeCount = face->edgeCount();
	nvDebugCheck(edgeCount == 3);

	PmTriangleAccPatch * patch = new PmTriangleAccPatch(face);

	convertToP3patch(patch);

	return patch;
}

void AccMeshBuilder::convertToP3patch(PmTriangleAccPatch * patch) const
{
	float sigma[7]={0.41f, 0.5f, 0.55f, 0.5797f, 0.5985f, 0.6111f, 0.62f};
	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 3);

	Vector3 vertices[3], cornerPos[3];
	uint valence[3], index[6];
	const HalfEdge::Edge * edges[6];
	uint v = 0;
	
    const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), ++v)
	{
		// vertices in ccw order (use right-hand coordinate system)
		vertices[v] = it.current()->from()->pos();
		valence[v] = it.current()->from()->valence();
		edges[v] = it.current();
	}

	uint sides, maxValence;
    sides=3;

	maxValence= max(max(valence[0], valence[1]), valence[2]);
	nvCheck(maxValence - 3 < 7);

	Array<Vector3> tangentPos; tangentPos.resize(3 * maxValence);
	Array<Vector3> edgePos; edgePos.resize(3 * maxValence);
	Array<Vector3> directNbr; directNbr.resize(maxValence);
	Array<Vector3> diagNbr; diagNbr.resize(maxValence);

	// vertex-based computations
	for (v = 0; v < sides; ++v)
	{
	    // compute corner points
		cornerPos[v]=Vector3(zero);
		const HalfEdge::Vertex * vertex = edges[v]->from();
		uint i=0;
		for(HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance())
		{
			const HalfEdge::Edge * edge = eit.current();
			nvDebugCheck(vertices[v] == edge->from()->pos());
		    // directNbr --> diagNbr in cw order
			directNbr[i]=edge->to()->pos();
			diagNbr[i]=edge->pair()->prev()->from()->pos();
			if (diagNbr[i]==edge->pair()->next()->to()->pos())
               diagNbr[i]= (edge->to()->pos() + edge->pair()->next()->to()->pos())*0.5;

			cornerPos[v] += float(valence[v])*vertices[v] + 4*directNbr[i] + diagNbr[i];

			i++;
		}
		cornerPos[v] = cornerPos[v]/float(valence[v]*(valence[v]+5));

		// compute edge points
		for (uint i=0; i<valence[v]; ++i)
		{
			int next=(i+1)%valence[v];
		    int pre=(i+valence[v]-1)%valence[v];
		    edgePos[v*maxValence+i] = (8*vertices[v] + 2*(directNbr[pre] + directNbr[next]) 
			         + 4*directNbr[i] + (diagNbr[i] + diagNbr[pre]))/18;
			if (directNbr[i] == vertices[(v+1)%sides]) 
				index[2*v]=i;
			if (directNbr[i] == vertices[(v+sides-1)%sides]) 
				index[2*v+1]=i;
		}

 
		// compute tangent points
		for (uint i=0; i<valence[v]; ++i)
		{
			tangentPos[v*maxValence+i]=Vector3(zero);
			for (uint j=0; j<valence[v]; j++) {
				tangentPos[v*maxValence+i]   += cosf((2*PI/valence[v])*((int)i-(int)j))*edgePos[v*maxValence+j];
			}
			tangentPos[v*maxValence+i]=cornerPos[v]+(1.0f/(sigma[valence[v]-3]*valence[v]))*tangentPos[v*maxValence+i];
		}

	}


	Vector3 facePos[3], b211[3], b121[3], b112[3];
	//compute face points
	for (uint v=0; v<sides; ++v) 
	{
		uint next=(v+1)%sides;
		uint prev=(v+sides-1)%sides;

        facePos[v]=(4*vertices[v] + 2*(vertices[next] + vertices[prev]) + (vertices[next]+vertices[prev])/2)/9;
	}
    //facet-based computations: (2)compute b211 and b121
	Vector3 center = Vector3(zero);
	const float w = 2;
	for (uint v=0; v<sides; ++v) 
	{
		uint next=(v+1)%sides;
		uint currentIndex = index[2*v];
		uint preIndex = index[2*next+1];
		
		Vector3 U0=tangentPos[v*maxValence+currentIndex]-cornerPos[v];
		Vector3 U1=tangentPos[next*maxValence+preIndex]-tangentPos[v*maxValence+currentIndex];  
		Vector3 U2=cornerPos[next]-tangentPos[next*maxValence+preIndex];  
		float mu=1-cosf(2*PI/sides);
		Vector3 pert=3*(facePos[v]- edgePos[v*maxValence+currentIndex])/(4*mu*(sinf(2*PI/valence[v])+sinf(2*PI/valence[next])));
		b211[v] = (cornerPos[v]+3*tangentPos[v*maxValence+currentIndex])/4 + (1+cosf(2*PI/valence[v]))/4*U1
			        +(1-cosf(2*PI/valence[next]))/8*U0+ pert;	
		
		float s1=(1+cosf(2*PI/valence[v]))/(4*mu);
		float s2=(2*mu-(1+cosf(2*PI/valence[next])))/(8*mu);
		b211[v] = (cornerPos[v]+3*tangentPos[v*maxValence+currentIndex])/4 + s1*U1
			        +s2*U0+ pert;
		// reverse orientation
		pert=3*(facePos[next]- edgePos[next*maxValence+preIndex])/(4.0f*mu*(sinf(2*PI/valence[v])+sinf(2*PI/valence[next])));

		s1=(1.0f+cosf(2.0f*PI/valence[next]))/(4.0f*mu);
		s2=(2.0f*mu-(1.0f+cosf(2.0f*PI/valence[v])))/(8.0f*mu);
		b121[v]= (cornerPos[next]+3.0f*tangentPos[next*maxValence+preIndex])/4.0f - s1*U1
			       -s2*U2 + pert;

		// center: bi-cubic patch evaluated at (0.5, 0.5)
		center += w*cornerPos[v] +9.0f*facePos[v] + 3.0f*(edgePos[v*maxValence+currentIndex]+edgePos[next*maxValence+preIndex]);
	}
	center = center/(sides*(15+w));

	//facet-based computations: (3)compute b112
	// assign 24 coefficients per patch, 6 per sector as follows
	//         5
	//       3  4
	//    0--1--2--

		// triangle facet
	for (v=0; v<sides; ++v) {
		uint next=(v+1)%sides;
		uint prev=(v+sides-1)%sides;
		uint currentIndex = index[2*v];
        uint mIndex = index[2*prev];
		uint pmIndex = index[2*prev+1];
        uint nmIndex = index[2*next+1];
		b112[v] =
			  (3.0f/2.0f ) * center
			- (1.0f/12.0f) * cornerPos[prev]
			- (1.0f/24.0f) * ( tangentPos[prev*maxValence+mIndex] + tangentPos[prev*maxValence+pmIndex] )
			- (1.0f/6.0f ) * ( b211[prev] + b121[next] );

		patch->m_positionArray[6*v]=cornerPos[v];
		patch->m_positionArray[6*v+1]=tangentPos[v*maxValence+currentIndex];
		patch->m_positionArray[6*v+2]=tangentPos[next*maxValence+nmIndex];
		patch->m_positionArray[6*v+3]=b211[v];
		patch->m_positionArray[6*v+4]=b121[v];
		patch->m_positionArray[6*v+5]=b112[v];
	}
	patch->m_positionArray[18]=center;
    
}


//-------------------------------
//  PmRegularAccPatch
//-------------------------------

PmRegularAccPatch * AccMeshBuilder::buildPmRegularAccPatch(const FaceOneRing * face) const
{
	const uint edgeCount = face->edgeCount();
	nvDebugCheck(edgeCount == 4);

	PmRegularAccPatch * patch = new PmRegularAccPatch(face);

	convertToBicubic(patch);

	return patch;
}

void AccMeshBuilder::convertToBicubic(PmRegularAccPatch * patch) const
{
	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);
	nvDebugCheck(FaceOneRing::isRegular(face));

	Vector3 vertices[4], cornerPos[4], edgePos[16], tangentPos[8], directNbr[4], diagNbr[4];
	uint index[8];
	const HalfEdge::Edge * edges[4];
	uint v=0;
	
    const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), ++v)
	{
		// vertices in ccw order (use right-hand coordinate system)
		vertices[v] = it.current()->from()->pos();
		nvDebugCheck(it.current()->from()->valence() == 4);
		edges[v] = it.current();
	}

	// vertex-based computations
	for (v = 0; v < 4; ++v)
	{
	    // compute corner points
		cornerPos[v]=Vector3(zero);
		const HalfEdge::Vertex * vertex = edges[v]->from();
		uint i=0;
		for(HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance())
		{
			const HalfEdge::Edge * edge = eit.current();
			nvDebugCheck(vertices[v] == edge->from()->pos());
			// directNbr --> diagNbr in cw order
			directNbr[i]=edge->to()->pos();
			diagNbr[i]=edge->pair()->prev()->from()->pos();
            if (diagNbr[i]==edge->pair()->next()->to()->pos())
            {
                diagNbr[i]= (edge->to()->pos() + edge->pair()->next()->to()->pos())*0.5;
            }

			cornerPos[v] += 4*vertices[v] + 4*directNbr[i] + diagNbr[i];
			i++;
		}
		cornerPos[v] = cornerPos[v]/36;
        nvDebugCheck(i == 4);
		// compute edge points
		for (i=0; i<4; ++i)
		{
			int next=(i+1)%4;
		    int pre=(i+3)%4;
		    edgePos[v*4+i] = (8*vertices[v] + 2*(directNbr[pre] + directNbr[next]) 
			         + 4*directNbr[i] + (diagNbr[i] + diagNbr[pre]))/18;
			if (directNbr[i] == vertices[(v+1)%4]) 
				index[2*v]=i;               //always 3, 2
			if (directNbr[i] == vertices[(v+3)%4]) 
				index[2*v+1]=i;
		}
	}

	Vector3 facePos[4];
	//compute face points
	for (v=0; v<4; ++v) 
	{
		uint next=(v+1)%4;
		uint prev=(v+3)%4;
		uint op=(v+2)%4;
		uint currentIndex = index[2*v];
		uint preIndex = index[2*v+1];
		facePos[v]=(4*vertices[v] + 2*(vertices[next] + vertices[prev]) + vertices[op])/9;
		tangentPos[2*v]=edgePos[4*v+currentIndex];
		tangentPos[2*v+1]=edgePos[4*v+preIndex];
	}
   
	//    12-13-14-15
	//    8--9--10-11   
	//    4--5--6--7 
	//    0--1--2--3
	patch->m_positionArray[0]=cornerPos[0];
	patch->m_positionArray[3]=cornerPos[1];
	patch->m_positionArray[15]=cornerPos[2];
	patch->m_positionArray[12]=cornerPos[3];
	patch->m_positionArray[1]=tangentPos[0];
	patch->m_positionArray[4]=tangentPos[1];
	patch->m_positionArray[7]=tangentPos[2];
	patch->m_positionArray[2]=tangentPos[3];
	patch->m_positionArray[14]=tangentPos[4];
	patch->m_positionArray[11]=tangentPos[5];
	patch->m_positionArray[8]=tangentPos[6];
	patch->m_positionArray[13]=tangentPos[7];
	patch->m_positionArray[5]=facePos[0];
	patch->m_positionArray[6]=facePos[1];
	patch->m_positionArray[10]=facePos[2];
	patch->m_positionArray[9]=facePos[3];
}


//-------------------------------
//  TexCoordPatch
//-------------------------------

TexCoordPatch * AccMeshBuilder::buildTexCoordPatch(const FaceOneRing * face) const
{
	TexCoordPatch * patch = new TexCoordPatch(face);

	computeTexcoords (patch);

	return patch;
}

void AccMeshBuilder::computeTexcoords(TexCoordPatch * patch) const
{
	nvDebugCheck(patch != NULL);

	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);

	const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
	nvDebugCheck(firstEdge != NULL);

	// texture coordinates
	uint v  = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Edge * edge = it.current();
		const HalfEdge::Vertex * vertex = edge->from();

		// Interior texcoord.
		patch->texCoord(v, 0) = vertex->tex();

		// Corner texcoord.
		uint minId = vertex->id();
		patch->texCoord(v, 3) = patch->texCoord(v, 0);

		for (HalfEdge::Vertex::ConstVertexIterator cit(vertex->colocals()); !cit.isDone(); cit.advance())
		{
			const HalfEdge::Vertex * colocal = cit.current();
			if (colocal->id() < minId)
			{
				minId = colocal->id();
				patch->texCoord(v, 3) = colocal->tex();
			}
		}

		// Boundary texcoords.
		const HalfEdge::Edge * pair = edge->pair();
		nvDebugCheck(pair != NULL);
		nvDebugCheck(edge->id() != pair->id());

		if (edge->id() < pair->id())
		{
			patch->texCoord(v, 1) = patch->texCoord(v, 0);
		}
		else
		{
			patch->texCoord(v, 1) = pair->to()->tex();
		}

		edge = edge->prev();
		pair = edge->pair();

		nvDebugCheck(pair != NULL);
		nvDebugCheck(edge->id() != pair->id());

		if (edge->id() < pair->id())
		{
			patch->texCoord(v, 2) = patch->texCoord(v, 0);
		}
		else
		{
			patch->texCoord(v, 2) = pair->from()->tex();
		}
	}
	//nvCheck (v == patch->edgeCount());
}


//-------------------------------
//  TangentSpacePatch
//-------------------------------

// @@ Is it necessary to do something this complicated?
static Vector2 computeLinearDependence(Vector3::Arg C, Vector3::Arg A, Vector3::Arg B)
{
	// Solve this equation in the least squares sense:

	// alpha A + beta B = C

	// [A B] [alpha beta]t = [C]
	
	// [A B]t[A B] [alpha beta]t = [A B]t[C]

	// [AA AB] [alpha] = [AC]
	// [AB BB] [beta ]   [BC]

	// [alpha] = [ BB -AB] [AC] / (AA * BB - AB * AB)
	// [beta ]   [-BA  AA] [BC]

	float AA = dot(A, A);
	float AB = dot(A, B);
	float BB = dot(B, B);

	float AC = dot(A, C);
	float BC = dot(B, C);

	float idet = 1.0f / (AA * BB - AB * AB);

	float alpha = (BB * AC - AB * BC) * idet;
	float beta  = (AA * BC - AB * AC) * idet;

	return Vector2(alpha, beta);
}

void AccMeshBuilder::computeChartTangents(BezierAccPatch * patch, const Array<Basis> & basisArray) const
{
	nvDebugCheck(patch != NULL);

	const static int tangentIdx[] = {0, 8, 11, 3};
	const static int bitangentIdx[] = {0, 3, 11, 8};

	const HalfEdge::Face * face = patch->face();
	nvDebugCheck(face != NULL);
	nvDebugCheck(face->edgeCount() == 4);

    const HalfEdge::Edge * firstEdge = patch->faceOneRing()->firstEdge();
    nvDebugCheck(firstEdge != NULL);

	uint v  = 0;
	for (HalfEdge::Face::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance(), v++)
	{
		const HalfEdge::Edge * edge = it.current();
		const HalfEdge::Vertex * vertex = edge->from();

		Vector3 tangent = basisArray[vertex->id()].tangent;
		Vector3 bitangent = basisArray[vertex->id()].bitangent;

	//	m_chartTangentArray[v] = computeLinearDependence(tangent, normalize(patch->tangent(tangentIdx[v])), normalize(patch->bitangent(bitangentIdx[v]));
	//	m_chartBitangentArray[v] = computeLinearDependence(bitangent, normalize(patch->tangent(tangentIdx[v])), normalize(patch->bitangent(bitangentIdx[v]));
	}
}


//-------------------------------
//  Shared methods
//-------------------------------

void AccMeshBuilder::computeGlobalTangentSpace(Array<Basis> & basisArray) const
{
	// @@ None of these methods is perfect, but this is the best one that I've found so far.

	//geometry::computeMeshTangents(m_mesh, basisArray);
	geometry::computeCatmullClarkTangents(m_mesh, basisArray);
	//geometry::computeConformalMeshTangents(m_mesh, basisArray);
	//geometry::computeConformalCatmullClarkTangents(m_mesh, basisArray);

	// @@ CharlesB suggests several ideas:

	// My main avenue of approach would be to do a global lsqr optimization of the tangents
	// on each chart.  You can minimize some error metric.  It's possible you could just 
	// completely ignore the UV's and optimize the tangents directly.  Or perhaps seed the 
	// solution with the tangents from the UV's to start and then just optimize from there.

	// You can put a few different terms in your error.  For example you could make error 
	// terms from the limit surface directly.  eg. you can do error terms that are per-triangle
	// rather than per-vertex.

	// This is just some random ideas but maybe a term for the 3 tangents on the triangle 
	// vertices to be as similar as possible, eg. something with the dot product of the 
	// tangents and also of the bitangents.

	// Another idea would be to take the average of the 3 tangents and compare that to the 
	// true limit surface normal at the center of the triangle and you want them to match.

	// In addition to that you could do the same thing per edge, take the average of the two 
	// tangents across each edge - make an error term for how well it matches the normal at 
	// the middle of that edge.  Or an error for term the tangents having good dot products 
	// across an edge.

	// Per triangle terms can be weighted by triangle area.  Per edge terms I'm not sure how to 
	// weight correctly.

	// Another idea :

	// Make the tangent space from the UV's.  Presumably this should be pretty close to good.

	// Do an iterative relaxation on the tangents to try to get them to relax more.  There are 
	// a few options here, maybe some kind of simulated annealing type relaxation with randomized
	// changes that get smaller over time.

	// Another total random idea I just had :

	// You could use some kind of semi-physical system.

	// Think of each tangent space as a magnet.  They apply magnetic forces to other nearby 
	// tangents to try to make them align the same way.  Just seed the sim and then let the 
	// physics evolve.

	// Another one would be that the tangents are sort of like the flow vectors for a fluid.  
	// Maybe imagine that some fluid is flowing over the surface of the mesh, the tangent is 
	// the direction of the velocity of the fluid.  Laminar flow will be lower energy so direction 
	// changes should go away.
}

void AccMeshBuilder::computeBoundaryTangentStencils(const AccPatch * patch, const HalfEdge::Vertex * vertex, StencilMask & r0, StencilMask & r1) const
{
	nvCheck(vertex->isBoundary());
	nvDebugCheck(r0.count() == patch->vertexCount());
	nvDebugCheck(r1.count() == patch->vertexCount());

	const HalfEdge::Edge * edge0 = vertex->edge();
	nvCheck(edge0->face() == NULL);
	nvCheck(edge0->to() != vertex);

	const HalfEdge::Edge * edgek = vertex->edge()->prev();
	nvCheck(edgek->face() == NULL);
	nvCheck(edgek->from() != vertex);

	const uint valence = vertex->valence();

	const int k = valence - 1;
	const float omega = PI / k;
	const float s = sinf(omega);
	const float c = cosf(omega);

	const float factor = 1.0f / (3 * k + c);

	const float gamma = -4 * s * factor;
	r0[patch->vertexIndex(vertex)] = gamma;
	//r1[patch->vertexIndex(vertex)] = 0;

	const float salpha0 = -((1 + 2 * c) * sqrtf(1 + c)) * factor / sqrtf(1 - c);
	const float calpha0 = 1.0f / 2.0f;

	r0[patch->vertexIndex(edge0->to())] = salpha0;
	r1[patch->vertexIndex(edge0->to())] = calpha0;

	const float salphak = salpha0;
	const float calphak = -1.0f / 2.0f;

	r0[patch->vertexIndex(edgek->from())] = salphak;
	r1[patch->vertexIndex(edgek->from())] = calphak;

	int j = 0;
	for (HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), j++)
	{
		const HalfEdge::Edge * edge = it.current();

		if (j == k) break;

		const HalfEdge::Vertex * p = edge->to();
		const HalfEdge::Vertex * q = edge->pair()->prev()->from();

		const float alphaj = 4 * sinf(j * omega) * factor;
		const float betaj = (sinf(j * omega) + sinf((j + 1) * omega)) * factor;

		if (j != 0)
		{
			r0[patch->vertexIndex(p)] += alphaj;
		}

		if (edge->pair()->prev()->prev()->prev() == edge->pair())
		{
			r0[patch->vertexIndex(vertex)] += betaj / 3.0f;
			r0[patch->vertexIndex(edge->pair()->from())] += betaj / 3.0f;
			r0[patch->vertexIndex(q)] += betaj / 3.0f;
		}
		else
		{
			r0[patch->vertexIndex(q)] += betaj;
		}
	}

    if (valence == 2)
    {
        // r0 perpendicular to r1
	    r0[patch->vertexIndex(vertex)] = -4.0f / 3.0f;
	    r0[patch->vertexIndex(edgek->from())] = 1.0f / 2.0f;
	    r0[patch->vertexIndex(edge0->to())] = 1.0f / 2.0f;
	    r0[patch->vertexIndex(edge0->next()->to())] = 1.0f / 3.0f;
    }
}
