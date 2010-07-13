// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include "AccPatch.h"
#include "AccMesh.h"
//#include "RemapFaces.h"

#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>

#include <nvmath/Basis.h>
#include <nvmath/Bezier.h>

#include <nvcore/Radix.h>


using namespace nv;

namespace
{
	// Set bits corresponding to the edges and corners that this face owns.
	// @@ Ideally, this should be done at the chart level instead of the face level.
	static uint computeBoundaryMask(const HalfEdge::Face * face)
	{
		nvDebugCheck(face->edgeCount() == 4);
		
		uint mask = 0;
		
		uint i = 0;
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++)
		{
			const HalfEdge::Edge * edge = it.current();
			
			bool bit = (edge->id() > edge->pair()->id());
			
			mask |= bit << i;
		}

		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++)
		{
			const HalfEdge::Vertex * vertex = it.current()->vertex();
			
			bool bit = true;
			
			for (HalfEdge::Vertex::ConstVertexIterator cit(vertex->colocals()); !cit.isDone(); cit.advance())
			{
				const HalfEdge::Vertex * colocal = cit.current();

				if (vertex != colocal && vertex->id() < colocal->id())
				{
					bit = false;
					break;
				}
			}
			
			mask |= bit << i;
		}
		
		nvDebugCheck(i == 8);
		
		return mask;
	}
}


//-------------------------------
//  FaceOneRing
//-------------------------------

FaceOneRing::FaceOneRing(const HalfEdge::Face * face, const HalfEdge::Edge * firstEdge, FaceTopology * faceTopology) : m_face(face), m_firstEdge(firstEdge), m_faceTopology(faceTopology)
{
	nvDebugCheck(m_face != NULL);

	m_edgeCount = face->edgeCount();
	nvCheck (m_edgeCount == 3 || m_edgeCount == 4);

	initVertices();

	// @@ This should not be an assertion. It shouldbe handled more gracefully.
	nvCheck (vertexCount() < 32);

	if (faceTopology->isRegular) {
		// @@ Fixme: Handle self-overlapping stencils properly.
		nvCheck (vertexCount() == 16);
	}
}

uint FaceOneRing::vertexCount() const
{
	return m_vertexArray.count();
}

const HalfEdge::Vertex * FaceOneRing::vertexAt(uint i) const
{
    return m_vertexArray[i];
}

uint FaceOneRing::vertexIndexAt(uint i) const
{
	return m_vertexIndexArray[i];
}

uint FaceOneRing::vertexIndex(const HalfEdge::Vertex * vertex) const
{
	nvDebugCheck(vertex != NULL);

	// @@ This code need to handle overlapping stencils. How?

	const uint count = this->vertexCount();
	for (uint i = 0; i < count; i++)
	{
		if (/*m_vertexArray[i] == vertex ||*/ m_vertexArray[i]->pos() == vertex->pos()) 
		{
			return i;
		}
	}
	
	nvCheck(0);
	return 0;
}

Vector3 FaceOneRing::evalutePositionStencil(const StencilMask & mask) const
{
	nvDebugCheck(mask.count() == vertexCount());
	
	Vector3 p(zero);
	
	const uint count = vertexCount();
	for (uint i = 0; i < count; i++)
	{
		// Always add vertices in the same order.
		uint idx = m_vertexIndexArray[i];
		p += m_vertexArray[idx]->pos() * mask[idx];
	}
	
	return p;
}

Vector3 FaceOneRing::evaluteNormalStencil(const StencilMask & mask) const
{
	nvDebugCheck(mask.count() == vertexCount());
	
	Vector3 normal(zero);
	
	const uint count = vertexCount();
	for (uint i = 0; i < count; i++)
	{
		// Always add vertices in the same order.
		uint idx = m_vertexIndexArray[i];
		normal += m_vertexArray[idx]->nor() * mask[idx];
	}

    // @@ Implement linear interpolation using exponential maps?

	return normalizeSafe(normal, Vector3(zero), 0.0f);
}

bool FaceOneRing::isTriangle() const
{
	return m_edgeCount == 3;
}

bool FaceOneRing::isQuad() const
{
	return m_edgeCount == 4;
}

uint FaceOneRing::edgeCount() const
{
	return m_edgeCount;
}

/*static*/ bool FaceOneRing::isRegular(const HalfEdge::Face * face)
{
	if (!isQuad(face)) return false;

	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();

		// Regular faces don't have boundary edges.
		if (edge->isBoundary()) return false;

		// Regular faces only have ordinary vertices.
		if (edge->vertex()->valence() != 4) return false;

		// Regular faces are only adjacent to quads
		for (HalfEdge::Vertex::ConstEdgeIterator eit(edge); !eit.isDone(); eit.advance())
		{
			if (eit.current()->face() == NULL || eit.current()->face()->edgeCount() != 4) return false;
		}
	}

	return true;
}

/*static*/ bool FaceOneRing::isPolarTriangle(const HalfEdge::Face * face)
{
	if (!isTriangle(face)) return false;

	const HalfEdge::Vertex * pole = NULL;

	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();

		uint valence = edge->vertex()->valence();

		if (valence != 4)	// @@ != 4 or > 4?
		{
			if (pole != NULL) return false;

			pole = edge->vertex();
		}
	}

	if (pole == NULL) return false;

	// Make sure the pole is really a pole.
	for (HalfEdge::Vertex::ConstEdgeIterator it(pole->edges()); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();

		// Surrounding edges must have valence 4.
		if (edge->to()->valence() != 4) return false;

		// Edge is not boundary.
		if (edge->face() == NULL) return false;

		// Surrounding faces are triangles.
		if (!isTriangle(edge->face())) return false;

		// Surrounding faces are not on the boundaries.
		if (edge->face()->isBoundary()) return false;
	}

	return true;
}

/*static*/ bool FaceOneRing::isTriangle(const HalfEdge::Face * face)
{
	nvDebugCheck(face != NULL);
	return face->edgeCount() == 3;
}
/*static*/ bool FaceOneRing::isQuad(const HalfEdge::Face * face)
{
	nvDebugCheck(face != NULL);
	return face->edgeCount() == 4;
}

/*static*/ bool FaceOneRing::isPent(const HalfEdge::Face * face)
{
	nvDebugCheck(face != NULL);
	return face->edgeCount() == 5;
}

/*static*/ bool FaceOneRing::isBoundary(const HalfEdge::Face * face)
{
	nvDebugCheck(face != NULL);
	// Note that face->isBoundary() returns a different result. That function returns true when any of the *edges* are on the boundary.
	// However, this function returns true if any of the face *vertices* are on the boundary. 

	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();
		const HalfEdge::Vertex * vertex = edge->vertex();

		if (vertex->isBoundary()) return true;
	}

	return false;
}


void FaceOneRing::initVertices()
{
	const uint edgeCount = m_face->edgeCount();
	
	m_vertexArray.reserve(16);

	// Add face vertices.
	for (HalfEdge::Face::ConstEdgeIterator it(m_firstEdge); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();
		
		m_vertexArray.append(edge->from());
	}
	
	// @@ Add support for non manifold surfaces!
	// The fix: 
	// - not all colocals should point to the same edge.
	// - multiple colocals could belong to different boundaries, make sure they point to the right one.

	// @@ When the face neighborhood wraps that could result in overlapping stencils. 

	// Add surronding vertices.
	for (HalfEdge::Face::ConstEdgeIterator it(m_firstEdge); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * firstEdge = it.current();
		const HalfEdge::Vertex * vertex = firstEdge->from();
	
		const uint valence = vertex->valence();
		uint i = 0;

		// Traverse edges around vertex
		for (HalfEdge::Vertex::ReverseConstEdgeIterator eit(firstEdge); !eit.isDone(); eit.advance(), i++)
		{
			const HalfEdge::Edge * edge = eit.current();

			nvCheck(edge->from()->pos() == vertex->pos());

			appendVertex(edge->to());

			if (edge->face() != NULL)
			{
				appendVertex(edge->next()->to());
			}
		}
		nvDebugCheck(i == valence);
	}

	const uint vertexCount = m_vertexArray.count();

	// @@ You can only sort by id when the vertices do not have colocals.
#if 0
	// Sort vertices by id. Create index array per patch.
	Array<uint> vertexIdArray;
	vertexIdArray.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		vertexIdArray[v] = m_vertexArray[v]->id();
	}

	RadixSort radix;
	radix.sort(vertexIdArray);
#else
	// Sort vertices lexycographically. Create index array per patch.
	Array<float> vertexXArray;
	Array<float> vertexYArray;
	Array<float> vertexZArray;

	vertexXArray.resize(vertexCount);
	vertexYArray.resize(vertexCount);
	vertexZArray.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		vertexXArray[v] = m_vertexArray[v]->pos().x();
		vertexYArray[v] = m_vertexArray[v]->pos().y();
		vertexZArray[v] = m_vertexArray[v]->pos().z();
	}

	RadixSort radix;
	radix.sort(vertexXArray).sort(vertexYArray).sort(vertexZArray);
#endif

	m_vertexIndexArray.resize(vertexCount);

	const uint * indices = radix.ranks();
	for (uint v = 0; v < vertexCount; v++)
	{
		m_vertexIndexArray[v] = indices[v];
	}
}


void FaceOneRing::appendVertex(const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(vertex != NULL);
	
	// Avoid having duplicate verts.
	if (!hasVertex(vertex))
	{
		m_vertexArray.append(vertex);
	}
}


bool FaceOneRing::hasVertex(const HalfEdge::Vertex * vertex) const
{
	nvDebugCheck(vertex != NULL);
	
	const uint vertexCount = m_vertexArray.count();
	for (uint i = 0; i < vertexCount; i++)
	{
		if (m_vertexArray[i]->pos() == vertex->pos())
		{
			return true;
		}
	}

	return false;
}


//-------------------------------
//  AccPatch 
//-------------------------------

AccPatch::AccPatch(const FaceOneRing * face) : m_face(face)
{
	nvDebugCheck(face != NULL);
}

/*virtual*/ AccPatch::~AccPatch()
{
}


//-------------------------------
//  SurfacePatch 
//-------------------------------

SurfacePatch::SurfacePatch(const FaceOneRing * face) : AccPatch(face)
{
}


//-------------------------------
//  TexCoordPatch 
//-------------------------------

TexCoordPatch::TexCoordPatch(const FaceOneRing * face) : AccPatch(face)
{
}

const Vector2 TexCoordPatch::texCoord(uint v, uint i) const
{
	nvDebugCheck(v < m_face->edgeCount() && i < 4);
	return m_texCoordArray[4 * v + i];
}

Vector2 & TexCoordPatch::texCoord(uint v, uint i)
{
	nvDebugCheck(v < m_face->edgeCount() && i < 4);
	return m_texCoordArray[4 * v + i];
}


//-------------------------------
//  TangentSpacePatch 
//-------------------------------

Vector3 TangentSpacePatch::chartTangent(uint idx) const
{
	nvDebugCheck(idx < 4);
	/*
	const static uint tangentIndex[4] = { 0, 3, 8, 11 };
	const static uint bitangentIndex[4] = { 0, 8, 3, 11 };
	
	Vector3 T = m_tangentArray[tangentIndex[idx]];
	Vector3 B = m_bitangentArray[bitangentIndex[idx]];

	Vector2 t = m_chartTangentArray[idx];
	return normalize(T * t.x() + B * t.y());
	*/
	return Vector3(zero);
}

Vector3 TangentSpacePatch::chartBitangent(uint idx) const
{
	nvDebugCheck(idx < 4);
	/*
	const static uint tangentIndex[4] = { 0, 3, 8, 11 };
	const static uint bitangentIndex[4] = { 0, 8, 3, 11 };
	
	Vector3 T = m_tangentArray[tangentIndex[idx]];
	Vector3 B = m_bitangentArray[bitangentIndex[idx]];

	Vector2 t = m_chartBitangentArray[idx];
	return normalize(T * t.x() + B * t.y());
	*/
	return Vector3(zero);
}

Vector3 TangentSpacePatch::evaluateChartTangent(float u, float v) const
{
	Vector2 L0(1-u,1-v);
	Vector2 L1(u,v);

#if 0
	// @@ This interpolation mode is more efficient, but is not continuous.
	Vector2 t =
		(L0.x() * m_chartTangent[0] + L1.x() * m_chartTangent[1]) * L0.y() + 
		(L0.x() * m_chartTangent[2] + L1.x() * m_chartTangent[3]) * L1.y();

	Vector3 T = normalize(evaluateTangent(u, v));
	Vector3 B = normalize(evaluateBitangent(u, v));

	return normalize(T * t.x() + B * t.y());
#else
	Vector3 tangent = 
		(L0.x() * chartTangent(0) + L1.x() * chartTangent(1)) * L0.y() + 
		(L0.x() * chartTangent(2) + L1.x() * chartTangent(3)) * L1.y();

	return tangent;
#endif
}

Vector3 TangentSpacePatch::evaluateChartBitangent(float u, float v) const
{
	Vector2 L0(1-u,1-v);
	Vector2 L1(u,v);

#if 0 
	// @@ This interpolation mode is more efficient, but is not continuous.
	Vector2 t =
		(L0.x() * m_chartBitangent[0] + L1.x() * m_chartBitangent[1]) * L0.y() + 
		(L0.x() * m_chartBitangent[2] + L1.x() * m_chartBitangent[3]) * L1.y();

	Vector3 T = normalize(evaluateTangent(u, v));
	Vector3 B = normalize(evaluateBitangent(u, v));

	return normalize(T * t.x() + B * t.y());
#else
	Vector3 bitangent = 
		(L0.x() * chartBitangent(0) + L1.x() * chartBitangent(1)) * L0.y() + 
		(L0.x() * chartBitangent(2) + L1.x() * chartBitangent(3)) * L1.y();

	return bitangent;
#endif
}
/*
Vector3 GregoryAccPatch::chartTangent(uint idx) const
{
	nvDebugCheck(idx < 4);
	
	const static float u[4] = { 0, 1, 0, 1 };
	const static float v[4] = { 0, 0, 1, 1 };
	
	Vector3 T = evaluateTangent(u[idx], v[idx]);
	Vector3 B = evaluateBitangent(u[idx], v[idx]);

	Vector2 t = m_chartTangentArray[idx];	
	return normalize(T * t.x() + B * t.y());	
}

Vector3 GregoryAccPatch::chartBitangent(uint idx) const
{
	nvDebugCheck(idx < 4);

	const static float u[4] = { 0, 1, 0, 1 };
	const static float v[4] = { 0, 0, 1, 1 };
	
	Vector3 T = evaluateTangent(u[idx], v[idx]);
	Vector3 B = evaluateBitangent(u[idx], v[idx]);

	Vector2 t = m_chartBitangentArray[idx];	
	return normalize(T * t.x() + B * t.y());	
}
*/


//-------------------------------
//  BezierAccStencil
//-------------------------------

BezierAccStencil::BezierAccStencil(const FaceOneRing * face) : m_face(face)
{
	const uint vertexCount = m_face->vertexCount();

	// Resize stencils.
	for (int i = 0; i < 16; i++)
	{
		m_positionStencil[i].resize(vertexCount);
	}

	for (int i = 0; i < 12; i++)
	{
		m_tangentStencil[i].resize(vertexCount);
		m_bitangentStencil[i].resize(vertexCount);
	}
/*
	for (int i = 0; i < 8; i++)
	{
		m_symmetricTangentStencil[i].resize(vertexCount);
		m_symmetricBitangentStencil[i].resize(vertexCount);
	}

	for (int i = 0; i < 8; i++)
	{
		m_boundaryTangentStencil[i].resize(vertexCount);
		m_boundaryBitangentStencil[i].resize(vertexCount);
	}
	
	for (int i = 0; i < 4; i++)
	{
		m_cornerTangentStencil[i].resize(vertexCount);
		m_cornerBitangentStencil[i].resize(vertexCount);
	}
*/
}

const StencilMask & BezierAccStencil::positionStencil(uint i) const
{
	nvDebugCheck(i < 16);
	return m_positionStencil[i];
}
StencilMask & BezierAccStencil::positionStencil(uint i)
{
	nvDebugCheck(i < 16);
	return m_positionStencil[i];
}
float & BezierAccStencil::positionStencil(uint i, const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(i < 16);
	nvDebugCheck(vertex != NULL);
	uint vi = m_face->vertexIndex(vertex);
	return m_positionStencil[i][vi];
}

const StencilMask & BezierAccStencil::tangentStencil(uint i) const
{
	nvDebugCheck(i < 12);
	return m_tangentStencil[i];
}
StencilMask & BezierAccStencil::tangentStencil(uint i)
{
	nvDebugCheck(i < 12);
	return m_tangentStencil[i];
}
float & BezierAccStencil::tangentStencil(uint i, const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(i < 12);
	nvDebugCheck(vertex != NULL);
	uint vi = m_face->vertexIndex(vertex);
	return m_tangentStencil[i][vi];
}

const StencilMask & BezierAccStencil::bitangentStencil(uint i) const
{
	nvDebugCheck(i < 12);
	return m_bitangentStencil[i];
}
StencilMask & BezierAccStencil::bitangentStencil(uint i)
{
	nvDebugCheck(i < 12);
	return m_bitangentStencil[i];
}
float & BezierAccStencil::bitangentStencil(uint i, const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(i < 12);
	nvDebugCheck(vertex != NULL);
	uint vi = m_face->vertexIndex(vertex);
	return m_bitangentStencil[i][vi];
}

/*
const StencilMask & BezierAccStencil::symmetricTangentStencil(uint i) const
{
	nvDebugCheck(i < 8);
	return m_symmetricTangentStencil[i];
}
StencilMask & BezierAccStencil::symmetricTangentStencil(uint i)
{
	nvDebugCheck(i < 8);
	return m_symmetricTangentStencil[i];
}
float & BezierAccStencil::symmetricTangentStencil(uint i, const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(vertex != NULL);
	uint vi = vertexIndex(vertex);
	return symmetricTangentStencil(i)[vi];
}

const StencilMask & BezierAccStencil::symmetricBitangentStencil(uint i) const
{
	nvDebugCheck(i < 8);
	return m_symmetricBitangentStencil[i];
}
StencilMask & BezierAccStencil::symmetricBitangentStencil(uint i)
{
	nvDebugCheck(i < 8);
	return m_symmetricBitangentStencil[i];
}
float & BezierAccStencil::symmetricBitangentStencil(uint i, const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(vertex != NULL);
	uint vi = vertexIndex(vertex);
	return symmetricBitangentStencil(i)[vi];
}
*/

bool BezierAccStencil::operator==(const BezierAccStencil & other) const
{
	for (int i = 0; i < 16; i++)
	{
		if (this->m_positionStencil[i] != other.m_positionStencil[i]) 
			return false;
	}

	for (int i = 0; i < 12; i++)
	{
		if (this->m_tangentStencil[i] != other.m_tangentStencil[i])
			return false;
		if (this->m_bitangentStencil[i] != other.m_bitangentStencil[i])
			return false;
	}

	return true;
}


//-------------------------------
//  BezierAccPatch 
//-------------------------------

BezierAccPatch::BezierAccPatch(const FaceOneRing * face) : SurfacePatch(face)
{
	m_stencil = new BezierAccStencil(face);
}

BezierAccPatch::~BezierAccPatch()
{
	delete m_stencil;
}

void BezierAccPatch::evaluateControlPoints()
{
	// Evaluate control points from stencils.
	for (uint i = 0; i < 16; i++)
	{
		m_positionArray[i] = m_face->evalutePositionStencil(m_stencil->positionStencil(i));
        m_normalArray[i] = m_face->evaluteNormalStencil(m_stencil->positionStencil(i));
	}
	for (uint i = 0; i < 12; i++)
	{
		m_tangentArray[i] = m_face->evalutePositionStencil(m_stencil->tangentStencil(i));
		m_bitangentArray[i] = m_face->evalutePositionStencil(m_stencil->bitangentStencil(i));
	}
	/*for (uint i = 0; i < 8; i++)
	{
		m_symmetricTangentArray[i] = m_face->evalutePositionStencil(m_stencil->symmetricTangentStencil(i));
		m_symmetricBitangentArray[i] = m_face->evalutePositionStencil(m_stencil->symmetricBitangentStencil(i));
	}
	for (uint i = 0; i < 8; i++)
	{
		m_boundaryTangentArray[i] = m_face->evalutePositionStencil(m_stencil->boundaryTangentStencil(i));
		m_boundaryBitangentArray[i] = m_face->evalutePositionStencil(m_stencil->boundaryBitangentStencil(i));
	}
	for (uint i = 0; i < 4; i++)
	{
		m_cornerTangentArray[i] = m_face->evalutePositionStencil(m_stencil->cornerTangentStencil(i));
		m_cornerBitangentArray[i] = m_face->evalutePositionStencil(m_stencil->cornerBitangentStencil(i));
	}*/
	/*for (uint i = 0; i < 4; i++)
	{
		m_normalArray[i] = normalizeSafe(cross(m_stencil->cornerTangentArray(i), m_stencil->cornerBitangentArray(i)), Vector3(zero), 0.0f);
	}*/
}

Vector3 BezierAccPatch::position(uint idx) const
{
	nvDebugCheck(idx < 16);
	return m_positionArray[idx];
}

Vector3 BezierAccPatch::tangent(uint idx) const
{
	nvDebugCheck(idx < 12);
	return m_tangentArray[idx];
}

Vector3 BezierAccPatch::bitangent(uint idx) const
{
	nvDebugCheck(idx < 12);
	return m_bitangentArray[idx];
}

Vector3 BezierAccPatch::normal(uint idx) const
{
	nvDebugCheck(idx < 16);
	return m_normalArray[idx];
}

/*
Vector3 BezierAccPatch::boundaryTangent(uint idx) const
{
	nvDebugCheck(idx < 8);
	return m_boundaryTangentArray[idx];
}

Vector3 BezierAccPatch::boundaryBitangent(uint idx) const
{
	nvDebugCheck(idx < 8);
	return m_boundaryBitangentArray[idx];
}
*/
/*Vector3 BezierAccPatch::normal(uint idx) const
{
	nvDebugCheck(idx < 4);
	
	const static uint tangentIndex[4] = { 0, 3, 8, 11 };
	const static uint bitangentIndex[4] = { 0, 8, 3, 11 };
	
	// @@ Use chart or corner tangents to obtain consistent results?

	Vector3 T = m_tangentArray[tangentIndex[idx]];
	Vector3 B = m_bitangentArray[bitangentIndex[idx]];

	return normalizeSafe(cross(T, B), Vector3(zero), 0.0f);
}*/

Vector3 BezierAccPatch::evaluateTangent(float u, float v) const
{
	Vector2 L0(1-u,1-v);
	Vector2 L1(u,v);
	
	Vector2 Q0 =     L0 * L0;
	Vector2 Q1 = 2 * L0 * L1;
	Vector2 Q2 =     L1 * L1;

	Vector2 B0 =     L0 * Q0;
	Vector2 B1 = 3 * Q0 * L1;
	Vector2 B2 = 3 * Q2 * L0;
	Vector2 B3 =     Q2 * L1;

	const Vector3 * pu = m_tangentArray;
	return 
		(pu[0 ] * B0.y() + pu[1 ] * B1.y() + pu[2 ] * B2.y() + pu[3 ] * B3.y()) * Q0.x() +
		(pu[4 ] * B0.y() + pu[5 ] * B1.y() + pu[6 ] * B2.y() + pu[7 ] * B3.y()) * Q1.x() +
		(pu[8 ] * B0.y() + pu[9 ] * B1.y() + pu[10] * B2.y() + pu[11] * B3.y()) * Q2.x();
}

Vector3 BezierAccPatch::evaluateBitangent(float u, float v) const
{
	Vector2 L0(1-u,1-v);
	Vector2 L1(u,v);
	
	Vector2 Q0 =     L0 * L0;
	Vector2 Q1 = 2 * L0 * L1;
	Vector2 Q2 =     L1 * L1;

	Vector2 B0 =     L0 * Q0;
	Vector2 B1 = 3 * Q0 * L1;
	Vector2 B2 = 3 * Q2 * L0;
	Vector2 B3 =     Q2 * L1;
	
	const Vector3 * pv = m_bitangentArray;
	return 
		(pv[0 ] * B0.x() + pv[1 ] * B1.x() + pv[2 ] * B2.x() + pv[3 ] * B3.x()) * Q0.y() +
		(pv[4 ] * B0.x() + pv[5 ] * B1.x() + pv[6 ] * B2.x() + pv[7 ] * B3.x()) * Q1.y() +
		(pv[8 ] * B0.x() + pv[9 ] * B1.x() + pv[10] * B2.x() + pv[11] * B3.x()) * Q2.y();
}

/*virtual*/ void BezierAccPatch::evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const
{
	nvDebugCheck(pos != NULL);

	evaluatePosition(u, v, pos);

	*du = evaluateTangent(u,v);
	*dv = evaluateBitangent(u,v);
}

/*virtual*/ void BezierAccPatch::evaluatePosition(float u, float v, Vector3 * pos) const
{
	nvDebugCheck(pos != NULL);

	Vector2 L0(1-u,1-v);
	Vector2 L1(u,v);

	Vector2 B0 =     (L0 * L0) * L0;
	Vector2 B1 = 3 * (L0 * L0) * L1;
	Vector2 B2 = 3 * (L1 * L1) * L0;
	Vector2 B3 =     (L1 * L1) * L1;

	*pos = 
		(B0.x() * p(0,0) + B1.x() * p(1,0) + B2.x() * p(2,0) + B3.x() * p(3,0)) * B0.y() +
		(B0.x() * p(0,1) + B1.x() * p(1,1) + B2.x() * p(2,1) + B3.x() * p(3,1)) * B1.y() +
		(B0.x() * p(0,2) + B1.x() * p(1,2) + B2.x() * p(2,2) + B3.x() * p(3,2)) * B2.y() +
		(B0.x() * p(0,3) + B1.x() * p(1,3) + B2.x() * p(2,3) + B3.x() * p(3,3)) * B3.y();
}

/*virtual*/ void BezierAccPatch::evaluatePatchFrame(float u, float v, Basis * frame) const
{
	nvDebugCheck(frame != NULL);

	frame->tangent = evaluateTangent(u,v);
	frame->bitangent = evaluateBitangent(u,v);

#if 1
	// Accurate surface normal:
	frame->normal = normalizeSafe(cross(frame->tangent, frame->bitangent), Vector3(zero), 0);
#else
	// Linearly interpolated normal:
	Vector2 L0(1-u,1-v);
	Vector2 L1(u,v);

	frame->normal = 
		(L0.x() * normal(0) + L1.x() * normal(1)) * L0.y() + 
		(L0.x() * normal(2) + L1.x() * normal(3)) * L1.y();
#endif
}


Vector3 BezierAccPatch::p(uint u, uint v) const
{
	nvDebugCheck(u < 4 && v < 4);
	return m_positionArray[4 * v + u];
}

Vector3 BezierAccPatch::tu(uint u, uint v) const
{
	nvDebugCheck(u < 3 && v < 4);
	return m_tangentArray[4 * v + u]; // @@ Not sure this is correct.
}

Vector3 BezierAccPatch::tv(uint u, uint v) const
{
	nvDebugCheck(u < 4 && v < 3);
	return m_bitangentArray[4 * u + v]; // @@ Not sure this is correct.
}


//-------------------------------
//  GregoryAccStencil
//-------------------------------

GregoryAccStencil::GregoryAccStencil(const FaceOneRing * face) : m_face(face)
{
	const uint vertexCount = m_face->vertexCount();

	// Resize stencils.
	for (int i = 0; i < 20; i++)
	{
		m_positionStencil[i].resize(vertexCount);
	}
}

const StencilMask & GregoryAccStencil::positionStencil(uint i) const
{
	nvDebugCheck(i < 20);
	return m_positionStencil[i];
}

StencilMask & GregoryAccStencil::positionStencil(uint i)
{
	nvDebugCheck(i < 20);
	return m_positionStencil[i];
}

float & GregoryAccStencil::positionStencil(uint i, const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(i < 20);
	nvDebugCheck(vertex != NULL);
	uint vi = m_face->vertexIndex(vertex);
	return m_positionStencil[i][vi];
}

bool GregoryAccStencil::operator==(const GregoryAccStencil & other) const
{
	for (int i = 0; i < 20; i++)
	{
		if (this->m_positionStencil[i] != other.m_positionStencil[i])
			return false;
	}

	return true;
}


//-------------------------------
//  GregoryAccPatch 
//-------------------------------

GregoryAccPatch::GregoryAccPatch(const FaceOneRing * face) : SurfacePatch(face)
{
	m_stencil = new GregoryAccStencil(face);
}

GregoryAccPatch::~GregoryAccPatch()
{
	delete m_stencil;
}

void GregoryAccPatch::evaluateControlPoints()
{
	// Evaluate control points from stencils.
	for (uint i = 0; i < 20; i++)
	{
		m_positionArray[i] = m_face->evalutePositionStencil(m_stencil->positionStencil(i));
	}
}
	
Vector3 GregoryAccPatch::position(uint idx) const
{
	nvDebugCheck(idx < 20);
	return m_positionArray[idx];
}

/*virtual*/ void GregoryAccPatch::evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const
{
	if (isQuad())
	{
        evaluateQuadGregoryPatch(u, v, m_positionArray, pos, du, dv);
	}
	else
	{
        nvDebugCheck(isTriangle());
        evaluateTriangleGregoryPatch(u, v, m_positionArray, pos, du, dv);
	}
}


//-------------------------------
//  PolarPatch 
//-------------------------------

PolarPatch::PolarPatch(const FaceOneRing * face) : SurfacePatch(face)
{
}

/*virtual*/ void PolarPatch::evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const
{
    const Vector3 * p = m_positionArray;
    Vector3 q[16];

    for (int i = 0; i < 4; i++) {
        q[i] = p[0];
    }
    for (int i = 0; i < 12; i++) {
        q[4+i] = p[1+i];
    }

    evaluateCubicBezierPatch(u, v, q, pos, du, dv);

    if (u == 0 && v == 0) {
        // @@ Handle corner case.
        //if (du) *du = ...;
    }
}


//-------------------------------
//  PmAccPatch
//-------------------------------

PmQuadAccPatch::PmQuadAccPatch(const FaceOneRing * face) : SurfacePatch(face)
{
}

Vector3 PmQuadAccPatch::position(uint idx) const
{
	nvDebugCheck(idx < 24);
	return m_positionArray[idx];
}

/*virtual*/ void PmQuadAccPatch::evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const
{
    evaluateQuadPmPatch(u, v, m_positionArray, pos, du, dv);
}


//-------------------------------
//  RegularAccPatch
//-------------------------------

PmRegularAccPatch::PmRegularAccPatch(const FaceOneRing * face) : SurfacePatch(face)
{
}

Vector3 PmRegularAccPatch::position(uint idx) const
{
	nvDebugCheck(idx < 16);
	return m_positionArray[idx];
}

/*virtual*/ void PmRegularAccPatch::evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const
{
    evaluateCubicBezierPatch(u, v, m_positionArray, pos, du, dv);
}


//-------------------------------
//  PmTriangleAccPatch
//-------------------------------

PmTriangleAccPatch::PmTriangleAccPatch(const FaceOneRing * face) : SurfacePatch(face)
{
}
Vector3 PmTriangleAccPatch::position(uint idx) const
{
	nvDebugCheck(idx < 19);
	return m_positionArray[idx];
}

/*virtual*/ void PmTriangleAccPatch::evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const
{
    evaluateTrianglePmPatch(u, v, m_positionArray, pos, du, dv);
}
