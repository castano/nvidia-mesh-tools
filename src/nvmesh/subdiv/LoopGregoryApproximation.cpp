// Copyright NVIDIA Corporation 2007 -- Denis Kovacs <den.kovacs@gmail.com>

#include <nvmesh/subdiv/LoopGregoryApproximation.h>

#include <nvcore/Ptr.h>

#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>

#include <nvmesh/subdiv/CubicGregoryMesh.h>

using namespace nv;


LoopGregoryApproximation::LoopGregoryApproximation(const HalfEdge::Mesh * mesh, BoundaryMode boundaryMode)
{
	m_mesh = (HalfEdge::Mesh *)mesh;
	m_boundaryMode = boundaryMode;	// @@ Only mode supported so far!
}

LoopGregoryApproximation::~LoopGregoryApproximation()
{
}

static float cos_table[5] = { 1.0f, 1.0f/sqrtf(2.0f), 0.0f, -1.0f/sqrtf(2.0f), -1.0f };
static float GetCosAlpha(float a, float valence)
{
	a = fabsf(a);
	while (a > 2 * valence)
	{ a -= 2 * valence; }
	if (a > valence)
	{ a = 2 * valence - a; }
	a /= valence;
	float cos_pi_a;
	int quarter = int(a / 0.25f);
	if (quarter * 0.25f == a)
	{
		quarter = abs(quarter) % 8;
		if (quarter > 4) quarter = 8 - quarter;
		cos_pi_a = cos_table[quarter];
	}
	else
	{
		cos_pi_a = cosf(PI * a);
	}
	if (valence == 2)
		return cos_pi_a / 6;
	if (valence == 4)
		return cos_pi_a / 9;
	return cos_pi_a * (1.0f + 1 / sqrtf(4.0f / (cosf(PI / valence) * cosf(PI / valence)) + 1)) / (3.0f * valence);
}
static float GetCosBeta(float a, float valence)
{
	a = fabsf(a);
	while (a > 2 * valence)
	{ a -= 2 * valence; }
	if (a > valence)
	{ a = 2 * valence - a; }
	a /= valence;
	float cos_pi_a;
	int quarter = int(a / 0.25);
	if (quarter * 0.25 == a)
	{
		quarter = abs(quarter) % 8;
		if (quarter > 4) quarter = 8 - quarter;
		if (valence == 4)
		{
			if (quarter == 1)
				return  1.f/36;
			if (quarter == 3)
				return -1.f/36;
		}
		cos_pi_a = cos_table[quarter];
	}
	else
	{
		cos_pi_a = cosf(PI * a);
	}
	if (valence == 2)
		return cos_pi_a / 12;
	return cos_pi_a / (3.0f * valence * sqrtf(4.0f + cosf(PI/valence)*cosf(PI/valence)));
}

CubicGregoryMesh * LoopGregoryApproximation::buildMesh(LoopGregoryApproximation::Mode mode)
{
	const uint faceCount = m_mesh->faceCount();

	// Count patches.
	uint patchCount = 0;
	for(uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = m_mesh->faceAt(f);

		const uint edgeCount = face->edgeCount();
		nvCheck(edgeCount == 3 || edgeCount == 4);

		if ((edgeCount == 3 && mode == QuadsOnly) || (edgeCount == 4 && mode == TrianglesOnly)) {
			continue;
		}

		if (m_boundaryMode != BoundaryMode_None || !face->isBoundary()) {
			patchCount++;
		}
	}

	AutoPtr<CubicGregoryMesh> gregoryMesh(new CubicGregoryMesh(patchCount));

	//                                         /--- quads ---/   /--- tris ---/
	static const uint cornerIndices[7]      = {  8, 9, 10, 11,     6,  7,  8};
	static const uint corner2Indices[7]     = { 9, 10, 11,  8,     7,  8,  6};
	static const uint edge1Indices[7]       = {13, 15, 17, 19,    10, 12, 14};
	static const uint cornerEdge2Indices[7] = {12, 14, 16, 18,     9, 11, 13};
	static const uint interior1Indices[7]   = { 1,  3,  6,  4,     1,  3,  5};
	static const uint interior2Indices[7]   = { 2,  7,  5,  0,     2,  4,  0};

	// Compute control points.
	uint p = 0;
	for(uint f = 0; f < faceCount; f++)
	{
		bool bWrongPatch = false;
		const HalfEdge::Face * face = m_mesh->faceAt(f);

		const uint edgeCount = face->edgeCount();
		if ((edgeCount == 3 && mode == QuadsOnly) || (edgeCount == 4 && mode == TrianglesOnly)) {
			continue;
		}
		if (face->isBoundary())
		{
			printf("Face %d is on the boundary - skipping\n", f);
			continue;
		}

		CubicGregoryMesh::Patch &patch = gregoryMesh->patchAt(p++);
		patch.stencils.Clear();
		patch.face_id = face->id();

		patch.isQuad = (edgeCount == 4);
		uint primitiveOffset = patch.isQuad ? 0 : 4;

		// find an edge from which to start
		int nVerts = patch.isQuad ? 4 : 3;
		const HalfEdge::Edge *pFirstEdge = NULL;
		patch.stencils.bRegular = patch.isQuad; // at first we assume patch is regular
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			if (pFirstEdge == NULL)
			{
				pFirstEdge = it.current();
				continue;
			}
			const HalfEdge::Edge *e0 = pFirstEdge;
			const HalfEdge::Edge *e1 = it.current();
			for (int ie = 0; ie < nVerts; ++ie, e0 = e0->next(), e1 = e1->next())
			{
				if (e0->from()->valence() > e1->from()->valence())
					break;
				if (e0->from()->valence() < e1->from()->valence())
				{
					pFirstEdge = it.current();
					break;
				}
			}
		}

		// in stencils we want to have corner vertices first (so that we can compute LOD in tcs)
		patch.stencils[pFirstEdge->from()->pos()].AddWeight(0, 0);
		patch.stencils[pFirstEdge->next()->from()->pos()].AddWeight(0, 0);
		patch.stencils[pFirstEdge->next()->next()->from()->pos()].AddWeight(0, 0);
		if (patch.isQuad)
		{ patch.stencils[pFirstEdge->next()->next()->next()->from()->pos()].AddWeight(0, 0); }

		// Compute corner / edge control points
		int v = 0;
		for(HalfEdge::Face::ConstEdgeIterator it(pFirstEdge); !it.isDone(); it.advance(), v++)
		{
			const HalfEdge::Vertex * vertex = it.current()->from();
			patch.valenceArray[v] = vertex->valence() == 4 ? 0 : cos(2 * PI / vertex->valence());
			patch.iValence[v] = vertex->valence();
			if (patch.iValence[v] != 4)
			{
			   patch.stencils.bRegular = false;
			}

			const Vector3 v0 = vertex->pos();
			const float valence = float(vertex->valence());

			float total = (valence * valence);
			Vector3 pos = vertex->pos() * (valence * valence);
			patch.stencils[vertex->pos()].AddWeight(cornerIndices[primitiveOffset+v], valence * valence);

			// corner control points
			int i1 = 0, i2 = 0, j = 0;
			for(HalfEdge::Vertex::ConstEdgeIterator eit(it.current()); !eit.isDone(); eit.advance(), j++)
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(vertex->pos() == edge->from()->pos());

				const Vector3 v1 = edge->to()->pos();
				Vector3 v2 = edge->next()->to()->pos();

				// do the centroid trick for triangles.
				bool bIsTriangle = edge->next()->next()->next() == edge;
				if (bIsTriangle)
				{
					// @@ IC: why vertex->pos instead of v0?
					v2 = edge->face()->trg_centroid();
				}
				pos += v1 * 4;
				pos += v2;
				patch.stencils[v1].AddWeight(cornerIndices[primitiveOffset+v], 4);
				patch.stencils[v2].AddWeight(cornerIndices[primitiveOffset+v], 1);

				total += 5;

				// find index of "our" edge for edge control points
				if (edge == it.current()) {
					i1 = j;
				}
				if (edge == it.current()->prev()->pair()) {
					i2 = j;
				}
			}

			pos /= total;
			patch.stencils.Divide(cornerIndices[primitiveOffset+v], total);
			patch.positionArray[cornerIndices[primitiveOffset+v]] = pos;

			Vector3 edgepos1(0, 0, 0);
			Vector3 edgepos2(0, 0, 0);

			// edge control points (two per edge)
			j=0;
			for(HalfEdge::Vertex::ConstEdgeIterator eit(it.current()); !eit.isDone(); eit.advance(), j++)
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(vertex->pos() == edge->from()->pos());

				Vector3 v1 = edge->to()->pos();
				Vector3 v2 = edge->next()->to()->pos();

				// do the centroid trick for triangles.
				bool bIsTriangle = edge->next()->next()->next() == edge;
				if (bIsTriangle)
				{
					v2 = edge->face()->trg_centroid();
				}
				float costerm_a = GetCosAlpha(float(2*(j-i1)), valence);
				float costerm_b = GetCosBeta(float(2*(j-i1)-1), valence); // -1 instead of +1 b/c of edge->next()->to()
				edgepos1 += costerm_a * v1;
				edgepos1 += costerm_b * v2;
				if (!patch.stencils[v1].AddWeight(edge1Indices[primitiveOffset+v], costerm_a))
					bWrongPatch = true;
				if (!patch.stencils[v2].AddWeight(edge1Indices[primitiveOffset+v], costerm_b))
					bWrongPatch = true;
				costerm_a = GetCosAlpha(float(2*(j-i2)), valence);
				costerm_b = GetCosBeta(float(2*(j-i2)-1), valence); // -1 instead of +1 b/c of edge->next()->to()
				edgepos2 += costerm_a * v1;
				edgepos2 += costerm_b * v2;
				if (!patch.stencils[v1].AddWeight(cornerEdge2Indices[primitiveOffset+v], costerm_a))
					bWrongPatch = true;
				if (!patch.stencils[v2].AddWeight(cornerEdge2Indices[primitiveOffset+v], costerm_b))
					bWrongPatch = true;
			}

			patch.positionArray[edge1Indices[primitiveOffset+v]] = edgepos1;
			patch.positionArray[cornerEdge2Indices[primitiveOffset+v]] = edgepos2;
		}

		if (bWrongPatch || patch.stencils.stencils.count() > 32)
		{
			printf("Face %d has wrong topology - skipping\n", f);
			--p;
			continue;
		}

		Vector3 yArray[20];
		// interior control points
		v = 0;
		for(HalfEdge::Face::ConstEdgeIterator it(pFirstEdge); !it.isDone(); it.advance(), v++)
		{
			const HalfEdge::Edge * edge = it.current();

			const HalfEdge::Vertex *v0 = edge->from();
			const float costerm0 = cos( 2*PI / (float) v0->valence() );

			const HalfEdge::Vertex *v1 = it.current()->to();
			const float costerm1 = cos( 2*PI / (float) v1->valence() );

			//  p0 +------+ q0
			//     |      |
			//  e0 +======+ f0 <=== current edge
			//     |      |
			//  p1 +------+ q1

			const Vector3 e0 = edge->from()->pos();
			const Vector3 f0 = edge->to()->pos();

			const Vector3 q0 = edge->next()->to()->pos();
			Vector3 p0 = edge->prev()->from()->pos();
			Vector3 p1 = edge->pair()->next()->to()->pos();
			const Vector3 q1 = edge->pair()->prev()->from()->pos();

			const Vector3 midedgeA1 = (e0 + p0) / 2.0f;
			const Vector3 midedgeA2 = (f0 + q0) / 2.0f;
			const Vector3 midedgeB1 = (e0 + p1) / 2.0f;
			const Vector3 midedgeB2 = (f0 + q1) / 2.0f;
			patch.stencils[p0].AddWeight(interior1Indices[primitiveOffset+v], 1.0f/18);
			patch.stencils[p1].AddWeight(interior1Indices[primitiveOffset+v], -1.0f/18);
			patch.stencils[q0].AddWeight(interior2Indices[primitiveOffset+v], 1.0f/18);
			patch.stencils[q1].AddWeight(interior2Indices[primitiveOffset+v], -1.0f/18);

			bool bIsTriangleA = edge->next()->next()->next() == edge;
			if (bIsTriangleA)  // do the centroid trick for triangles.
			{
				p0 = edge->face()->trg_centroid();
			}
			const Vector3 centroidA = (p0 + e0 + f0 + q0) / 4;
			patch.stencils[e0].AddWeight(interior1Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);
			patch.stencils[f0].AddWeight(interior1Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);
			patch.stencils[q0].AddWeight(interior1Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);
			patch.stencils[p0].AddWeight(interior1Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);
			patch.stencils[e0].AddWeight(interior2Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);
			patch.stencils[f0].AddWeight(interior2Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);
			patch.stencils[q0].AddWeight(interior2Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);
			patch.stencils[p0].AddWeight(interior2Indices[primitiveOffset+v], 1.0f/4 * 2.0f/9);

			bool bIsTriangleB = edge->pair()->next()->next()->next() == edge->pair();
			if (bIsTriangleB)  // do the centroid trick for triangles.
			{
				p1 = edge->pair()->face()->trg_centroid();
			} 
			const Vector3 centroidB = (p1 + e0 + f0 + q1) / 4;
			patch.stencils[e0].AddWeight(interior1Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);
			patch.stencils[f0].AddWeight(interior1Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);
			patch.stencils[p1].AddWeight(interior1Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);
			patch.stencils[q1].AddWeight(interior1Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);
			patch.stencils[e0].AddWeight(interior2Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);
			patch.stencils[f0].AddWeight(interior2Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);
			patch.stencils[p1].AddWeight(interior2Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);
			patch.stencils[q1].AddWeight(interior2Indices[primitiveOffset+v], -1.0f/4 * 2.0f/9);

			//Vector3 y = ( 2*(p0 - p1) + (q0 - q1) ) / 18.0f;
			//instead of the above, use an expression on the vertices of subdivided faces:
			Vector3 y = (2 * (midedgeA1 - midedgeB1) + 4 * (centroidA - centroidB)) / 18.0f;

			patch.positionArray[interior1Indices[primitiveOffset+v]] = y;

			//y = ( (p0 - p1) + 2*(q0 - q1) ) / 18.0f;
			//instead of the above, use an expression on the vertices of subdivided faces:
			y = (2 * (midedgeA2 - midedgeB2) + 4 * (centroidA - centroidB)) / 18.0f;

			patch.positionArray[interior2Indices[primitiveOffset+v]] = y;
		}

		// texture coordinates
		v = 0;
		for (HalfEdge::Face::ConstEdgeIterator it(pFirstEdge); !it.isDone(); it.advance(), v++)
		{
			const HalfEdge::Edge * edge = it.current();
			const HalfEdge::Vertex * vertex = edge->from();

			// Interior texcoord.
			patch.texcoords[4 * v + 0] = vertex->tex();

			// Corner texcoord.
			uint minId = vertex->id();
			patch.texcoords[4 * v + 3] = patch.texcoords[4 * v + 0];

			for (HalfEdge::Vertex::ConstVertexIterator cit(vertex->colocals()); !cit.isDone(); cit.advance())
			{
				const HalfEdge::Vertex * colocal = cit.current();
				if (colocal->id() < minId)
				{
					minId = colocal->id();
					patch.texcoords[4 * v + 3] = colocal->tex();
				}
			}

			// Boundary texcoords.
			const HalfEdge::Edge * pair = edge->pair();
			nvDebugCheck(pair != NULL);
			nvDebugCheck(edge->id() != pair->id());

			if (edge->id() < pair->id())
			{
				patch.texcoords[4 * v + 1] = patch.texcoords[4 * v + 0];
			}
			else
			{
				patch.texcoords[4 * v + 1] = pair->to()->tex();
			}

			edge = edge->prev();
			pair = edge->pair();

			nvDebugCheck(pair != NULL);
			nvDebugCheck(edge->id() != pair->id());

			if (edge->id() < pair->id())
			{
				patch.texcoords[4 * v + 2] = patch.texcoords[4 * v + 0];
			}
			else
			{
				patch.texcoords[4 * v + 2] = pair->from()->tex();
			}
		}

		for (int iv = 2 + patch.isQuad; iv >= 0; --iv)
		{
			Vector3 v0 = patch.positionArray[patch.iEdge0(iv)];
			Vector3 v1 = patch.positionArray[patch.iEdge1(iv)];
			// v1 = v0 * cosWeight + w * sinWeight;
			patch.vW[iv] = (v1 - v0 * cos(2 * PI / patch.iValence[iv])) / sin(2 * PI / patch.iValence[iv]);
			patch.vSeed[iv] = v0;
			patch.iOwnershipMask = 0;
			patch.SetupIndices(iv, patch, iv);
		}
	}
	gregoryMesh->m_patchArray.resize(p);

	return gregoryMesh.release();
}

