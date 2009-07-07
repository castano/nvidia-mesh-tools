// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include "AccMesh.h"
#include "AccPatch.h"
#include "RemapFaces.h"

#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>

#include <nvmesh/MeshBuilder.h>

#include <nvmesh/render/MeshOptimizer.h>

#include <nvcore/Radix.h>
#include <nvcore/Ptr.h>

#include <limits.h> // UINT_MAX

using namespace nv;


//-------------------------------
//  VertexTopology
//-------------------------------

VertexTopology::VertexTopology()
{
}

VertexTopology::VertexTopology(const VertexTopology & vt)
{
	memcpy(this, &vt, sizeof(VertexTopology));
}

void VertexTopology::operator=(const VertexTopology & vt)
{
	memcpy(this, &vt, sizeof(VertexTopology));
}

bool VertexTopology::operator==(const VertexTopology & vt) const
{
	return memcmp (this, &vt, sizeof(VertexTopology)) == 0;
}

VertexTopology::VertexTopology(const HalfEdge::Vertex * vertex, const HalfEdge::Edge * firstEdge)
{
	// This should be prevented in advance.
	nvCheck (vertex->valence() <= 16);

	// Clear all fields.
	memset(this, 0, sizeof(VertexTopology));

	valence = vertex->valence();
	isBoundary = vertex->isBoundary();

    // Regular vertices must have valence 4 and be surrounded by quad faces.
	if (valence == 4 && !isBoundary)
	{
		bool allQuads = true;
		for (HalfEdge::Vertex::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance())
		{
			const HalfEdge::Edge * edge = it.current();
			const HalfEdge::Face * face = edge->face();
			if (face != NULL)
			{
				if (face->edgeCount() != 4) {
					allQuads = false;
					break;
				}
			}
		}

		isRegular = allQuads;
	}

	// Traverse edges, compute one ring flags.
	for (HalfEdge::Vertex::ConstEdgeIterator it(firstEdge); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();
		const HalfEdge::Face * face = edge->face();

		faceBitField <<= 2;
		if (face == NULL) faceBitField |= 3;
		else if (face->edgeCount() == 3) faceBitField |= 1;
	}
}


//-------------------------------
//  AccMesh
//-------------------------------

AccMesh::AccMesh()
{
}

AccMesh::~AccMesh()
{
	foreach (i, m_faceTopologyArray)
	{
	//	delete m_faceTopologyArray[i]->bezierAccStencil;
	//	delete m_faceTopologyArray[i]->gregoryAccStencil;
		delete m_faceTopologyArray[i];
	}
	foreach (i, m_regularPatchArray)
	{
		delete m_regularPatchArray[i].faceOneRing;
        delete m_regularPatchArray[i].bezierAccPatch;
        delete m_regularPatchArray[i].gregoryAccPatch;
        delete m_regularPatchArray[i].pmAccPatch;
		delete m_regularPatchArray[i].texCoordPatch;
	}
	foreach (i, m_quadPatchArray)
	{
		delete m_quadPatchArray[i].faceOneRing;
        delete m_quadPatchArray[i].bezierAccPatch;
        delete m_quadPatchArray[i].gregoryAccPatch;
        delete m_quadPatchArray[i].pmAccPatch;
		delete m_quadPatchArray[i].texCoordPatch;
	}
	foreach (i, m_triPatchArray)
	{
		delete m_triPatchArray[i].faceOneRing;
        delete m_triPatchArray[i].gregoryAccPatch;
        delete m_triPatchArray[i].pmAccPatch;
		delete m_triPatchArray[i].texCoordPatch;
	}
}

uint AccMesh::evaluateTopologyId(const HalfEdge::Face * face, const HalfEdge::Edge * edge)
{
	uint topologyId = 0;
	
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges(edge)); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();

		const uint vertexTopologyId = this->addVertexTopology(edge->vertex(), edge);
		nvCheck(vertexTopologyId < 256);

		topologyId <<= 8;
		topologyId |= vertexTopologyId;
	}

	return topologyId;
}

FaceTopology * AccMesh::addFaceTopology(const HalfEdge::Face * face, const HalfEdge::Edge ** firstEdge)
{
	AutoPtr<FaceTopology> faceTopology(new FaceTopology);
	
	faceTopology->bezierAccStencil = NULL;
	faceTopology->gregoryAccStencil = NULL;

	faceTopology->edgeCount = face->edgeCount();
	faceTopology->isRegular = false;

	// Evaluate face topology id for each edge, pick the smallest and record the edge that produced it.
	uint topologyId = UINT_MAX;
	*firstEdge = NULL;

	for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();
	
		uint t = evaluateTopologyId(face, edge);
		if (t < topologyId)
		{
			topologyId = t;
			*firstEdge = edge;
		}
	}

	faceTopology->topologyId = topologyId;

	if ((topologyId & 0xFF) == ((topologyId >> 8) & 0xFF) && 
		(topologyId & 0xFF) == ((topologyId >> 16) & 0xFF) && 
		(topologyId & 0xFF) == ((topologyId >> 24) & 0xFF))
	{
		uint idx = (topologyId & 0xFF);
		faceTopology->isRegular = m_vertexTopologyArray[idx].isRegular;
	}

	// @@ edgeCount and topologyId could be const members.

	// Add face topology.
	uint index = 0;
	if (!m_faceTopologyMap.get(topologyId, &index))
	{
		index = m_faceTopologyArray.count();
		faceTopology->index = index;
		m_faceTopologyArray.append(faceTopology.release());
		m_faceTopologyMap.add(topologyId, index);
	}

	return m_faceTopologyArray[index];
}

uint AccMesh::faceTopologyCount() const
{
	return m_faceTopologyArray.count();
}

const FaceTopology * AccMesh::faceTopologyAt(uint i) const
{
	return m_faceTopologyArray[i];
}

uint AccMesh::addVertexTopology(const HalfEdge::Vertex * vertex, const HalfEdge::Edge * firstEdge)
{
	VertexTopology vertexTopology(vertex, firstEdge);

	uint index = 0;
	if (m_vertexTopologyMap.get(vertexTopology, &index))
	{
		nvCheck(memcmp (&m_vertexTopologyArray[index], &vertexTopology, sizeof(VertexTopology)) == 0);
	}
	else
	{
		index = m_vertexTopologyArray.count();
		m_vertexTopologyArray.append(vertexTopology);
		m_vertexTopologyMap.add(vertexTopology, index);
	}

	return index;
}

uint AccMesh::vertexTopologyCount() const
{
	return m_vertexTopologyArray.count();
}

uint AccMesh::addRegularPatch(FaceOneRing * faceOneRing)
{
    uint idx = m_regularPatchArray.count();
    m_regularPatchArray.resize(idx + 1);
	m_regularPatchArray[idx].faceType = FaceType_Regular;
	m_regularPatchArray[idx].faceOneRing = faceOneRing;
    return idx;
}

uint AccMesh::addQuadPatch(FaceOneRing * faceOneRing)
{
    uint idx = m_quadPatchArray.count();
    m_quadPatchArray.resize(idx + 1);
	m_quadPatchArray[idx].faceType = FaceType_Quad;
	m_quadPatchArray[idx].faceOneRing = faceOneRing;
    return idx;
}

uint AccMesh::addTriPatch(FaceOneRing * faceOneRing)
{
    uint idx = m_triPatchArray.count();
    m_triPatchArray.resize(idx + 1);
	m_triPatchArray[idx].faceType = FaceType_Tri;
	m_triPatchArray[idx].faceOneRing = faceOneRing;
    return idx;
}


void AccMesh::setRegularPatch(uint p, AccPatch * patch)
{
	nvDebugCheck(patch != NULL);

    PatchType patchType = patch->patchType();

    if (patchType == PatchType_Bezier) m_regularPatchArray[p].bezierAccPatch = static_cast<BezierAccPatch *>(patch);
    else if (patchType == PatchType_Gregory) m_regularPatchArray[p].gregoryAccPatch = static_cast<GregoryAccPatch *>(patch);
    else if (patchType == PatchType_PmRegular) m_regularPatchArray[p].pmAccPatch = static_cast<SurfacePatch *>(patch);
	else if (patchType == PatchType_TexCoord) m_regularPatchArray[p].texCoordPatch = static_cast<TexCoordPatch *>(patch);
	else if (patchType == PatchType_TangentSpace) m_regularPatchArray[p].tangentSpacePatch = static_cast<TangentSpacePatch *>(patch);
}

void AccMesh::setQuadPatch(uint p, AccPatch * patch)
{
	nvDebugCheck(patch != NULL);

    PatchType patchType = patch->patchType();

    if (patchType == PatchType_Bezier) m_quadPatchArray[p].bezierAccPatch = static_cast<BezierAccPatch *>(patch);
    else if (patchType == PatchType_Gregory) m_quadPatchArray[p].gregoryAccPatch = static_cast<GregoryAccPatch *>(patch);
    else if (patchType == PatchType_PmQuad) m_quadPatchArray[p].pmAccPatch = static_cast<SurfacePatch *>(patch);
	else if (patchType == PatchType_TexCoord) m_quadPatchArray[p].texCoordPatch = static_cast<TexCoordPatch *>(patch);
	else if (patchType == PatchType_TangentSpace) m_quadPatchArray[p].tangentSpacePatch = static_cast<TangentSpacePatch *>(patch);
}

void AccMesh::setTriPatch(uint p, AccPatch * patch)
{
	nvDebugCheck(patch != NULL);

    PatchType patchType = patch->patchType();

    if (patchType == PatchType_Gregory) m_triPatchArray[p].gregoryAccPatch = static_cast<GregoryAccPatch *>(patch);
    else if (patchType == PatchType_PmTriangle) m_triPatchArray[p].pmAccPatch = static_cast<SurfacePatch *>(patch);
	else if (patchType == PatchType_TexCoord) m_triPatchArray[p].texCoordPatch = static_cast<TexCoordPatch *>(patch);
	else if (patchType == PatchType_TangentSpace) m_triPatchArray[p].tangentSpacePatch = static_cast<TangentSpacePatch *>(patch);
}

uint AccMesh::patchCount() const
{
    return regularPatchCount() + quadPatchCount() + triPatchCount();
}

const AccMesh::Patch & AccMesh::patchAt(uint i) const
{
    if (i < regularPatchCount()) return regularPatchAt(i);

    i -= regularPatchCount();
    if (i < quadPatchCount()) return quadPatchAt(i);

    i -= quadPatchCount();
    nvCheck (i < triPatchCount());
    
    return triPatchAt(i);
}


void AccMesh::sort(const Sort * order, uint orderCount)
{
	nvDebugCheck (order != NULL);

	sort(order, orderCount, m_regularPatchArray);
	sort(order, orderCount, m_quadPatchArray);
	sort(order, orderCount, m_triPatchArray);
}


void AccMesh::sort(const Sort * order, uint orderCount, Array<Patch> & patchArray)
{
	nvDebugCheck (order != NULL);

	const uint patchCount = patchArray.count();
	if (patchCount == 0)
	{
		// Nothing to do.
		return;
	}

	RadixSort radix;

	for (uint i = 0; i < orderCount; i++)
	{
		Array<uint> patchValue;
		patchValue.resize(patchCount);

		if (order[i] == Sort_ByVertexCount)
		{
			for (uint p = 0; p < patchCount; p++)
			{
				patchValue[p] = patchArray[p].faceOneRing->vertexCount();
			}
		}
		else if (order[i] == Sort_ByTopology)
		{
			for (uint p = 0; p < patchCount; p++)
			{
				patchValue[p] = patchArray[p].faceOneRing->faceTopology()->index;
			}
		}
		else if (order[i] == Sort_ByFaceId)
		{
			for (uint p = 0; p < patchCount; p++)
			{
				patchValue[p] = patchArray[p].faceOneRing->face()->id();
			}
		}

		radix.sort(patchValue);

		uint32 * ranks = radix.ranks();

		Array<Patch> tmpPatchArray;
		tmpPatchArray.resize(patchArray.size());

		for (uint p = 0; p < patchCount; p++)
		{
			tmpPatchArray[p] = patchArray[ranks[p]];
		}

		swap(tmpPatchArray, patchArray);
	}
}


enum Granularity
{
	Granularity_Single,
	Granularity_PerTopology,
	Granularity_PerVertexCount,
};

void AccMesh::optimize(Granularity granularity/*= Granularity_Single*/)
{
	optimize(granularity, m_regularPatchArray);
	optimize(granularity, m_quadPatchArray);
	optimize(granularity, m_triPatchArray);
}

void AccMesh::optimize(Granularity granularity/*= Granularity_Single*/, Array<Patch> & patchArray)
{
	const uint patchCount = patchArray.count();
	if (patchCount == 0)
	{
		// Nothing to do.
		return;
	}

	Array<Patch> tmp(patchCount);

	uint maxVertexCount = 0;
	for	(uint p = 0; p < patchCount; p++)
	{
		maxVertexCount = max(maxVertexCount, patchArray[p].faceOneRing->vertexCount());
	}
	nvCheck (maxVertexCount <= 32);

	uint firstPatch = 0;
	Array<uint>	indexArray(maxVertexCount * patchCount);

	for	(uint p = 0; p < patchCount;)
	{
		uint currentVertexCount = patchArray[firstPatch].faceOneRing->vertexCount();
		uint currentTopologyIdx = patchArray[firstPatch].faceOneRing->faceTopology()->index;

		uint minVertexIndex = ~0;
		uint maxVertexIndex = 0;

		for	(; p < patchCount; p++)
		{
			const FaceOneRing * faceOneRing = patchArray[p].faceOneRing;

			const uint vertexCount = faceOneRing->vertexCount();
			nvDebugCheck (vertexCount <= maxVertexCount);

			const uint topologyIdx = faceOneRing->faceTopology()->index;
			nvDebugCheck (topologyIdx < faceTopologyCount());

			if (((granularity == Granularity_PerVertexCount) && currentVertexCount != vertexCount) ||
				((granularity == Granularity_PerTopology) && currentTopologyIdx != topologyIdx))
			{
				break;
			}

			// Update index array to optimize and compute vertex range.
			uint v = 0, id;
			for (v = 0; v < vertexCount; v++)
			{
				id = faceOneRing->vertexAt(v)->id();

				indexArray.append(id);

				minVertexIndex = min(minVertexIndex, id);
				maxVertexIndex = max(maxVertexIndex, id);
			}
			for (; v < maxVertexCount; v++)
			{
				indexArray.append(id);
			}
		}

		foreach(i, indexArray)
		{
			indexArray[i] -= minVertexIndex;
		}

		MeshOptimizer batcher(maxVertexIndex - minVertexIndex + 1, indexArray, maxVertexCount, 32, VertexCache::Mode_Batch);

		Array<uint> patchIndices;
		batcher.optimize(MeshOptimizer::Method_FermiBatcher, NULL, &patchIndices);

		uint count = p - firstPatch;
		tmp.resize(count);

		for (uint i = 0; i < count; i++)
		{
			tmp[i] = patchArray[firstPatch + patchIndices[i]];
		}
		for (uint i = 0; i < count; i++)
		{
			patchArray[firstPatch + i] = tmp[i];
		}

		firstPatch = p;
		indexArray.clear();
	}

	// @@ Debug print stats?
}


HalfEdge::Mesh * AccMesh::tessellate(uint level)
{
    MeshBuilder builder;

    Array<int> positionIndices;
    positionIndices.resize((level + 1) * (level + 1));
    
    Array<int> normalIndices;
    normalIndices.resize((level + 1) * (level + 1));

    const uint patchCount = this->patchCount();
	for	(uint p = 0; p < patchCount; p++)
	{
        const Patch & patch = patchAt(p);

#pragma message(NV_FILE_LINE "Add support for triangle tessellation")
        if (patch.faceType == AccMesh::FaceType_Tri) 
            continue;

        // Add vertices.
	    for	(uint y = 0; y <= level; y++)
	    {
            const float v = float(y) / level;

	        for	(uint x = 0; x <= level; x++)
	        {
                const float u = float(x) / level;

                Vector3 pos, du, dv;
                if (patch.bezierAccPatch != NULL) patch.bezierAccPatch->evaluateSurface(u, v, &pos, &du, &dv);
                if (patch.gregoryAccPatch != NULL) patch.gregoryAccPatch->evaluateSurface(u, v, &pos, &du, &dv);
                if (patch.pmAccPatch != NULL) patch.pmAccPatch->evaluateSurface(u, v, &pos, &du, &dv);

                positionIndices[y * (level + 1) + x] = builder.addPosition(pos);

                Vector3 normal = normalizeSafe(cross(du, dv), Vector3(zero), 0);
                
                normalIndices[y * (level + 1) + x] = builder.addNormal(normal);
            }
        }

        // Add triangles.
	    for	(uint y = 0; y < level; y++)
	    {
	        for	(uint x = 0; x < level; x++)
	        {
                int p0 = positionIndices[(y + 0) * (level + 1) + x + 0];
                int p1 = positionIndices[(y + 0) * (level + 1) + x + 1];
                int p2 = positionIndices[(y + 1) * (level + 1) + x + 0];
                int p3 = positionIndices[(y + 1) * (level + 1) + x + 1];

                int n0 = normalIndices[(y + 0) * (level + 1) + x + 0];
                int n1 = normalIndices[(y + 0) * (level + 1) + x + 1];
                int n2 = normalIndices[(y + 1) * (level + 1) + x + 0];
                int n3 = normalIndices[(y + 1) * (level + 1) + x + 1];

                builder.beginPolygon();
                builder.addVertex(p0, n0);
                builder.addVertex(p1, n1);
                builder.addVertex(p2, n2);
                builder.endPolygon();

                builder.beginPolygon();
                builder.addVertex(p2, n2);
                builder.addVertex(p1, n1);
                builder.addVertex(p3, n3);
                builder.endPolygon();
            }
        }
    }

    builder.done();

    //builder.optimize();

    return builder.buildHalfEdgeMesh();
}