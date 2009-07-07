// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include <nvcore/BitArray.h>

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>

#include <nvmesh/MeshBuilder.h>
#include <nvmesh/MeshTopology.h>

#include "Atlas.h"

using namespace nv;

namespace 
{
	// Determine if the given mesh is a	quad mesh.
	bool isQuadMesh(const HalfEdge::Mesh * mesh)
	{
		nvDebugCheck(mesh != NULL);

		const uint faceCount = mesh->faceCount();
		for(uint i = 0; i < faceCount; i++) {
			const HalfEdge::Face * face = mesh->faceAt(i);
			if (face->edgeCount() != 4) {
				return false;
			}
		}

		return true;
	}

} // namespace


/// Ctor.
Atlas::Atlas(const HalfEdge::Mesh * mesh) : m_mesh(mesh)
{
	nvCheck(mesh != NULL);
}

// Dtor.
Atlas::~Atlas()
{
	deleteAll(m_chartArray);
}


uint Atlas::chartCount() const
{
	return m_chartArray.count();
}

const Chart * Atlas::chartAt(uint i) const
{
	return m_chartArray[i];
}
Chart * Atlas::chartAt(uint i)
{
	return m_chartArray[i];
}


/// Compute a seamless texture atlas.
bool Atlas::computeSeamlessTextureAtlas(uint w/*=1024*/, uint h/*=1024*/, Atlas::SeamlessTextureAtlas mode/*=Default*/)
{
	switch(mode)
	{
		case Tiles:
			return staTiles();
		case GroupTiles:
			return staGTiles();
		case AdaptiveTiles:
			return staATiles();
		default:
			return false;
	}
}


/// Straightforward atlas similar to what ZBrush does. This is 
/// close to the algorithm described in:
///
/// Meshed Atlases for Real-Time Procedural Solid Texturing
/// http://graphics.cs.uiuc.edu/~jch/papers/rtpst.pdf
bool Atlas::staTiles()
{
/*	if (!isQuadMesh(m_mesh)) {
		// Only handle quads for now.
		return false;
	}

	// Each face is a chart.
	const uint faceCount = m_mesh->faceCount();
	m_chartArray.resize(faceCount);

	for(uint f = 0; f < faceCount; f++) {
		m_chartArray[f].faceArray.clear();
		m_chartArray[f].faceArray.append(f);
	}
	
	// Map each face to a separate square.
	
	// Determine face layout according to width and height.
	float aspect = float(m_width) / float(m_height);

	uint i = 2;
	uint total = (m_width / (i+1)) * (m_height / (i+1));
	while(total > faceCount) {
		i *= 2;
		total = (m_width / (i+1)) * (m_height / (i+1));
	}
	
	uint tileSize = i / 2;

	int x = 0;
	int y = 0;

	m_result = new HalfEdge::Mesh();

	// Once you have that it's just matter of traversing the faces.
	for(uint f = 0; f < faceCount; f++) {
		
		// Compute texture coordinates.
		Vector2 tex[4];
		tex[0] = Vector2(float(x), float(y));
		tex[1] = Vector2(float(x+tileSize), float(y));
		tex[2] = Vector2(float(x+tileSize), float(y+tileSize));
		tex[3] = Vector2(float(x), float(y+tileSize));

		Array<uint> indexArray(4);

		const HalfEdge::Face * face = m_mesh->faceAt(f);

		int i = 0;
		for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++) {
			const HalfEdge::Edge * edge = it.current();
			const HalfEdge::Vertex * vertex = edge->from();

			HalfEdge::Vertex * newVertex = m_result->addVertex(vertex->id(), vertex->pos());

			newVertex->setTex(Vector3(tex[i], 0));
			newVertex->setNor(vertex->nor());

			indexArray.append(m_result->vertexCount() + 1);
		}

		m_result->addFace(indexArray);

		// Move to the next tile.
		x += tileSize + 1;
		if (x + tileSize > m_width) {
			x = 0;
			y += tileSize + 1;
		}
	}
	*/
	return true;
}

bool Atlas::staGTiles()
{
	if (!isQuadMesh(m_mesh)) {
		// Only handle quads for now.
		return false;
	}

	// Search for extraordinary vertices.

	return false;
}

bool Atlas::staATiles()
{
	if (!isQuadMesh(m_mesh)) {
		// Only handle quads for now.
		return false;
	}

	// Search for extraordinary vertices.

	return false;
}


// Other methods that we should experiment with:
// 
// Seamless Texture Atlases:
// http://www.cs.jhu.edu/~bpurnomo/STA/index.html
// 
// Rectangular Multi-Chart Geometry Images:
// http://graphics.cs.uiuc.edu/~jch/papers/rmcgi.pdf
// 
// Discrete differential geometry also provide a way of constructing  
// seamless quadrangulations as shown in:
// http://www.geometry.caltech.edu/pubs/TACD06.pdf
// 


void Atlas::extractCharts()
{
	const uint faceCount = m_mesh->faceCount();

	int first = 0;
	Array<uint> queue(faceCount);

	BitArray bitFlags(faceCount);
	bitFlags.clearAll();

	for (uint f = 0; f < faceCount; f++)
	{
		if (bitFlags.bitAt(f) == false)
		{
			// Start new patch. Reset queue.
			first = 0;
			queue.clear();
			queue.append(f);
			bitFlags.setBitAt(f);
			
			while (first != queue.count())
			{
				const HalfEdge::Face * face = m_mesh->faceAt(queue[first]);
				
				// Visit face neighbors of queue[first]
				for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
				{
					const HalfEdge::Edge * edge = it.current();
					nvDebugCheck(edge->pair() != NULL);
					
					if (!edge->isBoundary() && /*!edge->isSeam()*/ 
						!(edge->from()->tex() != edge->pair()->to()->tex() || edge->to()->tex() != edge->pair()->from()->tex()))
					{ 
						const HalfEdge::Face * neighborFace = edge->pair()->face();
						nvDebugCheck(neighborFace != NULL);
					
						if (bitFlags.bitAt(neighborFace->id()) == false)
						{
							queue.append(neighborFace->id());
							bitFlags.setBitAt(neighborFace->id());
						}
					}
				}
				
				first++;
			}
			
			Chart * chart = new Chart();
			chart->build(queue, m_mesh);

			m_chartArray.append(chart);
		}
	}

}


// Append all the chart meshes, and return new mesh.
HalfEdge::Mesh * Atlas::mergeCharts() const
{
	// @@ Ideally we would build the half edge mesh directly... using the mesh builder is overkill!!
	/*HalfEdge::Mesh * mesh = new HalfEdge::Mesh();

	foreach(i, m_chartArray)
	{
		Chart * chart = m_chartArray[i];
		nvCheck(chart != NULL);

		//mesh->append(chart->mesh());
	}

	//mesh->linkColocals();

	mesh->linkBoundary();

	return mesh;*/

	MeshBuilder builder;

	uint baseVertex = 0;

	foreach(i, m_chartArray)
	{
		const Chart * chart = m_chartArray[i];
		nvCheck(chart != NULL);

		const HalfEdge::Mesh * mesh = chart->mesh();

		const uint vertexCount = mesh->vertexCount();
		builder.hintVertexCount(vertexCount);

		for (uint v = 0; v < vertexCount; v++)
		{
			const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
			
			builder.addPosition(vertex->pos());
			builder.addNormal(vertex->nor());
			builder.addTexCoord(vertex->tex());
		}

		const uint faceCount = mesh->faceCount();
		builder.hintTriangleCount(faceCount);

		for (uint f = 0; f < faceCount; f++)
		{
			const HalfEdge::Face * face = mesh->faceAt(f);

			builder.beginPolygon();
			
			for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
			{
				uint v = baseVertex + it.current()->vertex()->id();
				builder.addVertex(v, v, v);
			}

			builder.endPolygon();
		}

		baseVertex += vertexCount;
	}

	builder.done();

	return builder.buildHalfEdgeMesh();
}



Chart::Chart() : m_mesh(NULL)
{
}

void Chart::build(const Array<uint> & faceArray, const HalfEdge::Mesh * originalMesh)
{
	// Copy face indices.
	m_faceArray = faceArray;

#if 0
	// Create chart mesh.
	MeshBuilder builder;
	
	const uint faceCount = faceArray.count();
	builder.hintTriangleCount(faceCount);

	for (uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = originalMesh->faceAt(faceArray[f]);

		builder.beginPolygon();
		
		for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			const HalfEdge::Vertex * vertex = it.current()->vertex();
			
			builder.addVertex(vertex->pos(), vertex->nor(), vertex->tex());
		}

		builder.endPolygon();
	}

	builder.optimize();
	builder.done();
	
	m_mesh = builder.buildHalfEdgeMesh();

	/*
	const uint chartVertexCount = m_mesh->vertexCount();

	HashMap<Vector3, uint> map(chartVertexCount);

	for (uint i = 0; i < chartVertexCount; i++)
	{
		const HalfEdge::Vertex * vertex = m_mesh->vertexAt(i);

		nvCheck(vertex->colocalCount() == 1);

		map.add(vertex->pos(), i);
	}

	m_vertexIndexArray.resize(chartVertexCount);

	for (uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = originalMesh->faceAt(faceArray[f]);

		for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			const HalfEdge::Vertex * vertex = it.current()->vertex();

			uint idx = 0;
			map.get(vertex->pos(), &idx);
			nvCheck(idx != 0);

			m_vertexIndexArray[idx] = vertex->id();
		}
	}
	*/
#else // @@ Do not use builder...

	const uint meshVertexCount = originalMesh->vertexCount();

	m_mesh = new HalfEdge::Mesh();

	Array<uint> meshIndices;
	meshIndices.resize(meshVertexCount, NIL);

	// Add vertices.
	const uint faceCount = faceArray.count();
	for (uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = originalMesh->faceAt(faceArray[f]);
		nvDebugCheck(face != NULL);

		for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			const HalfEdge::Vertex * vertex = it.current()->vertex();
			nvDebugCheck(vertex != NULL);

			if (meshIndices[vertex->id()] == NIL)
			{
				meshIndices[vertex->id()] = m_vertexIndexArray.count();
				m_vertexIndexArray.append(vertex->id());

				HalfEdge::Vertex * v = m_mesh->addVertex(Vector3(vertex->tex(), 0)); // Build mesh in 2D
				v->setNor(vertex->nor());
				v->setTex(vertex->tex());
			}
		}
	}

	m_mesh->linkColocals();

	Array<uint> faceIndices(7);

	// Add faces.
	for (uint f = 0, v = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = originalMesh->faceAt(faceArray[f]);
		nvDebugCheck(face != NULL);

		faceIndices.clear();

		for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			const HalfEdge::Vertex * vertex = it.current()->vertex();
			nvDebugCheck(vertex != NULL);

			faceIndices.append(meshIndices[vertex->id()]);
		}

		m_mesh->addFace(faceIndices);
	}

	m_mesh->linkBoundary();

	// Restore 3D positions.
	for (uint i = 0; i < m_vertexIndexArray.count(); i++)
	{
		uint idx = m_vertexIndexArray[i];
		
		const HalfEdge::Vertex * vertex = originalMesh->vertexAt(idx);
		HalfEdge::Vertex * chartVertex = m_mesh->vertexAt(i);

		chartVertex->setPos(vertex->pos());
	}

#endif
}


static void getBoundaryEdges(HalfEdge::Mesh * mesh, Array<HalfEdge::Edge *> & boundaryEdges)
{
	nvDebugCheck(mesh != NULL);

	const uint edgeCount = mesh->edgeCount();

	BitArray bitFlags(edgeCount);
	bitFlags.clearAll();

	boundaryEdges.clear();

	// Search for boundary edges. Create array of boundaries. Compute lengths and areas.
	for (uint e = 0; e < edgeCount; e++)
	{
		HalfEdge::Edge * startEdge = mesh->edgeAt(e);

		if (startEdge->isBoundary() && bitFlags.bitAt(e) == false)
		{
			nvDebugCheck(startEdge->face() != NULL);
			nvDebugCheck(startEdge->pair()->face() == NULL);

			startEdge = startEdge->pair();

			const HalfEdge::Edge * edge = startEdge;
			do {
				bitFlags.setBitAt(edge->id() / 2);
				edge = edge->next();
			} while(startEdge != edge);

			boundaryEdges.append(startEdge);
		}
	}
}


void Chart::closeHoles()
{
	Array<HalfEdge::Edge *> boundaryEdges;
	getBoundaryEdges(m_mesh.ptr(), boundaryEdges);

	const uint boundaryCount = boundaryEdges.count();
	if (boundaryCount <= 1)
	{
		// Nothing to close.
		return;
	}

	// Compute lengths and areas.
	Array<float> boundaryLengths;
	Array<float> boundaryAreas;
	Array<Vector3> boundaryCentroids;

	for (uint i = 0; i < boundaryCount; i++)
	{
		const HalfEdge::Edge * startEdge = boundaryEdges[i];
		nvCheck(startEdge->face() == NULL);

		float boundaryEdgeCount = 0;
		float boundaryLength = 0.0f;
		float boundaryArea = 0.0f;
		Vector3 boundaryCentroid(zero);

		const HalfEdge::Edge * edge = startEdge;
		do {
			Vector2 t0 = edge->from()->tex();
			Vector2 t1 = edge->to()->tex();

			boundaryEdgeCount++;
			boundaryLength += length(t1 - t0);
			boundaryArea += fabsf(t0.x() * t1.y() - t0.y() * t1.x());
			boundaryCentroid += edge->vertex()->pos();

			edge = edge->next();
		} while(edge != startEdge);

		boundaryLengths.append(boundaryLength);
		boundaryAreas.append(boundaryArea);
		boundaryCentroids.append(boundaryCentroid / boundaryEdgeCount);
	}


	// Find disk boundary.
	uint diskBoundary = 0;
	float maxLength = boundaryLengths[0];

	for (uint i = 1; i < boundaryCount; i++)
	{
		if (boundaryLengths[i] > maxLength)
		{
			maxLength = boundaryLengths[i];
			diskBoundary = i;
		}
	}


	// Close holes.
	for (uint i = 0; i < boundaryCount; i++)
	{
		if (diskBoundary == i)
		{
			// Skip disk boundary.
			continue;
		}

		HalfEdge::Edge * startEdge = boundaryEdges[i];
		nvCheck(startEdge->face() == NULL);

		// Add face and connect boundary edges.
		HalfEdge::Face * face = m_mesh->addFace();
		face->setEdge(startEdge);

		HalfEdge::Edge * edge = startEdge;
		do {
			edge->setFace(face);

			edge = edge->next();
		} while(edge != startEdge);
	}
}

bool Chart::isDisk() const
{
	nvDebugCheck(m_mesh != NULL);
	
	MeshTopology topology(m_mesh.ptr());

	// This is a requirement of the extraction/creation process.
	nvCheck(topology.connectedCount() == 1);

	return topology.isDisk();
}

