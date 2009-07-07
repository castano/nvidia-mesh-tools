// This code is in the public domain -- castanyo@yahoo.es

#include "HalfEdgeMesh.h"
#include "HalfEdge.h"
#include "HalfEdgeVertex.h"
#include "HalfEdgeFace.h"

#include <nvmesh/TriMesh.h>
#include <nvmesh/QuadTriMesh.h>
#include <nvmesh/MeshBuilder.h>

using namespace nv;
using namespace HalfEdge;

Mesh::Mesh() : m_colocalVertexCount(0)
{
}

Mesh::Mesh(const Mesh * mesh)
{
	// Copy mesh vertices.
	const uint vertexCount = mesh->vertexCount();
	m_vertexArray.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		const Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex->id() == v);

		m_vertexArray[v] = new Vertex(v);
		m_vertexArray[v]->setPos(vertex->pos());
		m_vertexArray[v]->setNor(vertex->nor());
		m_vertexArray[v]->setTex(vertex->tex());
	}

	m_colocalVertexCount = vertexCount;


	// Copy mesh faces.
	const uint faceCount = mesh->faceCount();

	Array<uint> indexArray(3);

	for (uint f = 0; f < faceCount; f++)
	{
		const Face * face = mesh->faceAt(f);

		for(Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const Vertex * vertex = it.current()->from();
			indexArray.append(vertex->id());
		}

		addFace(indexArray);
		indexArray.clear();
	}
}

Mesh::~Mesh()
{
   clear();
}


void Mesh::clear()
{
	deleteAll(m_vertexArray); 
	m_vertexArray.clear();

	foreach(i, m_edgeMap)
	{
		delete m_edgeMap[i].value;
	}
	//deleteAll(m_edgeArray);	// edgeArray only contains 1/2 of the edges!
	m_edgeArray.clear();
	m_edgeMap.clear();

	deleteAll(m_faceArray);
	m_faceArray.clear();
}


Vertex * Mesh::addVertex(const Vector3 & pos)
{
	return addVertex(m_vertexArray.count(), pos);
}

Vertex * Mesh::addVertex(uint id, const Vector3 & pos)
{
	nvDebugCheck(isValid(pos));
	
	Vertex * v = new Vertex(id);
	v->setPos(pos);
	m_vertexArray.append(v);

	return v;
}

/*void Mesh::addVertices(const Mesh * mesh)
{
	nvCheck(mesh != NULL);

	// Add mesh vertices
	for (uint v = 0; v < vertexCount; v++)
	{
		const Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		Vertex * v = addVertex(vertex->pos());
		nvDebugCheck(v != NULL);

		v->setNor(vertex->nor());
		v->setTex(vertex->tex());
	}
}*/


/// Link colocal vertices.
void Mesh::linkColocals()
{
	nvDebug("--- Linking colocals:\n");
	
	const uint vertexCount = this->vertexCount();
	HashMap<Vector3, Vertex *> vertexMap(vertexCount);
	
	for (uint v = 0; v < vertexCount; v++)
	{
		Vertex * vertex = vertexAt(v);
		
		Vertex * colocal;
		if (vertexMap.get(vertex->pos(), &colocal))
		{
			colocal->linkColocal(vertex);
		}
		else
		{
			vertexMap.add(vertex->pos(), vertex);
		}
	}
	
	m_colocalVertexCount = vertexMap.count();

	nvDebug("---   %d vertex positions.\n", m_colocalVertexCount);

	// @@ Remove duplicated vertices? or just leave them as colocals?

}

Face * Mesh::addFace()
{
	Face * f = new Face(m_faceArray.count());
	m_faceArray.append(f);
	return f;
}

Face * Mesh::addFace(uint v0, uint v1, uint v2)
{
	Array<uint> indexArray(3);
	indexArray << v0 << v1 << v2;
	return addFace(indexArray, 0, 3);
}

Face * Mesh::addFace(uint v0, uint v1, uint v2, uint v3)
{
	Array<uint> indexArray(4);
	indexArray << v0 << v1 << v2 << v3;
	return addFace(indexArray, 0, 4);
}

Face * Mesh::addFace(const Array<uint> & indexArray)
{
	return addFace(indexArray, 0, indexArray.count());
}


Face * Mesh::addFace(const Array<uint> & indexArray, uint first, uint num)
{
	nvDebugCheck(first < indexArray.count());
	nvDebugCheck(num <= indexArray.count()-first);
	nvDebugCheck(num > 2);

	if (!canAddFace(indexArray, first, num)) {
		// @@ Try to add the face in the opposite winding.
		nvDebug("Warning: non manifold mesh, invalid face '%u'.\n", m_faceArray.count());
		return NULL;
	}
	
	Face * f = new Face(m_faceArray.count());
	
	Edge * firstEdge = NULL;
	Edge * last = NULL;
	Edge * current = NULL;

	for(uint i = 0; i < num-1; i++)
	{
		current = addEdge(indexArray[first+i], indexArray[first+i+1]);
		nvCheck(current != NULL);
		
		current->setFace(f);
		
		if (last != NULL) last->setNext(current);
		else firstEdge = current;
		
		last = current;
	}

	current = addEdge(indexArray[first+num-1], indexArray[first]);
	nvCheck(current != NULL);
	
	current->setFace(f);

	last->setNext(current);
	current->setNext(firstEdge);

	f->setEdge(firstEdge);
	m_faceArray.append(f);

	return f;
}

/*void Mesh::addFaces(const Mesh * mesh)
{
	nvCheck(mesh != NULL);

	Array indexArray;
	// Add faces

}*/


// Return true if the face can be added to the manifold mesh.
bool Mesh::canAddFace(const Array<uint> & indexArray, uint first, uint num) const
{
	for(uint i = 0; i < num-1; i++)
	{
		if (!canAddEdge(indexArray[first+i], indexArray[first+i+1])) {
			return false;
		}
	}

	return canAddEdge(indexArray[first+num-1], indexArray[first]);
}

// Return true if the edge doesn't exist or doesn't have any adjacent face. 
bool Mesh::canAddEdge(uint i, uint j) const
{
	if (i == j) {
		// Skip degenerate edges.
		return false;
	}

	// Same check, but taking into account colocal vertices.
	const Vertex * v0 = vertexAt(i);
	const Vertex * v1 = vertexAt(j);

	for(Vertex::ConstVertexIterator it(v0->colocals()); !it.isDone(); it.advance())
	{
		if (it.current() == v1)
		{
			// Skip degenerate edges.
			return false;
		}
	}

	// Make sure edge has not been added yet.
	Edge * edge = findEdge(i, j);
	return edge == NULL;
}

Edge * Mesh::addEdge(uint i, uint j)
{
	nvCheck(i != j);
	
	// Make sure edge is not there already.
	Edge * edge = NULL; //findEdge(i, j);
	//nvDebugCheck(edge == NULL);

	// Lookup pair.
	Edge * pair = findEdge(j, i);


	if (pair != NULL)
	{
		// Create edge with same id.
		edge = new Edge(pair->id() + 1);
		
		// Link edge pairs.
		edge->setPair(pair);
		pair->setPair(edge);
		
		// @@ I'm not sure this is necessary!
		pair->vertex()->setEdge(pair);
	}
	else
	{
		// Create edge.
		edge = new Edge(2*m_edgeArray.count());
		
		// Add only unpaired edges.
		m_edgeArray.append(edge);
	}
	
	edge->setVertex(m_vertexArray[i]);
	m_edgeMap.add(Key(i,j), edge);
	
	// Face and Next are set by addFace.
	
	return edge;
}


/// Find edge, test all colocals.
Edge * Mesh::findEdge(uint i, uint j) const
{
	Edge * edge = NULL;

	const Vertex * v0 = vertexAt(i);
	const Vertex * v1 = vertexAt(j);

	// Test all colocal pairs.
	for(Vertex::ConstVertexIterator it0(v0->colocals()); !it0.isDone(); it0.advance())
	{
		for(Vertex::ConstVertexIterator it1(v1->colocals()); !it1.isDone(); it1.advance())
		{
			Key key(it0.current()->id(), it1.current()->id());

			if (edge == NULL) {
				m_edgeMap.get(key, &edge);
#if !defined(_DEBUG)
				if (edge != NULL) return edge;
#endif
			}
			else {
				// Make sure that only one edge is found.
				nvDebugCheck(!m_edgeMap.get(key));
			}
		}
	}

	return edge;
}

/// Link boundary edges once the mesh has been created.
bool Mesh::linkBoundary()
{
	nvDebug("--- Linking boundaries:\n");
	
	int num = 0;
	
	// Create boundary edges.
	uint edgeCount = this->edgeCount();
	for(uint e = 0; e < edgeCount; e++)
	{
		Edge * edge = edgeAt(e);
		if (edge->pair() == NULL) {
			Edge * pair = new Edge(edge->id() + 1);

			uint i = edge->from()->id();
			uint j = edge->to()->id();

			Key key(j,i);
			nvCheck(!m_edgeMap.get(key));

			pair->setVertex(m_vertexArray[j]);
			m_edgeMap.add(key, pair);
			
			edge->setPair(pair);
			pair->setPair(edge);
			
			num++;
		}
	}

	// Link boundary edges.
	for (uint e = 0; e < edgeCount; e++) {
		Edge * edge = edgeAt(e);
		if (edge->pair()->face() == NULL) {
			linkBoundaryEdge(edge->pair());
		}
	}
	
	nvDebug("---   %d boundary edges.\n", num);

	// Detect boundary intersections.
	uint boundaryIntersections = 0;
	
	const uint vertexCount = this->vertexCount();
	for (uint v = 0; v < vertexCount; v++)
	{
		const Vertex * vertex = vertexAt(v);

		uint boundaryEdges = 0;
		for (Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance())
		{
			if (it.current()->isBoundary())
				boundaryEdges++;
		}

		if (boundaryEdges > 2)
		{
			nvDebugCheck((boundaryEdges & 1) == 0);
			boundaryIntersections++;
		}
	}

	if (boundaryIntersections != 0)
	{
		nvDebug("---   Invalid mesh, boundary intersections found!\n");
		return false;
	}

	return true;
}

/// Link this boundary edge.
void Mesh::linkBoundaryEdge(Edge * edge)
{
	nvCheck(edge->face() == NULL);
	
	// Make sure next pointer has not been set.
	nvCheck(edge->next() == NULL);
		
	Edge * next = edge;
	while(next->pair()->face() != NULL) {
		// Get pair prev
		Edge * e = next->pair()->next();
		while (e->next() != next->pair()) {
			e = e->next();
		}
		next = e;
	}
	edge->setNext(next->pair());
	
	// Adjust vertex edge, so that it's the boundary edge. (required for isBoundary())
	if (edge->vertex()->edge() != edge)
	{
		// Multiple boundaries in the same edge.
		//nvCheck( edge->vertex()->edge() == NULL || edge->vertex()->edge()->face() != NULL );
		edge->vertex()->setEdge(edge);
	}
}


/// Convert to tri mesh.
TriMesh * Mesh::toTriMesh() const
{
	uint triangleCount = 0;
	
	// Count triangle faces.
	const uint faceCount = this->faceCount();
	for(uint f = 0; f < faceCount; f++)
	{
		const Face * face = faceAt(f);
		triangleCount += face->edgeCount() - 2;
	}
	
	TriMesh * triMesh = new TriMesh(triangleCount, vertexCount());

	// Add vertices.
	Array<TriMesh::Vertex> & vertices = triMesh->vertices();

	const uint vertexCount = this->vertexCount();
	for(uint v = 0; v < vertexCount; v++)
	{
		const Vertex * vertex = vertexAt(v);
		
		TriMesh::Vertex triVertex;
		triVertex.id = vertices.count();
		triVertex.pos = vertex->pos();
		triVertex.nor = vertex->nor();
		triVertex.tex = vertex->tex();

		vertices.append(triVertex);
	}
	
	// Add triangles.
	Array<TriMesh::Face> & triangles = triMesh->faces();

	for(uint f = 0; f < faceCount; f++)
	{
		const Face * face = faceAt(f);

		// @@ Triangulate arbitrary polygons correctly.
		const uint v0 = face->edge()->vertex()->id();
		uint v1 = face->edge()->next()->vertex()->id();
		
		for(Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			uint v2 = it.current()->vertex()->id();

			// Skip the first two vertices.
			if (v2 == v0 || v2 == v1) continue;

			TriMesh::Face triangle;
			triangle.id = triangles.count();
			triangle.v[0] = v0;
			triangle.v[1] = v1;
			triangle.v[2] = v2;

			v1 = v2;

			triangles.append(triangle);
		}
	}
	
	return triMesh;
}

QuadTriMesh * Mesh::toQuadTriMesh() const
{
	MeshBuilder builder;
	
	const uint vertexCount = this->vertexCount();
	builder.hintVertexCount(vertexCount);

	for(uint v = 0; v < vertexCount; v++)
	{
		const Vertex * vertex = vertexAt(v);
		
		builder.addPosition(vertex->pos());
		builder.addNormal(vertex->nor());
		builder.addTexCoord(vertex->tex());
	}

	const uint faceCount = this->faceCount();
	builder.hintTriangleCount(faceCount);

	for(uint f = 0; f < faceCount; f++)
	{
		const Face * face = faceAt(f);

		builder.beginPolygon();
		
		for(Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			uint v = it.current()->vertex()->id();
			builder.addVertex(v, v, v);
		}

		builder.endPolygon();
	}

	builder.done();

	return builder.buildQuadTriMesh();
}
