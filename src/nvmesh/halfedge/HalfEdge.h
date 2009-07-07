// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_HALFEDGE_EDGE_H
#define NV_MESH_HALFEDGE_EDGE_H

#include <nvcore/Debug.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>

namespace nv
{

	/// Half edge edge. 
	class HalfEdge::Edge
	{
		NV_FORBID_COPY(Edge);
	public:
		
		// Default constructor.
		Edge(uint id) : m_id(id), m_next(NULL), m_prev(NULL), m_pair(NULL), m_vertex(NULL), m_face(NULL)
		{
		}
		
		uint id() const { return m_id; }
		
		// Vertex queries.
		const Vertex * vertex() const { return m_vertex; }
		Vertex * vertex() { return m_vertex; }
		void setVertex(Vertex * v) { m_vertex = v; }

		const Vertex * from() const { return m_vertex; }
		Vertex * from() { return m_vertex; }

		const Vertex * to() const { return m_next->m_vertex; }
		Vertex * to() { return m_next->m_vertex; }

		const Edge * pair() const { return m_pair; }
		Edge * pair() { return m_pair; }
		void setPair(Edge * e) { m_pair = e; }

		const Face * face() const { return m_face; }
		Face * face() { return m_face; }
		void setFace(Face * f) { m_face = f; }

		// Assumes it's a triangle.
		//Vertex * oppositeVertex() { return m_next->m_next->m_vertex; }
		
		// Edge queries.
		const Edge * next() const { return m_next; }
		Edge * next() { return m_next; }
		void setNext(Edge * e) { m_next = e; e->m_prev = this; }
		
		const Edge * prev() const { return m_prev; }
		Edge * prev() { return m_prev; }
		
		
		// Face queries. 
		bool isBoundary() const { return !(m_face && m_pair->m_face); }
		
		// @@ This is not exactly accurate, we should compare the texture coordinates...
		bool isSeam() const { return from() != pair()->to() || to() != pair()->from(); }
		
		// Geometric queries.
		Vector3 midPoint() const;
		float length() const;
		
	private:
		uint m_id;
		
		Edge * m_next;
		Edge * m_prev;	// This is not strictly half-edge, but makes algorithms easier and faster.
		Edge * m_pair;
		Vertex * m_vertex;
		Face * m_face;
		
		// no geometry data
		
	};

} // nv namespace
	

#endif // NV_MESH_HALFEDGE_EDGE_H
