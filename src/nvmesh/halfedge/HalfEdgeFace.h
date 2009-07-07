// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_HALFEDGE_FACE_H
#define NV_MESH_HALFEDGE_FACE_H

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdge.h>

namespace nv
{

	/// Face of a half-edge mesh.
	class HalfEdge::Face
	{
		NV_FORBID_COPY(Face);
	public:
		
		Face(uint id) : m_id(id), m_edge(NULL)
		{
		}
		
		uint id() const { return m_id; }
		
		Edge * edge() { return m_edge; }
		const Edge * edge() const { return m_edge; }
		void setEdge(Edge * e) { m_edge = e; }
		
		
		// The iterator that visits the edges of this face in clockwise order.
		class EdgeIterator //: public Iterator<Edge *>
		{
		public:
			EdgeIterator(Edge * e) : m_end(NULL), m_current(e) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->next();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual Edge * current() const { return m_current; }
			Vertex * vertex() const { return m_current->vertex(); }

		private:
			Edge * m_end;
			Edge * m_current;
		};

		EdgeIterator edges() { return EdgeIterator(m_edge); }
		EdgeIterator edges(Edge * edge)
		{ 
			nvDebugCheck(contains(edge));
			return EdgeIterator(edge); 
		}

		// The iterator that visits the edges of this face in clockwise order.
		class ConstEdgeIterator //: public Iterator<const Edge *>
		{
		public:
			ConstEdgeIterator(const Edge * e) : m_end(NULL), m_current(e) { }
			ConstEdgeIterator(const EdgeIterator & it) : m_end(NULL), m_current(it.current()) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->next();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual const Edge * current() const { return m_current; }
			const Vertex * vertex() const { return m_current->vertex(); }

		private:
			const Edge * m_end;
			const Edge * m_current;
		};

		ConstEdgeIterator edges() const { return ConstEdgeIterator(m_edge); }
		ConstEdgeIterator edges(const Edge * edge) const
		{ 
			nvDebugCheck(contains(edge));
			return ConstEdgeIterator(edge); 
		}

		/// Determine if this face contains the given edge.
		bool contains(const Edge * e) const
		{
			for(ConstEdgeIterator it(edges()); !it.isDone(); it.advance())
			{
				if(it.current() == e) return true;
			}
			return false;
		}

		/// Count the number of edges in this face.
		uint edgeCount() const
		{
			uint count = 0;
			for(ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) { ++count; }
			return count;
		}

		/// Determine if this is a boundary face.
		bool isBoundary() const
		{
			for(ConstEdgeIterator it(edges()); !it.isDone(); it.advance())
			{
				const Edge * edge = it.current();
				nvDebugCheck(edge->pair() != NULL);

				if (edge->pair()->face() == NULL) {
					return true;
				}
			}
			return false;
		}

		/// Count the number of boundary edges in the face.
		uint boundaryCount() const
		{
			uint count = 0;
			for(ConstEdgeIterator it(edges()); !it.isDone(); it.advance())
			{
				const Edge * edge = it.current();
				nvDebugCheck(edge->pair() != NULL);

				if (edge->pair()->face() == NULL) {
					count++;
				}
			}
			return count;
		}
		
		// Get face area.
		float area() const;
		
		// Get face normal.
		Vector3 normal() const;

		Vector3 trg_centroid() const;
		
	private:
		uint m_id;
		Edge * m_edge;
	};

} // nv namespace

#endif // NV_MESH_HALFEDGE_FACE_H
