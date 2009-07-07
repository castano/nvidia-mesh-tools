// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_HALFEDGE_VERTEX_H
#define NV_MESH_HALFEDGE_VERTEX_H

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdge.h>

namespace nv
{

	/// Half edge vertex.
	class HalfEdge::Vertex
	{
		NV_FORBID_COPY(Vertex);
	public:
		
		/// Default constructor.
		Vertex(uint id) : m_id(id), m_edge(NULL), m_pos(zero), m_nor(zero), m_tex(zero)
		{
			m_next = this;
			m_prev = this;
		}
		
		
		uint id() const { return m_id; }
		
		Edge * edge() { return m_edge; }
		const Edge * edge() const { return m_edge; }
		void setEdge(Edge * e) {
			//m_edge = e;
			for(VertexIterator it(colocals()); !it.isDone(); it.advance()) { 
				it.current()->m_edge = e;
			}
		}

		Vertex * next() { return m_next; }
		const Vertex * next() const { return m_next; }
		Vertex * prev() { return m_prev; }
		const Vertex * prev() const { return m_prev; }
		
		void linkColocal(Vertex * v) {
			m_next->m_prev = v;
			v->m_next = m_next; 
			m_next = v;
			v->m_prev = this;
		}
		void unlink()
		{
			m_next->m_prev = m_prev;
			m_prev->m_next = m_next;
			m_next = this;
			m_prev = this;
		}
		
		
		/// @note This only works if linkBoundary has been called.
		bool isBoundary() const
		{
			return (m_edge && !m_edge->face());
		}

		Vector3 pos() const { return m_pos; }
		void setPos(Vector3::Arg value) { m_pos = value; }
		
		Vector3 nor() const { return m_nor; }
		void setNor(Vector3::Arg value) { m_nor = value; }

		Vector2 tex() const { return m_tex; }
		void setTex(Vector2::Arg value) { m_tex = value; }


		//	for(EdgeIterator it(iterator()); !it.isDone(); it.advance()) { ... }
		//
		//	EdgeIterator it(iterator());
		//	while(!it.isDone()) {
		//		...
		//		id.advance(); 
		//	}

		// Iterator that visits the edges around this vertex in counterclockwise order.
		class EdgeIterator //: public Iterator<Edge *>
		{
		public:
			EdgeIterator(Edge * e) : m_end(NULL), m_current(e) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->pair()->next();
				//m_current = m_current->prev()->pair();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual Edge * current() const { return m_current; }
			Vertex * vertex() const { return m_current->vertex(); }

		private:
			Edge * m_end;
			Edge * m_current;
		};

		EdgeIterator edges() { return EdgeIterator(m_edge); }
		EdgeIterator edges(Edge * edge) { return EdgeIterator(edge); }

		// Iterator that visits the edges around this vertex in counterclockwise order.
		class ConstEdgeIterator //: public Iterator<Edge *>
		{
		public:
			ConstEdgeIterator(const Edge * e) : m_end(NULL), m_current(e) { }
			ConstEdgeIterator(EdgeIterator it) : m_end(NULL), m_current(it.current()) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->pair()->next();
				//m_current = m_current->prev()->pair();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual const Edge * current() const { return m_current; }
			const Vertex * vertex() const { return m_current->to(); }

		private:
			const Edge * m_end;
			const Edge * m_current;
		};

		ConstEdgeIterator edges() const { return ConstEdgeIterator(m_edge); }
		ConstEdgeIterator edges(const Edge * edge) const { return ConstEdgeIterator(edge); }


		// Iterator that visits the edges around this vertex in counterclockwise order.
		class ReverseEdgeIterator //: public Iterator<Edge *>
		{
		public:
			ReverseEdgeIterator(Edge * e) : m_end(NULL), m_current(e) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->prev()->pair();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual Edge * current() const { return m_current; }
			Vertex * vertex() const { return m_current->vertex(); }

		private:
			Edge * m_end;
			Edge * m_current;
		};

		// Iterator that visits the edges around this vertex in counterclockwise order.
		class ReverseConstEdgeIterator //: public Iterator<Edge *>
		{
		public:
			ReverseConstEdgeIterator(const Edge * e) : m_end(NULL), m_current(e) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->prev()->pair();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual const Edge * current() const { return m_current; }
			const Vertex * vertex() const { return m_current->to(); }

		private:
			const Edge * m_end;
			const Edge * m_current;
		};



		// Iterator that visits all the colocal vertices.
		class VertexIterator //: public Iterator<Edge *>
		{
		public:
			VertexIterator(Vertex * v) : m_end(NULL), m_current(v) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->next();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual Vertex * current() const { return m_current; }

		private:
			Vertex * m_end;
			Vertex * m_current;
		};

		VertexIterator colocals() { return VertexIterator(this); }

		// Iterator that visits all the colocal vertices.
		class ConstVertexIterator //: public Iterator<Edge *>
		{
		public:
			ConstVertexIterator(const Vertex * v) : m_end(NULL), m_current(v) { }

			virtual void advance()
			{
				if (m_end == NULL) m_end = m_current;
				m_current = m_current->next();
			}

			virtual bool isDone() const { return m_end == m_current; }
			virtual const Vertex * current() const { return m_current; }

		private:
			const Vertex * m_end;
			const Vertex * m_current;
		};

		ConstVertexIterator colocals() const { return ConstVertexIterator(this); }


		uint colocalCount() const
		{
			uint count = 0;
			for(ConstVertexIterator it(colocals()); !it.isDone(); it.advance()) { ++count; }
			return count;
		}
		
		uint valence() const
		{
			uint count = 0;
			for(ConstEdgeIterator it(edges()); !it.isDone(); it.advance()) { ++count; }
			return count;
		}

		bool isFirstColocal() const
		{
			uint minId = m_id;
			for(ConstVertexIterator it(colocals()); !it.isDone(); it.advance()) { 
				minId = min(minId, it.current()->id());
			}
			return m_id == minId;
		}

	private:
		uint m_id;
		
		Edge * m_edge;
		Vertex * m_next;
		Vertex * m_prev;
		
		// geometry data
		Vector3 m_pos;
		Vector3 m_nor;
		Vector2 m_tex;
		
	};

} // nv namespace

#endif // NV_MESH_HALFEDGE_VERTEX_H
