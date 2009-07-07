// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_HALFEDGE_H
#define NV_MESH_HALFEDGE_H

#include <nvcore/Containers.h>
#include <nvmath/Vector.h>
#include <nvmesh/nvmesh.h>

namespace nv
{
	class TriMesh;
	class QuadTriMesh;

	namespace HalfEdge
	{
		class Edge;
		class Face;
		class Vertex;

		/// Simple half edge mesh designed for dynamic mesh manipulation.
		/// Our half edge mesh has support for colocal vertices, but it's not possible to
		/// determine what vertices are being used by which faces.
		class Mesh
		{
		public:

			Mesh();
			Mesh(const Mesh * mesh);
			~Mesh();

			void clear();

			Vertex * addVertex(const Vector3 & pos);
			Vertex * addVertex(uint id, const Vector3 & pos);
			//void addVertices(const Mesh * mesh);
			
			void linkColocals();

			Face * addFace();
			Face * addFace(uint v0, uint v1, uint v2);
			Face * addFace(uint v0, uint v1, uint v2, uint v3);
			Face * addFace(const Array<uint> & indexArray);
			Face * addFace(const Array<uint> & indexArray, uint first, uint num);
			//void addFaces(const Mesh * mesh);

			bool linkBoundary();


			// Vertices
			uint vertexCount() const { return m_vertexArray.count(); }
			const Vertex * vertexAt(int i) const { return m_vertexArray[i]; }
			Vertex * vertexAt(int i) { return m_vertexArray[i]; }

			uint colocalVertexCount() const { return m_colocalVertexCount; }

			// Faces
			uint faceCount() const { return m_faceArray.count(); }
			const Face * faceAt(int i) const { return m_faceArray[i]; }
			Face * faceAt(int i) { return m_faceArray[i]; }
			
			// Edges
			uint edgeCount() const { return m_edgeArray.count();  }
			const Edge * edgeAt(int i) const { return m_edgeArray[i]; }
			Edge * edgeAt(int i) { return m_edgeArray[i]; }
			
            class ConstVertexIterator;

		    class VertexIterator
		    {
                friend class ConstVertexIterator;
		    public:
			    VertexIterator(Mesh * mesh) : m_mesh(mesh), m_current(0) { }

			    virtual void advance() { m_current++; }
			    virtual bool isDone() const { return m_current == m_mesh->vertexCount(); }
                virtual Vertex * current() const { return m_mesh->vertexAt(m_current); }

		    private:
                HalfEdge::Mesh * m_mesh;
			    uint m_current;
		    };
		    VertexIterator vertices() { return VertexIterator(this); }

		    class ConstVertexIterator
		    {
		    public:
			    ConstVertexIterator(const Mesh * mesh) : m_mesh(mesh), m_current(0) { }
                ConstVertexIterator(class VertexIterator & it) : m_mesh(it.m_mesh), m_current(it.m_current) { }

			    virtual void advance() { m_current++; }
			    virtual bool isDone() const { return m_current == m_mesh->vertexCount(); }
                virtual const Vertex * current() const { return m_mesh->vertexAt(m_current); }

		    private:
                const HalfEdge::Mesh * m_mesh;
			    uint m_current;
		    };
		    ConstVertexIterator vertices() const { return ConstVertexIterator(this); }

            class ConstFaceIterator;

		    class FaceIterator
		    {
                friend class ConstFaceIterator;
		    public:
			    FaceIterator(Mesh * mesh) : m_mesh(mesh), m_current(0) { }

			    virtual void advance() { m_current++; }
			    virtual bool isDone() const { return m_current == m_mesh->faceCount(); }
                virtual Face * current() const { return m_mesh->faceAt(m_current); }

		    private:
                HalfEdge::Mesh * m_mesh;
			    uint m_current;
		    };
		    FaceIterator faces() { return FaceIterator(this); }

		    class ConstFaceIterator
		    {
		    public:
			    ConstFaceIterator(const Mesh * mesh) : m_mesh(mesh), m_current(0) { }
                ConstFaceIterator(const FaceIterator & it) : m_mesh(it.m_mesh), m_current(it.m_current) { }

			    virtual void advance() { m_current++; }
			    virtual bool isDone() const { return m_current == m_mesh->faceCount(); }
                virtual const Face * current() const { return m_mesh->faceAt(m_current); }

		    private:
                const HalfEdge::Mesh * m_mesh;
			    uint m_current;
		    };
		    ConstFaceIterator faces() const { return ConstFaceIterator(this); }

            class ConstEdgeIterator;

		    class EdgeIterator
		    {
                friend class ConstEdgeIterator;
		    public:
			    EdgeIterator(Mesh * mesh) : m_mesh(mesh), m_current(0) { }

			    virtual void advance() { m_current++; }
			    virtual bool isDone() const { return m_current == m_mesh->edgeCount(); }
                virtual Edge * current() const { return m_mesh->edgeAt(m_current); }

		    private:
                HalfEdge::Mesh * m_mesh;
			    uint m_current;
		    };
		    EdgeIterator edges() { return EdgeIterator(this); }

		    class ConstEdgeIterator
		    {
		    public:
			    ConstEdgeIterator(const Mesh * mesh) : m_mesh(mesh), m_current(0) { }
                ConstEdgeIterator(const EdgeIterator & it) : m_mesh(it.m_mesh), m_current(it.m_current) { }

			    virtual void advance() { m_current++; }
			    virtual bool isDone() const { return m_current == m_mesh->edgeCount(); }
                virtual const Edge * current() const { return m_mesh->edgeAt(m_current); }

		    private:
                const HalfEdge::Mesh * m_mesh;
			    uint m_current;
		    };
		    ConstEdgeIterator edges() const { return ConstEdgeIterator(this); }

            // @@ Add half-edge iterator.



			// Convert to tri mesh.
			TriMesh * toTriMesh() const;
			QuadTriMesh * toQuadTriMesh() const;
			
		private:
			
			bool canAddFace(const Array<uint> & indexArray, uint first, uint num) const;
			bool canAddEdge(uint i, uint j) const;
			Edge * addEdge(uint i, uint j);

			Edge * findEdge(uint i, uint j) const;

			void linkBoundaryEdge(Edge * edge);
			
		private:
			
			Array<Vertex *> m_vertexArray;
			Array<Edge *> m_edgeArray;
			Array<Face *> m_faceArray;

			struct Key {
				Key() {}
				Key(const Key & k) : p0(k.p0), p1(k.p1) {}
				Key(uint v0, uint v1) : p0(v0), p1(v1) {}
				void operator=(const Key & k) { p0 = k.p0; p1 = k.p1; }
				bool operator==(const Key & k) const { return p0 == k.p0 && p1 == k.p1; }

				uint p0;
				uint p1;
			};
			friend struct hash<Mesh::Key>;

			HashMap<Key, Edge *> m_edgeMap;
			
			uint m_colocalVertexCount;

		};
		/*
		// This is a much better hash than the default and greatly improves performance!
		template <> struct hash<Mesh::Key>
		{
			uint operator()(const Mesh::Key & k) const { return k.p0 + k.p1; }
		};
		*/

	} // HalfEdge namespace

} // nv namespace

#endif // NV_MESH_HALFEDGE_H
