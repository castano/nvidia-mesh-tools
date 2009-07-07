// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_MESHBUILDER_H
#define NV_MESH_MESHBUILDER_H

#include <nvcore/StrLib.h>
#include <nvmath/Vector.h>

namespace nv
{
	class TriMesh;
	class QuadTriMesh;
	class PolyMesh;
	namespace HalfEdge
	{
		class Mesh;
	}
	class MeshAdjacency;
	class MeshGroup;
	class MeshMaterial;
	class VertexData;
	class EdgeFeatures;


	/// Mesh builder is a helper class for importers.
	/// Ideally it should handle any vertex data, but for now it only accepts positions, 
	/// normals and texcoords.
	class MeshBuilder
	{
		NV_FORBID_COPY(MeshBuilder);
		NV_FORBID_HEAPALLOC();
	public:
		MeshBuilder();
		~MeshBuilder();

		// Builder methods.
		uint addPosition(const Vector3 & v);
		uint addNormal(const Vector3 & v);
		uint addTexCoord(const Vector2 & v);
		
		void beginGroup(uint id);
		void endGroup();
		
		uint beginMaterial(const String & name);
		void beginMaterial(uint id);
		void endMaterial();
		
		void beginPolygon();
		uint addVertex(uint v, uint n = NIL, uint t = NIL);
		uint addVertex(const Vector3 & v);
		uint addVertex(const Vector3 & v, const Vector3 & n, const Vector2 & t);
		void endPolygon();
		
		void optimize(); // eliminate duplicate components and duplicate vertices.
		
		void done();
		void reset();
		
		// Hints.
		void hintTriangleCount(uint count);
		void hintVertexCount(uint count);
		void hintPositionCount(uint count);
		void hintNormalCount(uint count);
		void hintTexCoordCount(uint count);

		// Helpers.
		void addTriangle(uint v0, uint v1, uint v2);
		void addQuad(uint v0, uint v1, uint v2, uint v3);

		// Get result.
		TriMesh * buildTriMesh(uint group = NIL, uint material = NIL) const;
		QuadTriMesh * buildQuadTriMesh(uint group = NIL, uint material = NIL) const;
		PolyMesh * buildPolyMesh(uint group = NIL, uint material = NIL) const;
		
		HalfEdge::Mesh * buildHalfEdgeMesh(uint group = NIL, uint material = NIL) const;
		MeshAdjacency * buildAdjacency(uint group = NIL, uint material = NIL) const;
		
		MeshGroup * buildMeshGroup() const;
		MeshMaterial * buildMeshMaterial() const;

		// Expose attribute indices of the unified vertex array.
		uint vertexCount() const;
		uint positionIndex(uint vertex) const;
		uint normalIndex(uint vertex) const;
		uint texcoordIndex(uint vertex) const;

	private:
		
		struct PrivateData;
		PrivateData * d;
		
	};

} // nv namespace

#endif // NV_MESH_MESHBUILDER_H
