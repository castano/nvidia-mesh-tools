// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_POLYMESH_H
#define NV_MESH_POLYMESH_H

#include <nvcore/Containers.h>
#include <nvmath/Vector.h>
#include <nvmesh/nvmesh.h>
#include <nvmesh/BaseMesh.h>

namespace nv
{
class TriMesh;
class QuadTriMesh;

	
/// Polygon mesh.
class PolyMesh : public BaseMesh
{
public:
	struct Face;
	typedef BaseMesh::Vertex Vertex;

	PolyMesh(uint indexCount, uint faceCount, uint vertexCount) :
		BaseMesh(vertexCount),
		m_indexArray(indexCount), m_faceArray(faceCount) {}

	// Face methods.
	uint faceCount() const { return m_faceArray.count(); }
	const Face & faceAt(uint i) const { return m_faceArray[i]; }
	const Array<Face> & faces() const { return m_faceArray; }
	
	// Conversion.
	TriMesh * toTriMesh(bool triangulate);
	QuadTriMesh * toQuadTriMesh();
	
	friend Stream & operator<< (Stream & s, PolyMesh & obj);

private:
	
	Array<uint> m_indexArray;
	Array<Face> m_faceArray;
	
};


/// PolyMesh face.
struct PolyMesh::Face
{
	uint id;
	uint first, num;
};

} // nv namespace

#endif // NV_MESH_POLYMESH_H
