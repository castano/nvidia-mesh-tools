// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_MESHADJACENCY_H
#define NV_MESH_MESHADJACENCY_H

#include <nvcore/Containers.h>
#include <nvmath/Vector.h>
#include <nvmesh/nvmesh.h>

namespace nv
{
class TriMesh;

/// Mesh adjacency information for triangle meshes.
/// This representation is very redundant, it uses fat edges and faces.
class TriMeshAdjacency {
public:

	/// Edge flags.
	enum EdgeFlags
	{
		EdgeFlags_None = 0x00,			///< No flags.
		EdgeFlags_Boundary = 0x01,		///< Mesh boundary (not interior).
		EdgeFlags_TextureSeam = 0x02,	///< Texture seam.
		EdgeFlags_NormalSeam = 0x04,	///< Normal seam.
		EdgeFlags_TangentSeam = 0x08	///< Tangent seam.
	};	
	
	/// Vertex with adjacency info.
	struct Vertex
	{
		uint id;		///< Vertex id.
		uint prev;		///< Prev colocal vertex.
		uint next;		///< Next colocal vertex.
	};

	/// Edge with adjacency info.
	struct Edge
	{
		uint id;		///< Edge id.
		uint v[2];		///< Vertex indices.
		uint f[2];		///< Face indices.
		uint flags;
	};

	/// Face with adjacency info.
	struct Face
	{
		uint id;		///< Face id.
		uint f[3];		///< Face indices.
		uint e[3];		///< Edge indices.
		uint v[3];		///< Vertex indices.
	};


public:

	TriMeshAdjacency(TriMesh * mesh) { buildAdjacency(mesh); }

	/// @name Vertex queries.
	//@{
	uint vertexCount() const { return m_vertexArray.count(); }
	const Vertex & vertexAt( uint v ) const { return m_vertexArray[v]; }
	Vertex & vertexAt( uint v ) { return m_vertexArray[v]; }
	const Array<Vertex> & vertexArray() const { return m_vertexArray; };
	
	// Get the vertices of the given edge of the given face.
	//void getFaceEdgeVertices(uint f, uint e, uint * v0, uint * v1);
	//@}

	/// @name Edge queries.
	//@{
	uint edgeCount() const { return m_edgeArray.count(); }
	const Edge & edge( uint e ) const { return m_edgeArray[e]; }
	Edge & edgeAt( uint e ) { return m_edgeArray[e]; }
	const Array<Edge> & edgeArray() const { return m_edgeArray; };

	// Return true if the given edge is a boundary.
	//bool isBoundaryEdge(const uint e) const;
	
	// Get the next edge of the given face.
	//uint nextEdgeOfFace(const uint f, const uint e) const;
	
	// Get the opposite face of the given edge.
	//uint oppositeFace(const uint e, const uint f) const;
	
	// Get the next boundary edge.
	//uint nextBoundaryEdge(uint e) const;
	
	/// Get texture seam flag.
	bool isTextureSeamEdge(const uint e) const { return (m_edgeArray[e].flags & EdgeFlags_TextureSeam) != 0; }
	
	/// Get normal seam flag.
	bool isNormalSeamEdge(const uint e) const { return (m_edgeArray[e].flags & EdgeFlags_NormalSeam) != 0; }
	
	/// Get tangent seam flag.
	bool isTangentSeamEdge(const uint e) const { return (m_edgeArray[e].flags & EdgeFlags_TangentSeam) != 0; }
	//@}


	/// @name Face queries.
	//@{
	uint faceCount() const { return m_faceArray.count(); }
	const Face & faceAt( uint f ) const { return m_faceArray[f]; }
	Face & faceAt( uint f ) { return m_faceArray[f]; }
	const Array<Face> & faceArray() const { return m_faceArray; };
	//@}

	
private:	

	void buildAdjacency(TriMesh * mesh);	
	
	// Proximal methods.
	void computeProximals();
	void addProximal(uint r, uint v);
	bool hasProximal(uint v) const;
	uint countProximals(uint v) const;
	uint firstProximal(uint v) const;

private:

	/// Array of vertices.
	Array<Vertex> m_vertexArray;
	
	/// Array of edges.
	Array<Edge> m_edgeArray;
	
	/// Array of faces.
	Array<Face> m_faceArray;

};

} // nv namespace

#endif // NV_MESH_MESHADJACENCY_H
