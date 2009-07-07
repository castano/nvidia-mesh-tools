// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_TANGENTS_H
#define NV_MESH_TANGENTS_H

#include <nvcore/Containers.h>
#include <nvmath/Basis.h>

namespace nv
{
	class TriMesh;
	namespace HalfEdge { class Mesh; class Vertex; }
	
	namespace geometry
	{
		void computeMeshTangents(const TriMesh * mesh, Array<Basis> & meshBasis);
		void computeLeastSquaresMeshTangents(const TriMesh * mesh, Array<Basis> & meshBasis);

		void computeMeshTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis);
		void computeCatmullClarkTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis);

		void computeConformalMeshTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis);
		void computeConformalCatmullClarkTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis);
		
		void computeLengyelTangentBasis(const HalfEdge::Vertex * vertex, Basis * basis);
		void computeLeastSquaresTangentBasis(const HalfEdge::Vertex * vertex, Basis * basis);

	} // MeshTangents namespace
	
} // nv namespace

#endif // NV_MESH_TANGENTS_H
