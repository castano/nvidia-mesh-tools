// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_NORMALS_H
#define NV_MESH_NORMALS_H

#include <nvmath/Vector.h>

namespace nv
{
	namespace HalfEdge
	{
		class Mesh;
		class Vertex;
		class Edge;
	}

	class MeshSmoothGroup;
	class BaseMesh;
	class TriMesh;
	class QuadTriMesh;

	enum NormalFlags
	{
		WeightFaceArea = 0x01,
		WeightFaceAngle = 0x02
	};

	namespace MeshNormals
	{
		bool hasNormals(const BaseMesh * mesh);

		void computeNormals(TriMesh * mesh, uint flags);
		void computeNormals(QuadTriMesh * mesh, uint flags);
	//	void computeNormals(HalfEdge::Mesh * mesh, NormalFlags flags, MeshSmoothGroup * meshSmoothGroup = NULL);
	//	void computeNormals(TriMesh * mesh, NormalFlags flags, MeshSmoothGroup * meshSmoothGroup = NULL);
	
		// @@ Not implemented:
		//void computeCreasedNormals(HalfEdge::Mesh * mesh, bool weightFaceArea, bool weightFaceAngle, MeshSmoothGroup * meshSmoothGroup = NULL);

		// @@ These functions should be moved to subdiv/LimitSurface.h

		void computeLoopNormals(HalfEdge::Mesh * mesh);
		void computeCatmullClarkNormals(HalfEdge::Mesh * mesh);
		
		Vector3 computeCatmullClarkNormal(const HalfEdge::Vertex * vertex);
		Vector3 computeCatmullClarkTangent(const HalfEdge::Edge * edge);

	} // MeshNormals namespace
	
} // nv namespace

#endif // NV_MESH_NORMALS_H
