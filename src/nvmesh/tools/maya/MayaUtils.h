// Copyright NVIDIA Corporation 2007 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_MAYAUTILS_H
#define NV_MESH_MAYAUTILS_H

#include <nvmesh/nvmesh.h>

#include <maya/MTime.h>
#include <maya/MAnimControl.h>

class MStatus;
class MDagPath;
class MFnMesh;
class MFnSkinCluster;

namespace nv
{
	namespace HalfEdge { class Mesh; }
	class Skeleton;
	

	class MayaTime
	{
	public:
		MayaTime() : m_initialTime(currentTime()) {}
		
		MayaTime(MTime t) : m_initialTime(currentTime())
		{
			setCurrentTime(t);
		}
		~MayaTime()
		{
			setCurrentTime(m_initialTime);
		}

		void setCurrentTime(MTime t)
		{
			MAnimControl::setCurrentTime(t);
		}

		MTime currentTime() const
		{
			return MAnimControl::currentTime();
		}

		const MTime m_initialTime;
	};


	namespace MayaUtils
	{
		// Build a half edge mesh from the given maya mesh.
		HalfEdge::Mesh * buildHalfEdgeMesh(const MDagPath & node);

		bool getMeshPositions(const MDagPath & dagPath, Array<Vector3> * points);

		// Build a skeleton from the given maya mesh.
		Skeleton * buildSkeleton(const MDagPath & node);

		MStatus getSkinAndMesh(/*ref*/ MFnMesh & meshFn, /*out*/MFnSkinCluster & skinCluster);
		MStatus getSkinAndMeshGMAN(/*ref*/ MFnMesh & meshFn, /*out*/MFnSkinCluster & skinCluster);

	} // MayaUtils namespace
}


#endif // NV_MESH_MAYAUTILS_H
