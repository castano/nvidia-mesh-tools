// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_MORPHTARGET_H
#define NV_MESH_MORPHTARGET_H

#include <nvcore/Ptr.h>
#include <nvcore/RefCounted.h>

#include <nvmath/Vector.h>

#include <nvmesh/nvmesh.h>

namespace nv
{
	class Stream;
	class BaseMesh;

	NV_DECLARE_PTR(MorphTarget);


	// Morph target.
	class MorphTarget : public RefCounted
	{
		NV_FORBID_COPY(MorphTarget);
	public:

		MorphTarget();
		explicit MorphTarget(uint count);
		explicit MorphTarget(const BaseMesh & mesh);

		const uint positionCount() const;
		void setPositionCount(uint count);

		Vector3 position(uint i) const;
		void setPosition(uint i, Vector3::Arg p);

		friend Stream & operator<< (Stream & s, MorphTarget & obj);

	private:

		/// Array of vertices.
		Array<Vector3> m_positionArray;
		
	};


	// A sequence of morph targets.
	class NVMESH_CLASS MorphTargetSequence
	{
		NV_FORBID_COPY(MorphTargetSequence);
	public:

		MorphTargetSequence() : m_totalTime(0.0f)
		{
		}

		void setLinearSequence(uint frameCount, float totalTime);

		bool isLooping() const;
		void setLooping(bool);

		friend Stream & operator<< (Stream & s, MorphTargetSequence & obj);

	private:

		enum { MaxFrameTargets = 4 };

		struct Frame
		{
			float time;
			int targetCount;
			int target[MaxFrameTargets];
			float weight[MaxFrameTargets];
		};

		friend Stream & operator<< (Stream & s, Frame & frame);

	private:

		float m_totalTime;
		bool m_isLooping;
		Array<Frame> m_frameArray;

	};


	// Do we need a class that maps vertices to positions?
	// Or should we just have one position per vertex?

	class NVMESH_CLASS MorphTargetMesh
	{
		NV_FORBID_COPY(MorphTargetMesh);
	public:

		MorphTargetMesh();

		void setReferenceTarget(const BaseMesh & mesh);
		bool computeVertexMap(const BaseMesh & mesh);
		bool addMorphTarget(const BaseMesh & mesh);

		
		uint positionCount() const;
		uint vertexCount() const;
		uint morphTargetCount() const;

		uint vertexIndex(uint v) const;
		const MorphTarget * morphTarget(uint t) const;
		Vector3 position(uint t, uint v) const;

		friend Stream & operator<< (Stream & s, MorphTargetMesh & obj);

	private:

		// Map vertices to positions in the morph targets.
		Array<uint> m_vertexArray;

		MorphTargetPtr m_referenceTarget;
		Array<MorphTargetPtr> m_morphTargetArray;

	};


} // nv namespace


#endif // NV_MESH_MORPHTARGET_H
