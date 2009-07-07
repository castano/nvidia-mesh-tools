// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvcore/Stream.h>

#include <nvmath/TypeSerialization.h>

#include <nvmesh/BaseMesh.h>

#include "MorphTarget.h"

using namespace nv;

namespace nv
{

	Stream & operator<< (Stream & s, MorphTarget & obj)
	{
		return s << obj.m_positionArray;
	}

	Stream & operator<< (Stream & s, MorphTargetSequence::Frame & frame)
	{
		s << frame.time << frame.targetCount;

		for (int i = 0; i < MorphTargetSequence::MaxFrameTargets; i++)
		{
			s << frame.target[i] << frame.weight[i];
		}

		return s;
	}

	Stream & operator<< (Stream & s, MorphTargetSequence & obj)
	{
		return s << obj.m_totalTime << obj.m_isLooping << obj.m_frameArray;
	}

	Stream & operator<< (Stream & s, MorphTargetMesh & morphTargetMesh)
	{
		uint morphTargetCount = morphTargetMesh.morphTargetCount();

		// Serialize count.
		s << morphTargetCount;

		if (s.isLoading())
		{
			morphTargetMesh.m_morphTargetArray.resize(morphTargetCount);
		}

		// Serialize morph targets.
		for (uint i = 0; i < morphTargetCount; i++)
		{
			s << *morphTargetMesh.m_morphTargetArray[i];
		}

		return s;
	}

}



// Morph Target.

MorphTarget::MorphTarget()
{
}

MorphTarget::MorphTarget(uint count)
{
	m_positionArray.resize(count);
}

MorphTarget::MorphTarget(const BaseMesh & mesh)
{
	// Assume no colocal vertices.
	const uint vertexCount = mesh.vertexCount();

	m_positionArray.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		m_positionArray[v] = mesh.vertexAt(v).pos;
	}
}

const uint MorphTarget::positionCount() const
{
	return m_positionArray.count();
}

void MorphTarget::setPositionCount(uint count)
{
	m_positionArray.resize(count);
}

Vector3 MorphTarget::position(uint i) const
{
	return m_positionArray[i];
}
void MorphTarget::setPosition(uint i, Vector3::Arg p)
{
	m_positionArray[i] = p;
}




// Morph Target Sequence.

// Morph from one target to the next sequentially.
void MorphTargetSequence::setLinearSequence(const uint frameCount, const float totalTime)
{
	nvDebugCheck(totalTime > 0.0f);

	m_totalTime = totalTime; 
	m_frameArray.resize(frameCount);

	for (uint i = 0; i < frameCount; i++)
	{
		m_frameArray[i].time = (i * totalTime) / frameCount;
		m_frameArray[i].targetCount = 1;
		
		m_frameArray[i].target[0] = i;
		m_frameArray[i].weight[0] = 1.0f;

		for (uint t = 1; t < MaxFrameTargets; t++)
		{
			m_frameArray[i].target[t] = -1;
			m_frameArray[i].weight[t] = 0.0f;
		}
	}
}


bool MorphTargetSequence::isLooping() const
{
	return m_isLooping;
}

void MorphTargetSequence::setLooping(bool b)
{
	m_isLooping = b;
}



// Morph Target Meshs.

MorphTargetMesh::MorphTargetMesh()
{
}

void MorphTargetMesh::setReferenceTarget(const BaseMesh & mesh)
{
	m_referenceTarget = new MorphTarget(mesh);
}

bool MorphTargetMesh::computeVertexMap(const BaseMesh & mesh)
{
	const uint positionCount = this->positionCount();

	HashMap<Vector3, uint> map(positionCount);

	for (uint i = 0; i < positionCount; i++)
	{
		map.add(m_referenceTarget->position(i), i);
	}

	
	const uint vertexCount = mesh.vertexCount();
	m_vertexArray.resize(vertexCount);

	for (uint i = 0; i < vertexCount; i++)
	{
		const Vector3 & pos = mesh.vertexAt(i).pos;

		uint idx = 0;
		if (!map.get(pos, &idx))
		{
			// This morph target cannot be used with the given mesh.
			return false;
		}

		m_vertexArray[i] = idx;
	}

	return true;
}

bool MorphTargetMesh::addMorphTarget(const BaseMesh & mesh)
{
	nvDebugCheck(m_referenceTarget != NULL);

	MorphTargetPtr target( new MorphTarget(mesh) );

	if (target->positionCount() != m_referenceTarget->positionCount())
	{
		return false;
	}
	
	m_morphTargetArray.append(target);

	return true;
}


uint MorphTargetMesh::positionCount() const
{
	nvDebugCheck(m_referenceTarget != NULL);
	return m_referenceTarget->positionCount();
}

uint MorphTargetMesh::vertexCount() const
{
	return m_vertexArray.count();
}

uint MorphTargetMesh::morphTargetCount() const
{
	return m_morphTargetArray.count();
}

uint MorphTargetMesh::vertexIndex(uint v) const
{
	return m_vertexArray[v];
}

const MorphTarget * MorphTargetMesh::morphTarget(uint t) const
{
	return m_morphTargetArray[t].ptr();
}

Vector3 MorphTargetMesh::position(uint t, uint v) const
{
	const MorphTarget * target = this->morphTarget(t);
	nvDebugCheck(target != NULL);

	const uint idx = this->vertexIndex(v);
	nvDebugCheck(idx < target->positionCount());

	return target->position(idx);
}
