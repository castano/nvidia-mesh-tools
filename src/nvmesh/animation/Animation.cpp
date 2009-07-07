// This code is in the public domain -- castanyo@yahoo.es

#include "Animation.h"
#include "Bone.h"
#include "MeshSkeleton.h"

using namespace nv;

Animation::Animation() :
	m_boneCount(0),
	m_frameRate(0),
	m_frameCount(0),
	m_sparse(false)
{
}


const Animation::Key & Animation::key(uint frame, uint bone) const
{
	nvCheck(m_sparse == false);
	return m_keyArray[frame * m_boneCount + bone];
}

Animation::Key & Animation::key(uint frame, uint bone)
{
	nvCheck(m_sparse == false);
	return m_keyArray[frame * m_boneCount + bone];
}


void Animation::allocate(uint boneCount, uint frameCount)
{
	m_boneCount = boneCount;
	m_frameCount = frameCount;
	m_keyArray.resize(keyCount());
}


// static 
Animation * Animation::importMD5(Stream * stream)
{
	// @@ Loading code should be in MD5File.
	return NULL;
}

// static
Animation * Animation::importPSA(Stream * stream)
{
	// @@ Loading code should be in UnrealFile.
	return NULL;
}




Pose::Pose(uint boneCount)
{
	m_boneArray.resize(boneCount);
}

Pose::Pose(Animation * animation)
{
	nvDebugCheck(animation != NULL);
	m_boneArray.resize(animation->boneCount());
}

Pose::Pose(Skeleton * skeleton)
{
	nvDebugCheck(skeleton != NULL);
	m_boneArray.resize(skeleton->boneCount());
}

void Pose::setAnimationFrame(Skeleton * skeleton, Animation * animation, uint frame)
{
	nvDebugCheck(skeleton != NULL);
	nvDebugCheck(animation != NULL);
	nvDebugCheck(skeleton->boneCount() == animation->boneCount());
	
	const int boneCount = (int)m_boneArray.count();
	nvDebugCheck(skeleton->boneCount() == boneCount);
	
	for (int i = 0; i < boneCount; i++)
	{
		Pose::Bone & bone = m_boneArray[i];
		const Skeleton::Bone & baseFrame = skeleton->boneAt(i);
		const Animation::Key & key = animation->key(frame, i);
		
		bone.rotation = key.rotation;
		bone.offset = key.offset;
		
		bone.localTransform = boneMatrix(bone.rotation, bone.offset);

		if (baseFrame.parent != -1)
		{
			nvCheck(baseFrame.parent < i);
			const Pose::Bone & parentBone = m_boneArray[baseFrame.parent];
			bone.objectTransform = mul(parentBone.objectTransform, bone.localTransform);
		}
		else
		{
			bone.objectTransform = bone.localTransform;
		}
		
		bone.linkTransform = mul(bone.objectTransform, baseFrame.inverse);
	}
}


void Pose::setAnimationTime(Skeleton * skeleton, Animation * animation, float time)
{
	nvDebugCheck(skeleton != NULL);
	nvDebugCheck(animation != NULL);
	nvDebugCheck(skeleton->boneCount() == animation->boneCount());

	const uint boneCount = m_boneArray.count();
	nvDebugCheck(skeleton->boneCount() == boneCount);
	
	// @@ For each bone, find previous and next key for the given time.
	
}


