// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_ANIMATION_H
#define NV_MESH_ANIMATION_H

#include <nvcore/Containers.h> // Array
#include <nvmath/Matrix.h>
#include <nvmath/Quaternion.h>

namespace nv
{
	class Stream;
	class Skeleton;
	

	/// Skeletal animation track.
	class Animation
	{
	public:
		
		struct Key
		{
			Quaternion rotation;
			Vector3 offset;
		};
		
		Animation();
		
		const Key & key(uint frame, uint bone) const;
		Key & key(uint frame, uint bone);
		
		uint boneCount() const;
		uint frameCount() const;
		float frameRate() const;
		
		float time() const;
		uint keyCount() const;
		
		void allocate(uint boneCount, uint frameCount);		
		
		static Animation * importMD5(Stream * stream);
		static Animation * importPSA(Stream * stream);
		
	private:
		
		uint m_boneCount;
		uint m_frameCount;
		float m_frameRate;
		
		bool m_sparse;	// = false
		
		Array<Key> m_keyArray;
		
	};

	inline uint Animation::boneCount() const 
	{
		return m_boneCount;
	}
	inline uint Animation::frameCount() const
	{
		return m_frameCount;
	}
	inline float Animation::frameRate() const
	{
		return m_frameRate;
	}
	
	inline float Animation::time() const
	{
		return m_frameRate * m_frameCount;
	}
	inline uint Animation::keyCount() const
	{
		return m_boneCount * m_frameCount;
	}
	
	
	/// Animation pose.
	class Pose
	{
	public:
		
		struct Bone
		{
			Quaternion rotation;
			Vector3 offset;
			Matrix localTransform;
			Matrix objectTransform;
			Matrix linkTransform;
		};
		
		Pose(uint boneCount);
		Pose(Animation * animation);
		Pose(Skeleton * skeleton);
		
		void setAnimationFrame(Skeleton * skeleton, Animation * animation, uint frame);
		void setAnimationTime(Skeleton * skeleton, Animation * animation, float time);
		
		const Array<Bone> & boneArray() const;
		
	private:
		
		Array<Bone> m_boneArray;
		
	};
	
	inline const Array<Pose::Bone> & Pose::boneArray() const
	{
		return m_boneArray;
	}
	
	
} // nv namespace


#endif // NV_MESH_ANIMATION_H
