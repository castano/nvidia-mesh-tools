// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_STENCILMASK_H
#define NV_MESH_STENCILMASK_H

#include <nvcore/Containers.h>

namespace nv
{

	class StencilMask
	{
	public:
		StencilMask();
		explicit StencilMask(uint size);
		StencilMask(const StencilMask & mask);
	
		void resize(uint size);

		StencilMask & operator = (const StencilMask & mask);
		StencilMask & operator = (float value);
	
		void operator += (const StencilMask & mask);
		void operator -= (const StencilMask & mask);
		void operator *= (float scale);
		void operator /= (float scale);

		bool operator == (const StencilMask & other) const;
		bool operator != (const StencilMask & other) const;

		uint count() const { return m_weightArray.count(); }

		float operator[] (uint i) const { return m_weightArray[i]; }
		float & operator[] (uint i) { return m_weightArray[i]; }

		float sum() const;
		bool isNormalized() const;
		void normalize();

		friend void swap(StencilMask & a, StencilMask & b) 
		{
			swap(a.m_weightArray, b.m_weightArray);
		}

	private:
		
		Array<float> m_weightArray;
	
	};
	
} // nv namespace

#endif // NV_MESH_STENCILMASK_H
