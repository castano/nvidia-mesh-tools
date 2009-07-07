// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmath/nvmath.h>

#include "StencilMask.h"

using namespace nv;

StencilMask::StencilMask()
{
}

StencilMask::StencilMask(uint count)
{
	// Initialize weights to zero.
	m_weightArray.resize(count, 0U);
}

StencilMask::StencilMask(const StencilMask & mask)
{
	nvDebugCheck(mask.count() == count());
	m_weightArray = mask.m_weightArray;
}

void StencilMask::resize(uint size)
{
	m_weightArray.resize(size, 0U);
}

StencilMask & StencilMask::operator = (const StencilMask & mask)
{
	nvDebugCheck(mask.count() == count());
	m_weightArray = mask.m_weightArray;
	return *this;
}

StencilMask & StencilMask::operator = (float value)
{
	foreach(i, m_weightArray)
	{
		m_weightArray[i] = value;
	}
	return *this;
}

void StencilMask::operator += (const StencilMask & mask)
{
	nvDebugCheck(mask.count() == count());
	
	const uint count = m_weightArray.count();
	for (uint i = 0; i < count; i++)
	{
		m_weightArray[i] += mask.m_weightArray[i];
	}
}

void StencilMask::operator -= (const StencilMask & mask)
{
	nvDebugCheck(mask.count() == count());
	
	const uint count = m_weightArray.count();
	for (uint i = 0; i < count; i++)
	{
		m_weightArray[i] -= mask.m_weightArray[i];
	}
}

void StencilMask::operator *= (float scale)
{
	const uint count = m_weightArray.count();
	for (uint i = 0; i < count; i++)
	{
		m_weightArray[i] *= scale;
	}
}

void StencilMask::operator /= (float scale)
{
	float rcp = 1.0f / scale;
	*this *= rcp;
}

bool StencilMask::operator == (const StencilMask & other) const
{
	const uint count = m_weightArray.count();

	if (other.count() != count) return false;

	for (uint i = 0; i < count; i++)
	{
		//if (m_weightArray[i] != other.m_weightArray[i])
		if (!equal(m_weightArray[i], other.m_weightArray[i], 0.01f))
			return false;
	}

	return true;
}

bool StencilMask::operator != (const StencilMask & other) const
{
	return !(*this == other);
}

float StencilMask::sum() const
{
	float total = 0.0f;

	const uint count = m_weightArray.count();
	for (uint i = 0; i < count; i++)
	{
		total += m_weightArray[i];
	}
	
	return total;
}

bool StencilMask::isNormalized() const
{
	return equal(sum(), 1.0f);
}

void StencilMask::normalize()
{
	*this /= sum();
}

