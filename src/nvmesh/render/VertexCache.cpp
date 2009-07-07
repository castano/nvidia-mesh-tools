// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include "VertexCache.h"

#include <nvcore/Debug.h>
#include <string.h> // memcpy


using namespace nv;

const uint WARP_SIZE = 32;


VertexCache::VertexCache(uint cacheSize/*= 32*/, Mode mode/*= Mode_Batch*/, uint maxPrimitivePerBatch/*= 0*/) :
	m_cacheSize(cacheSize), m_mode(mode), m_maxPrimitivePerBatch(maxPrimitivePerBatch)
{
	m_head = 0;
	m_vertexCount = 0;
	m_indices = new uint[m_cacheSize];

	resetStats();
}

VertexCache::VertexCache(const VertexCache & other) : 
	m_cacheSize(other.m_cacheSize), 
	m_mode(other.m_mode),
	m_maxPrimitivePerBatch(m_maxPrimitivePerBatch)
{
	m_head = other.m_head;
	m_vertexCount = other.m_vertexCount;

	m_indices = new uint[m_cacheSize];
	memcpy(m_indices, other.m_indices, sizeof(m_vertexCount * sizeof(uint)));

	m_stats = other.m_stats;
}

VertexCache::~VertexCache()
{
	delete [] m_indices;
}
	
void VertexCache::operator=(const VertexCache & other)
{
	nvCheck(m_cacheSize == other.m_cacheSize);
	nvCheck(m_mode == other.m_mode);
	nvCheck(m_maxPrimitivePerBatch == other.m_maxPrimitivePerBatch);

	m_head = other.m_head;
	m_vertexCount = other.m_vertexCount;

	m_indices = new uint[m_cacheSize];
	memcpy(m_indices, other.m_indices, sizeof(m_vertexCount * sizeof(uint)));

	m_stats = other.m_stats;
}

void VertexCache::clear()
{
	m_head = 0;
	m_vertexCount = 0;
	m_primitiveCount = 0;
}

uint VertexCache::vertexCount() const
{
	return m_vertexCount;
}

uint VertexCache::vertexAt(uint idx) const
{
	nvCheck(idx < m_vertexCount);
	return m_indices[(m_head + idx) % m_cacheSize];
}


// Return vertex position relative to the head.
int VertexCache::vertexPosition(uint v) const
{
	for (uint i = 0; i < m_vertexCount; i++)
	{
		const int idx = (m_head + i) % m_cacheSize;
		if (m_indices[idx] == v) {
			return i;
		}
	}
	return -1;
}

// Return normalized distance to exit. 0 = last, 1 = first.
float VertexCache::distanceToExit(uint v) const
{
	nvDebugCheck(hasVertex(v));
	return float(m_cacheSize - 1 - vertexPosition(v)) / float(m_cacheSize - 1);
}

bool VertexCache::hasVertex(uint v) const
{
	return vertexPosition(v) != -1;
}
	
bool VertexCache::addVertex(uint v)
{
	nvDebugCheck(m_vertexCount <= m_cacheSize);

	if (m_mode == Mode_Batch)
	{
		if (hasVertex(v)) return true;
		if (m_vertexCount == m_cacheSize) return false;

		m_indices[m_vertexCount++] = v;
	}
	else
	{
		int position = vertexPosition(v);
		if (position >= 0)
		{
			if (m_mode == Mode_LRU && position != 0)
			{
				// Shift elements right.
				for (int i = position; i > 0; i--)
				{
					const int fm = (m_head + i - 1 + m_cacheSize) % m_cacheSize;
					const int to = (m_head + i) % m_cacheSize;
					m_indices[to] = m_indices[fm];
				}

				// Insert vertex at the head.
				m_indices[m_head] = v;
			}
		}
		else
		{
			// Add vertex to the cache:
			if (m_vertexCount < m_cacheSize)
			{
				m_indices[m_vertexCount++] = v;
			}
			else
			{
				m_indices[m_head] = v;

				m_head = (m_head + 1) % m_cacheSize;
			}
		}
	}

	return true;
}

void VertexCache::addPrimitive(const uint * indices, uint primitiveSize, uint firstIdx/*= ~0*/)
{
	m_stats.primitiveCount++;

	uint newVertexCount = 0;
	for (uint i = 0; i < primitiveSize;	i++)
	{
		if (!hasVertex(indices[i])) newVertexCount++;
	}

	if (m_mode == Mode_Batch)
	{
		m_primitiveCount++;

		bool fullBatch = false;
		if (m_primitiveCount > m_maxPrimitivePerBatch)
		{
			fullBatch = true;
			m_stats.fullBatchCount++;
		}

		bool newBatch = false;
		if (m_vertexCount + newVertexCount > m_cacheSize || fullBatch)
		{
			newBatch = true;
			m_stats.batchCount++;
			m_stats.transformCount += m_vertexCount;
			m_vertexCount = 0;
			m_primitiveCount = 0;
		}
		
		for (uint i = 0; i < primitiveSize; i++)
		{
			bool result	= addVertex(indices[i]);
			nvCheck(result == true);
		}
	}
	else
	{
		uint evictedCount = 0;
		if (m_vertexCount + newVertexCount > m_cacheSize) 
		{
			evictedCount = m_vertexCount + newVertexCount - m_cacheSize;
		}
		nvCheck(evictedCount <= newVertexCount);
		
		m_stats.transformCount += evictedCount;

		if (firstIdx != ~0) addVertex(firstIdx);
		for (uint i = 0; i < primitiveSize; i++)
		{
			if (indices[i] != firstIdx) addVertex(indices[i]);
		}
	}
}

void VertexCache::flush()
{
	m_stats.transformCount += m_vertexCount;

	// This is only for batch mode:
	if (m_primitiveCount > 0)
	{
		nvCheck(m_vertexCount > 0);
		m_stats.batchCount++;
	}

	clear();
}


void VertexCache::resetStats()
{
	m_stats.primitiveCount = 0;
	m_stats.transformCount = 0;
	m_stats.batchCount = 0;
	m_stats.fullBatchCount = 0;
}
