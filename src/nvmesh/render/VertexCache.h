// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_VERTEXCACHE_H
#define NV_MESH_VERTEXCACHE_H

#include <nvmesh/nvmesh.h>

namespace nv
{
	struct VertexCache
	{
		enum Mode
		{
			Mode_Batch,
			Mode_FIFO,
			Mode_LRU,
		};

		VertexCache(uint cacheSize = 32, Mode mode = Mode_Batch, uint maxPrimitivePerBatch = 32);
		VertexCache(const VertexCache & other);
		~VertexCache();
		
		void operator=(const VertexCache & other);

		uint cacheSize() const { return m_cacheSize; }
		Mode cacheMode() const { return m_mode; }

		void clear();

		uint vertexCount() const;
		uint vertexAt(uint idx) const;

		int vertexPosition(uint v) const;
		float distanceToExit(uint v) const;

		bool hasVertex(uint	v) const;
		bool addVertex(uint	v);

		void addPrimitive(const uint * indices, uint primitiveSize, uint firstIndex = ~0);
		void flush();

		struct Stats
		{
			uint primitiveCount;
			uint transformCount;
			uint batchCount;
			uint fullBatchCount;
		};

		// Get stats
		void resetStats();
		const Stats & stats() { return m_stats; }

	private:
		
		// Invariants:
		const uint m_cacheSize;
		const Mode m_mode;
		const uint m_maxPrimitivePerBatch;

		// Cache state:
		uint m_head;
		uint m_vertexCount;
		uint m_primitiveCount;
		uint * m_indices;

		// Statistics:
		Stats m_stats;

	};

} // nv namespace

#endif // NV_MESH_VERTEXCACHE_H
