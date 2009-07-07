// Copyright NVIDIA Corporation 2009 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmesh/render/MeshOptimizer.h>
#include <nvmesh/render/NvTriStrip.h>

#include <math.h> // powf
#include <float.h> // FLT_MAX
#include <limits.h> // INT_MAX

#include <stdio.h> // printf

#ifdef HAVE_TOOTLE
#define _LINUX
#include <tootlelib.h>
#undef _LINUX
#endif

using namespace nv;

/* Batching strategies:

  So far I've implemented a greedy approach, that considers each neighboring primitive and selects the one
  with the highest score. Score is computed based on the number of vertices that are in the cache.

  - It makes sense to weight all vertices in the batch equally.
  - However, for some weird reason weighting duplicate vertices multiple times helps?!
  - Using negative weights for vertices not in the cache does not help.
  - Considering the valence of the vertices does not seem to help.

  The first primitive is chosen based on the sum of the valence of the vertices. The goal is to start at boundaries
  or near primitives that have already been removed. That is, the less connected primitive goes first.

  The problem is that we have to optimize two metrics simultaneously. One is the primitives per batch, the other is
  the vertex per batch. We want both to be as close to the limit as possible.

  The current strategy is to have as many vertices as possible without considering the primitive limits. That works
  well for quads, where the primitive limit is generally not reached, but does not work so well for triangles and
  results in many batches that do not use all of its vertex threads.

  Once we have the batches computed, we should:
  - Sort batches to minimize overdraw.
  - Sort batches by proximity, to improve coherence on FIFO caches.
  - Sort primitives within the batch to improve coherence on FIFO caches.

*/

namespace nv
{

	// Each connectivity level
	struct MeshConnectivity::Vertex
	{
		Vertex(uint id) : id(id) {}

		void add(const MeshConnectivity::Primitive * p)
		{
			primitiveArray.append(p);
		}

		const uint id;
		Array<const MeshConnectivity::Primitive *> primitiveArray;
	};


	// Generic primitive with arbitrary adjacency information.
	struct MeshConnectivity::Primitive
	{
		Primitive(uint id, const Array<uint> & indices, uint first, uint count) : id(id)
		{
			indexArray.resize(count);
			
			for (uint i = 0; i < count; i++)
			{
				indexArray[i] = indices[first + i];
			}
		}

		const uint id;
		Array<uint> indexArray;
	};

	struct MeshOptimizerState
	{
		MeshOptimizerState(const MeshConnectivity * mesh, uint cacheSize, VertexCache::Mode mode);

		// Primitive processed flags:
		uint processedPrimitiveCount() const { return m_primitiveIndices.count(); }
		uint primitiveIndex(uint p) const { return m_primitiveIndices[p]; }
		bool isPrimitiveProcessed(uint p) const {
			return m_primitiveProcessed[p];
		}

	//	const Array<uint> & primitiveIndices() const { return m_primitiveIndices; }
	//	Array<uint> & primitiveIndices() { return m_primitiveIndices; }

	//	const Array<uint> & indices() const { return m_primitiveIndices; }
	//	Array<uint> & indices() { return m_primitiveIndices; }

		const VertexCache & cache() const { return m_cache; }

		bool isVertexInCache(uint idx) const;
		bool isVertexUsed(uint idx) const;
		bool isVertexUnused(uint idx) const;

		uint valenceSum(const MeshConnectivity::Primitive * primitive) const;
		float primitiveScore(MeshOptimizer::Method method, const MeshConnectivity::Primitive * primitive) const;

		float vertexScore_Forsyth(uint idx) const;
		float vertexScore_Castano(uint idx) const;

		uint vertexValence(uint idx) const { return m_vertexValence[idx]; }

		void addPrimitive(const MeshConnectivity::Primitive * primitive, uint idx = ~0);

		bool isVisitedStackEmpty() const { return m_vertexVisitedStack.isEmpty(); }
		uint popVisited()
		{
			uint v = m_vertexVisitedStack.back();
			m_vertexVisitedStack.pop_back();
			return v;
		}

	private:

		float evalConnectivityScore(const MeshConnectivity::Primitive * primitive);

	private:
		const MeshConnectivity * m_mesh;

		VertexCache m_cache;

		Array<bool> m_primitiveProcessed;
		Array<float> m_primitiveScore;

		Array<uint> m_vertexValence;

	public:
		// @@ LOL. Accessors are not working in gcc!
		Array<uint> m_primitiveIndices;
		Array<uint> m_indices;

		Array<uint> m_vertexVisitedStack;
	};

} // nv namespace



// Builds meshconnectivity in o(n^2)
MeshConnectivity::MeshConnectivity(uint vertexCount, const Array<uint> & indices, uint primitiveSize) :
	m_primitiveSize(primitiveSize)
{
	m_vertexArray.resize(vertexCount);

	const uint primitiveCount = indices.count() / m_primitiveSize;
	m_primitiveArray.resize(primitiveCount);

	// Build primitives.
	for (uint i = 0; i < primitiveCount; i++)
	{
		m_primitiveArray[i] = new Primitive(i, indices, i * m_primitiveSize, m_primitiveSize);
	}

	// Build vertices.
	for (uint i = 0; i < vertexCount; i++)
	{
		m_vertexArray[i] = new Vertex(i);
	}

	// Build vertex adjacency.
	for (uint i = 0; i < primitiveCount; i++)
	{
		const Primitive * primitive = m_primitiveArray[i];

		const uint indexCount = primitive->indexArray.count();
		for (uint j = 0; j < indexCount; j++)
		{
			uint idx = primitive->indexArray[j];

			m_vertexArray[idx]->add(primitive);
		}
	}
}


MeshOptimizerState::MeshOptimizerState(const MeshConnectivity * mesh, uint cacheSize, VertexCache::Mode mode) :
	m_mesh(mesh),
	m_cache(cacheSize, mode)
{
	const uint primitiveCount = m_mesh->primitiveCount();

	// Allocate space for indices of processed primitives.
	m_primitiveIndices.reserve(primitiveCount);

	// Init flags to false.
	m_primitiveProcessed.resize(primitiveCount, false);

	const uint vertexCount = m_mesh->vertexCount();

	// Init vertex visited stack.
	m_vertexVisitedStack.reserve(vertexCount);

	// Compute vertex valences.
	m_vertexValence.resize(vertexCount, 0); // Init valences to zero.

	for (uint p = 0; p < primitiveCount; p++)
	{
		const MeshConnectivity::Primitive * primitive = m_mesh->primitiveAt(p);

		const uint indexCount = primitive->indexArray.count();
		for (uint i = 0; i < indexCount; i++)
		{
			uint idx = primitive->indexArray[i];
			m_vertexValence[idx]++;
		}
	}	
}


bool MeshOptimizerState::isVertexInCache(uint idx) const
{
	nvDebugCheck(idx < m_vertexValence.count());
	return m_cache.hasVertex(idx);
}

bool MeshOptimizerState::isVertexUsed(uint idx) const
{
	nvDebugCheck(idx < m_vertexValence.count());
	return m_vertexValence[idx] > 0;
}

bool MeshOptimizerState::isVertexUnused(uint idx) const
{
	nvDebugCheck(idx < m_vertexValence.count());
	return m_vertexValence[idx] == 0;
}

uint MeshOptimizerState::valenceSum(const MeshConnectivity::Primitive * primitive) const
{
	uint sum = 0;

	const uint vertexCount = primitive->indexArray.count();
	for (uint v = 0; v < vertexCount; v++)
	{
		const uint idx = primitive->indexArray[v];

		sum += m_vertexValence[idx];
	}

	return sum;
}


float MeshOptimizerState::primitiveScore(MeshOptimizer::Method method, const MeshConnectivity::Primitive * primitive) const
{
	nvCheck (method == MeshOptimizer::Method_Forsyth || method == MeshOptimizer::Method_Castano || method == MeshOptimizer::Method_FermiBatcher);

	if (method == MeshOptimizer::Method_Forsyth)
	{
		float score = 0.0f;

		const uint vertexCount = primitive->indexArray.count();
		for (uint v = 0; v < vertexCount; v++)
		{
			const uint idx = primitive->indexArray[v];

			score += vertexScore_Forsyth(idx);
		}

		return score;
	}
	else if (method == method == MeshOptimizer::Method_Castano)
	{
		float score = 0.0f;

		const uint vertexCount = primitive->indexArray.count();
		for (uint v = 0; v < vertexCount; v++)
		{
			const uint idx = primitive->indexArray[v];

			score += vertexScore_Castano(idx);
		}

		return score;
	}
	else if (method == MeshOptimizer::Method_FermiBatcher)
	{
		uint newVertexCount = 0;
		float score = 0.0f;

		const uint vertexCount = primitive->indexArray.count();
		for (uint v = 0; v < vertexCount; v++)
		{
			const uint idx = primitive->indexArray[v];

			bool isUnique = true;
			for (uint w = 0; w < v; w++)
			{
				if (primitive->indexArray[w] == idx) {
					isUnique = false;
					break;
				}
			}
			//if (isUnique) // This causes worse results!!???
			{
				if (!m_cache.hasVertex(idx))
				{
					newVertexCount++;
				}

				if (m_cache.hasVertex(idx))
				{
					score += 1.0f;
				}
			}
		}

		// @@ WTF! this doesn't help!
		/*if (newVertexCount > m_cache.cacheSize() - m_cache.vertexCount())
		{
			// No space available for this primitive.
			return 0.0f;
		}*/

		return score;
	}

	return 0.0f;
}


// This is especially tuned for triangles. Would have to be generalized for other primitives.
float MeshOptimizerState::vertexScore_Forsyth(uint idx) const
{
	const float FindVertexScore_CacheDecayPower = 1.5f;
	const float FindVertexScore_LastTriScore = 0.75f;
	const float FindVertexScore_ValenceBoostScale = 2.0f;
	const float FindVertexScore_ValenceBoostPower = 0.5f;

	nvCheck (isVertexUsed(idx));

	float score = 0.0f;

	int cachePosition = m_cache.vertexPosition(idx);
	if (cachePosition < 0)
	{
		// Vertex is not in FIFO cache - no score.
	}
	else
	{
		if (cachePosition < 3)
		{
			// This vertex was used in the last triangle,
			// so it has a fixed score, whichever of the three
			// it's in. Otherwise, you can get very different
			// answers depending on whether you add
			// the triangle 1,2,3 or 3,1,2 - which is silly.
			score = FindVertexScore_LastTriScore;
		}
		else
		{
			const uint MaxSizeVertexCache = m_cache.cacheSize();

			nvDebugCheck (uint(cachePosition) < MaxSizeVertexCache);

			// Points for being high in the cache.
			const float scaler = 1.0f / (MaxSizeVertexCache - 3);
			score = 1.0f - (cachePosition - 3) * scaler;
			score = powf (score, FindVertexScore_CacheDecayPower);
		}
	}

	// Bonus points for having a low number of tris still to
	// use the vert, so we get rid of lone verts quickly.
	float valence = float(m_vertexValence[idx]);
	float valenceBoost = powf (valence, -FindVertexScore_ValenceBoostPower);

	score += valenceBoost * FindVertexScore_ValenceBoostScale;

	return score;
}


float MeshOptimizerState::vertexScore_Castano(uint idx) const
{
	nvCheck (isVertexUsed(idx));

	float score = 0.0f;

	// c1 = k1 * number of vertices not in the cache.
	// c2 = k2 * number of faces rendered.
	// c3 = k3 * cache position when turning black?

	// Per primitive:
	// c1 = can be computed
	// c2 = always 1? it can be higher.
	// c3 = position of the vertices and their valences.

	return score;
}


void MeshOptimizerState::addPrimitive(const MeshConnectivity::Primitive * primitive, uint firstIdx/*= ~0*/)
{
	nvCheck(!isPrimitiveProcessed(primitive->id));

	// Record list of processed primitive ids.
	m_primitiveIndices.append(primitive->id);

	// Record indices of processed primitives.
	m_indices.append(primitive->indexArray);

	const uint indexCount = primitive->indexArray.count();

	// Add primitive indices to the cache.
	m_cache.addPrimitive(primitive->indexArray.buffer(), indexCount, firstIdx);

	if (firstIdx != ~0)
	{
		if (!isVertexInCache(firstIdx))
		{
			m_vertexVisitedStack.pushBack(firstIdx);
		}
	}

	for (uint i = 0; i < indexCount; i++)
	{
		const uint idx = primitive->indexArray[i];
		if (idx != firstIdx)
		{
			if (!isVertexInCache(idx))
			{
				m_vertexVisitedStack.pushBack(idx);
			}
		}
	}

	// Update flags.
	m_primitiveProcessed[primitive->id] = true;

	//  Update vertex valence.
	for (uint i = 0; i < indexCount; i++)
	{
		uint idx = primitive->indexArray[i];

		nvCheck (m_vertexValence[idx] > 0);
		m_vertexValence[idx]--;
	}
}


MeshOptimizer::MeshOptimizer(uint vertexCount, const Array<uint> & indices, uint primitiveSize, uint cacheSize, VertexCache::Mode cacheMode) :
	m_indices(indices),	
	m_primitiveSize(primitiveSize), 
	m_primitiveCount(indices.count()/primitiveSize), 
	m_cacheSize(cacheSize),
	m_cacheMode(cacheMode),
	m_mesh(vertexCount, indices, primitiveSize),
	m_state(NULL)
{
	nvCheck(primitiveSize <= m_cacheSize);
	nvCheck(indices.count()	% primitiveSize	== 0);
}

/*static*/
VertexCache::Stats MeshOptimizer::processIndices(const Array<uint> & indices, uint primitiveSize, uint cacheSize, VertexCache::Mode cacheMode, uint primitivePerBatch)
{
	const uint primitiveCount =	indices.count() / primitiveSize;
	
	uint count = 0;

	VertexCache cache(cacheSize, cacheMode, primitivePerBatch);

	for	(uint p	= 0; p < primitiveCount; p++)
	{
		cache.addPrimitive(&indices[p * primitiveSize], primitiveSize);
	}

	cache.flush();
	
	return cache.stats();
}

/*static*/
void MeshOptimizer::sortIndices(Array<uint> & indices, uint primitiveSize, const Array<uint> & primitiveIndices)
{
	const uint indexCount = indices.count();
	const uint primitiveCount = primitiveIndices.count();
	nvCheck(indexCount / primitiveSize == primitiveCount);

	Array<uint> newIndices;
	newIndices.resize(indexCount);

	for (uint p = 0; p < primitiveCount; p++)
	{
		for (uint i = 0; i < primitiveSize; i++)
		{
			newIndices[p * primitiveSize + i] = indices[primitiveIndices[p] * primitiveSize + i];
		}
	}

	swap(newIndices, indices);
}


void MeshOptimizer::optimize(Method method, Array<uint> * indices, Array<uint> * primitiveIndices)
{
	if (method == Method_NVTriStrip)
	{
		TriStrip::SetCacheSize(m_cacheSize);
		TriStrip::SetListsOnly(true);

		uint count = m_indices.count();

		Array<uint16> shortIndices;
		shortIndices.resize(count);

		for (uint i = 0; i < count; i++)
		{
			shortIndices[i] = uint16(indices->at(i));
		}

		TriStrip::PrimitiveGroup * groupArray;
		uint16 groupCount;

		TriStrip::GenerateStrips(shortIndices.buffer(), shortIndices.count(), &groupArray, &groupCount);

		nvCheck (groupCount == 1);

		count = groupArray[0].numIndices;

		indices->clear();
		indices->reserve(count);
		for (uint i = 0; i < count; i++)
		{
			indices->append(groupArray[0].indices[i]);
		}

		delete [] groupArray;

		if (primitiveIndices != NULL)
		{
			// @@ TODO!!
		}

		return;
	}
	else if(method == Method_AMDTootle)
	{
#ifdef HAVE_TOOTLE
		nvCheck(m_primitiveSize == 3);

		Array<uint> indexArray;
		Array<uint> faceRemapArray;

		if (indices != NULL) indexArray.resize(m_indices.count());
		if (primitiveIndices != NULL) faceRemapArray.resize(m_primitiveCount);

		TootleResult result = TootleOptimizeVCache(
				m_indices.buffer(),
				m_primitiveCount,
				m_mesh.vertexCount(),
				m_cacheSize,
				(indices != NULL) ? indexArray.mutableBuffer() : NULL,
				(primitiveIndices != NULL) ? faceRemapArray.mutableBuffer() : NULL,
				TOOTLE_VCACHE_TIPSY);

		if (indices != NULL) swap(*indices, indexArray);
		if (primitiveIndices != NULL) swap(*primitiveIndices, faceRemapArray);
#endif
		return;
	}

	// Create state.
	if (method == Method_Forsyth)
	{
		// Tom Forsyth says that he uses an LRU cache, but it seems that FIFO works best.
		m_state = new MeshOptimizerState(&m_mesh, m_cacheSize, VertexCache::Mode_FIFO);
	}
	else if (method == Method_Castano)
	{
		m_state = new MeshOptimizerState(&m_mesh, m_cacheSize, VertexCache::Mode_FIFO);
	}
	else if (method == Method_Tipsy)
	{
		// I think the Tipsy paper uses a FIFO, but it appears that simulating an LRU produces best results??
		m_state = new MeshOptimizerState(&m_mesh, m_cacheSize, VertexCache::Mode_FIFO);
	}
	else if (method == Method_KCacheReorder)
	{
		m_state = new MeshOptimizerState(&m_mesh, m_cacheSize, VertexCache::Mode_FIFO);
	}
	else if (method == Method_FermiBatcher)
	{
		m_state = new MeshOptimizerState(&m_mesh, m_cacheSize, VertexCache::Mode_Batch);
	}

	// Process primitives.
	if (method == Method_Forsyth || method == Method_Castano || method == Method_FermiBatcher)
	{
		while (m_state->processedPrimitiveCount() < primitiveCount())
		{
			const MeshConnectivity::Primitive * primitive = selectFirstPrimitive(method);

			//if (primitive != NULL) printf("First Primitive %u\n", primitive->id);

			while (primitive != NULL)
			{
				// Update state.
				m_state->addPrimitive(primitive);

				primitive = selectNextPrimitive(method, primitive);

				//if (primitive != NULL) printf("Next Primitive %u\n", primitive->id);
			};
		}
	}
	else if (method == Method_Tipsy || method == Method_KCacheReorder)
	{
		while (m_state->processedPrimitiveCount() < primitiveCount())
		{
			uint v = selectFirstVertex(method);

			while (v != ~0)
			{
				// Add all primitives that use that vertex.
				const MeshConnectivity::Vertex * vertex = m_mesh.vertexAt(v);

				const uint primitiveCount = vertex->primitiveArray.count();
				for (uint i = 0; i < primitiveCount; i++)
				{
					const MeshConnectivity::Primitive * primitive = vertex->primitiveArray[i];

					if (!m_state->isPrimitiveProcessed(primitive->id))
					{
						m_state->addPrimitive(primitive, v);
					}
				}

				v = selectNextVertex(method, v);
			};
		}
	}

#if _DEBUG
	// Make sure all the valences are zero.
	for (uint i = 0; i < m_mesh.vertexCount(); i++)
	{
		nvCheck (m_state->vertexValence(i) == 0);
	}
#endif // _DEBUG

	if (indices != NULL) swap(*indices, m_state->m_indices);
	if (primitiveIndices != NULL) swap(*primitiveIndices, m_state->m_primitiveIndices);
}

// Find face with the lowest valences.
const MeshConnectivity::Primitive * MeshOptimizer::selectFirstPrimitive(Method method)
{
	nvDebugCheck (method == Method_Forsyth || method == Method_Castano || method == Method_FermiBatcher);

	if (method == Method_Forsyth || method == Method_Castano || method == Method_FermiBatcher)
	{
		const MeshConnectivity::Primitive * bestPrimitive = NULL;
		int bestValenceSum = INT_MAX;

		for (uint p = 0; p < m_primitiveCount; p++)
		{
			if (m_state->isPrimitiveProcessed(p) == false)
			{
				const MeshConnectivity::Primitive * primitive = m_mesh.primitiveAt(p);

				int valenceSum = m_state->valenceSum(primitive);
				if (valenceSum < bestValenceSum)
				{
					bestValenceSum = valenceSum;
					bestPrimitive = primitive;
				}
			}
		}

		return bestPrimitive;
	}

	return NULL;
}

const MeshConnectivity::Primitive * MeshOptimizer::selectNextPrimitive(Method method, const MeshConnectivity::Primitive * previousPrimitive)
{
	nvDebugCheck (method == Method_Forsyth || method == Method_Castano || method == Method_FermiBatcher);

	if (method == Method_Forsyth || method == Method_Castano || method == Method_FermiBatcher)
	{
		// Traverse vertices in the cache.
		const uint cacheVertexCount = m_state->cache().vertexCount();

		float bestScore = 0.0f;
		const MeshConnectivity::Primitive * bestPrimitive = NULL;

		for (uint v = 0; v < cacheVertexCount; v++)
		{
			uint idx = m_state->cache().vertexAt(v);

			// Traverse primitives adjacent to vertices in cache.
			const MeshConnectivity::Vertex * vertex = m_mesh.vertexAt(idx);

			// @@ This will visit many primitives multiple times and recompute their scores. How can we skip that? Flag/test counters?

			const uint adjacentPrimitiveCount = vertex->primitiveArray.count();
			for (uint i = 0; i < adjacentPrimitiveCount; i++)
			{
				const MeshConnectivity::Primitive * primitive = vertex->primitiveArray[i];

				if (m_state->isPrimitiveProcessed(primitive->id) == false)
				{
					const float score = m_state->primitiveScore(method, primitive);
					if (score > bestScore)
					{
						bestScore = score;
						bestPrimitive = primitive;
					}
				}
			}
		}

		return bestPrimitive;
	}

	return NULL;
}


uint MeshOptimizer::selectFirstVertex(Method method)
{
	nvDebugCheck(method == Method_Tipsy || method == Method_KCacheReorder);

	if (method == Method_Tipsy)
	{
#if 1
		// Last visited vertex that is still used.
		while (!m_state->isVisitedStackEmpty())
		{
			uint v = m_state->popVisited();

			if (m_state->isVertexUsed(v))
			{
				break;
				//return v;
			}
		}
#else
		const uint cacheVertexCount = m_state->cache().vertexCount();

		// Pick last cache vertex that is still used.
		for (uint i = 0; i < cacheVertexCount; i++)
		{
			uint v_last = m_state->cache().vertexAt(cacheVertexCount - 1 - i);
			uint v_first = m_state->cache().vertexAt(i);
			uint v = v_first;

			if (m_state->isVertexUsed(v))
			{
				return v;
			}
		}
#endif
	}

	if (method == Method_Tipsy || method == Method_KCacheReorder)
	{
		// Or the vertex with the lowest valence.
		uint bestVertex = ~0;
		uint bestValence = UINT_MAX;

		const uint vertexCount = m_mesh.vertexCount();
		for (uint v = 0; v < vertexCount; v++)
		{
			uint valence = m_state->vertexValence(v);
			if (valence > 0 && valence < bestValence)
			{
				bestValence = valence;
				bestVertex = v;
			}
		}

		return bestVertex;
	}

	return ~0;
}


uint MeshOptimizer::selectNextVertex(Method method, uint v)
{
	nvDebugCheck(method == Method_Tipsy || method == Method_KCacheReorder);

	uint bestVertex = ~0;

	if (method == Method_Tipsy)
	{
		int bestPosition = -1;
		bool bestStillInCache = false;

		// For efﬁciency reasons, we select the next fanning vertex from among this subset, which we call the 1-ring candidates.
		const MeshConnectivity::Vertex * vertex = m_mesh.vertexAt(v);

		const uint primitiveCount = vertex->primitiveArray.count();
		for (uint i = 0; i < primitiveCount; i++)
		{
			const MeshConnectivity::Primitive * primitive = vertex->primitiveArray[i];

			const uint vertexCount = primitive->indexArray.count();
			for (uint j = 0; j < vertexCount; j++)
			{
				const uint idx = primitive->indexArray[j];

				if (m_state->isVertexUsed(idx))
				{
					int p = 0;

					//const uint notInCacheVertices = countNotInCacheAdjacentVertices(idx);
					const uint notInCacheVertices = 2 * m_state->vertexValence(idx);
					const bool stillInCache = (m_state->cache().vertexPosition(idx) + notInCacheVertices <= m_state->cache().cacheSize());
					if (stillInCache) p = m_state->cache().vertexPosition(idx);

					// If there are multiple options, we pick the candidate that entered the cache the earliest.
					if (p > bestPosition)
					{
						bestStillInCache = stillInCache;
						bestPosition = p;
						bestVertex = idx;
					}


					/*const int position = m_state->cache().vertexPosition(idx);
					//const uint notInCacheVertices = countNotInCacheAdjacentVertices(idx);
					const uint notInCacheVertices = 2 * m_state->vertexValence(idx);
					const bool stillInCache = (position + notInCacheVertices <= m_state->cache().cacheSize());

					// Among them, we consider those that would still be in the cache, even after being fanned themselves.
					// @@ Apparently this doesn't help.
					if (stillInCache || !bestStillInCache)
					{
						// If there are multiple options, we pick the candidate that entered the cache the earliest.
						if (position > bestPosition)
						{
							bestStillInCache = stillInCache;
							bestPosition = position;
							bestVertex = idx;
						}
					}*/
				}
			}
		}
	}
	else if (method == Method_KCacheReorder)
	{
		const uint cacheVertexCount = m_state->cache().vertexCount();

		uint bestVertex = ~0;
		float bestScore = FLT_MAX;

		// Pick cached vertex with the best score.
		for (uint i = 0; i < cacheVertexCount; i++)
		{
			uint v = m_state->cache().vertexAt(cacheVertexCount - 1 - i);

			if (m_state->isVertexUsed(v))
			{
				float score = evalKCacheScore(v);
				if (score < bestScore)
				{
					bestScore = score;
					bestVertex = v;
				}
			}
		}

		return bestVertex;
	}

	return bestVertex;
}

uint MeshOptimizer::countNotInCacheAdjacentVertices(uint v)
{
	const MeshConnectivity::Vertex * vertex = m_mesh.vertexAt(v);

	uint count = 0;

	const uint primitiveCount = vertex->primitiveArray.count();
	for (uint i = 0; i < primitiveCount; i++)
	{
		const MeshConnectivity::Primitive * primitive = vertex->primitiveArray[i];

		if (m_state->isPrimitiveProcessed(primitive->id))
		{
			continue;
		}

		const uint vertexCount = primitive->indexArray.count();
		for (uint j = 0; j < vertexCount; j++)
		{
			uint neighborVertex = primitive->indexArray[j];

			if (!m_state->cache().hasVertex(neighborVertex))
			{
				count++;
			}
		}
	}

	return count;
}

float MeshOptimizer::evalKCacheScore(uint v) const
{
	const float k1 = 1.0f;
	const float k2 = 0.5f;
	const float k3 = 1.3f;

	float C1 = 0.0f; // Number of vertex pushes needed.
	float C2 = 0.0f; // Number of faces rendered.

	const MeshConnectivity::Vertex * vertex = m_mesh.vertexAt(v);

	VertexCache cache(m_state->cache());
	cache.resetStats();

	const uint primitiveCount = vertex->primitiveArray.count();
	for (uint p = 0; p < primitiveCount; p++)
	{
		const MeshConnectivity::Primitive * primitive = vertex->primitiveArray[p];

		if (!m_state->isPrimitiveProcessed(primitive->id))
		{
			cache.addPrimitive(primitive->indexArray.buffer(), primitive->indexArray.count(), v);

			C2 += 1;
		}
	}

	C1 = float(cache.stats().transformCount);

	// The cache position when turning black.
	float C3 = 1.0f;
	if (cache.hasVertex(v))
	{
		C3 = 1.0f - cache.distanceToExit(v);
	}

	return C1 * k1 + C2 * k2 + C3 * k3;
}

