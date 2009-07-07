// Copyright NVIDIA Corporation 2009 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_MESHOPTIMIZER_H
#define NV_MESH_MESHOPTIMIZER_H

#include <nvmesh/render/VertexCache.h>

#include <nvcore/Containers.h>



namespace nv
{
	// Represent arbitrary mesh connectivities.
	class MeshConnectivity
	{
	public:
		struct Vertex;
		struct Primitive;

		MeshConnectivity(uint vertexCount, const Array<uint> & indices, uint primitiveSize);

		uint primitiveSize() const { return m_primitiveSize; }
		uint primitiveCount() const { return m_primitiveArray.count(); }
		const Primitive * primitiveAt(uint p) const { return m_primitiveArray[p]; }

		uint vertexCount() const { return m_vertexArray.count(); }
		const Vertex * vertexAt(uint v) const { return m_vertexArray[v]; }

	private:

		const uint m_primitiveSize;

		Array<Vertex *> m_vertexArray;
		Array<Primitive *> m_primitiveArray;
	};

	struct MeshOptimizerState;

	struct MeshOptimizer
	{
		MeshOptimizer(uint vertexCount, const Array<uint> & indices, uint primitiveSize, uint cacheSize, VertexCache::Mode cacheMode);

		enum Method
		{
			Method_Forsyth,
			Method_Castano,
			Method_KCacheReorder,
			Method_Tipsy,
			Method_FermiBatcher,
			Method_NVTriStrip,
			Method_AMDTootle,
		};

		uint primitiveCount() const
		{
			return m_primitiveCount;
		}

		void optimize(Method method, Array<uint> * indices, Array<uint> * primitiveIndices);

		static VertexCache::Stats processIndices(const Array<uint> & indices, uint primitiveSize, uint cacheSize, VertexCache::Mode cacheMode, uint primitivePerBatch);

		static void sortIndices(Array<uint> & indices, uint primitiveSize, const Array<uint> & primitiveIndices);


	private:

		// Find face with the lowest score.
		const MeshConnectivity::Primitive * selectFirstPrimitive(Method method);
		const MeshConnectivity::Primitive * selectNextPrimitive(Method method, const MeshConnectivity::Primitive * primitive);

		uint selectFirstVertex(Method method);
		uint selectNextVertex(Method method, uint v);

		uint countNotInCacheAdjacentVertices(uint v);
		float evalKCacheScore(uint v) const;

	private:

		// Input data.
		const Array<uint> m_indices;
		const uint m_primitiveSize;
		const uint m_primitiveCount;
		
		// Cache behaviour.
		const uint m_cacheSize;
		const VertexCache::Mode m_cacheMode;
		
		// Mesh representation.
		MeshConnectivity m_mesh;

		// Current state.
		MeshOptimizerState * m_state;
	};

} // nv namespace


#endif // NV_MESH_MESHOPTIMIZER_H
