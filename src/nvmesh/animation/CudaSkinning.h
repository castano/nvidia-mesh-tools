// Copyright NVIDIA Corporation 2007 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_CUDASKINNING_H
#define NV_MESH_CUDASKINNING_H

#include <nvmesh/nvmesh.h>

namespace nv
{
	class Skeleton;

	/// CUDA skinning data and setup code.
	class CudaSkinning
	{
	public:

		struct Link
		{
			uint boneId;
			Vector3 offset;
			float weight;
		};

		struct Vertex
		{
			uint link_first;
			uint link_num;
		};

		struct Block
		{
			uint link_first;
			uint link_num;
		};

		CudaSkinning();

		void setBlockDim(uint dim);
		uint blockDim() const;

		uint linkCount() const;
		uint vertexCount() const;
		uint blockCount() const;

		const Link & linkAt(uint i) const;
		const Vertex & vertexAt(uint i) const;
		const Block & blockAt(uint i) const;


		uint maxLinkPerBlock() const;
		float linkPerVertexRatio() const;
		float maxLinkPerVertexRatio() const;

		void build(const Skeleton * skeleton);


	private:

		uint m_blockDim;

		Array<Link> m_linkArray;
		Array<Vertex> m_vertexArray;
		Array<Block> m_blockArray;

	};


	inline void CudaSkinning::setBlockDim(uint dim)
	{
		m_blockDim = dim;
	}

	inline uint CudaSkinning::blockDim() const
	{
		return m_blockDim;
	}

	inline uint CudaSkinning::linkCount() const
	{
		return m_linkArray.count();
	}

	inline uint CudaSkinning::vertexCount() const
	{
		return m_vertexArray.count();
	}

	inline uint CudaSkinning::blockCount() const
	{
		return m_blockArray.count();
	}

	inline const CudaSkinning::Link & CudaSkinning::linkAt(uint i) const
	{
		return m_linkArray[i];
	}

	inline const CudaSkinning::Vertex & CudaSkinning::vertexAt(uint i) const
	{
		return m_vertexArray[i];
	}

	inline const CudaSkinning::Block & CudaSkinning::blockAt(uint i) const
	{
		return m_blockArray[i];
	}


} // nv namespace


#endif // NV_MESH_CUDASKINNING_H
