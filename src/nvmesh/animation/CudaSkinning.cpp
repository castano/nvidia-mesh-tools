// Copyright NVIDIA Corporation 2007 -- Ignacio Castano <icastano@nvidia.com>

#include "MeshSkeleton.h"
#include "CudaSkinning.h"

using namespace nv;

CudaSkinning::CudaSkinning() : 
	m_blockDim(256)
{
}

void CudaSkinning::build(const Skeleton * skeleton)
{
	nvCheck(skeleton != NULL);

	nvDebug("Building CUDA skinning data:\n");

	const uint vertexCount = skeleton->vertexArray().count();

	// Allocate arrays.
	m_linkArray.clear();
	m_linkArray.reserve(2 * vertexCount);
	m_vertexArray.resize(vertexCount);
	m_blockArray.resize((vertexCount + m_blockDim - 1) / m_blockDim);

	for (uint v = 0; v < vertexCount; v++)
	{
		const uint bid = v / m_blockDim;
		const uint tid = v % m_blockDim;
		
		Block & block = m_blockArray[bid];

		if (tid == 0)
		{
			block.link_first = m_linkArray.count();
		}


		const Skeleton::Vertex & skeletonVertex = skeleton->vertexAt(v);

		bool addVertexLinks = true;

		// Search for links in this block.
		for (uint i = 0; i < tid; i++)
		{
			const Skeleton::Vertex & otherSkeletonVertex = skeleton->vertexAt(m_blockDim * bid + i);

			if (otherSkeletonVertex.link_first == skeletonVertex.link_first)
			{
				nvDebugCheck(skeletonVertex.link_num == otherSkeletonVertex.link_num);
				nvDebugCheck(skeletonVertex.link_num == m_vertexArray[m_blockDim * bid + i].link_num);

				addVertexLinks = false;
				m_vertexArray[v].link_first = m_vertexArray[m_blockDim * bid + i].link_first;
				m_vertexArray[v].link_num = skeletonVertex.link_num;
			}
		}

		if (addVertexLinks)
		{
			m_vertexArray[v].link_first = m_linkArray.count() - block.link_first;
			m_vertexArray[v].link_num = skeletonVertex.link_num;

			for (uint i = 0; i < skeletonVertex.link_num; i++)
			{
				const Skeleton::Link & skeletonLink = skeleton->linkAt(skeletonVertex.link_first + i);

				Link link;
				link.boneId = skeletonLink.bone;
				link.offset = skeletonLink.offset;
				link.weight = skeletonLink.weight;
				m_linkArray.append(link);
			}
		}
		
		if (tid == m_blockDim - 1 || v == vertexCount-1)
		{
			block.link_num = m_linkArray.count() - block.link_first;
		}
	}


	// Display info.
	nvDebug("  block dim: %d\n", this->blockDim());
	nvDebug("  link count: %d\n", this->linkCount());
	nvDebug("  vertex count: %d\n", this->vertexCount());
	nvDebug("  block count: %d\n", this->blockCount());
	nvDebug("  max link per block: %d\n", this->maxLinkPerBlock());
	nvDebug("  link per vertex ratio: %f\n", this->linkPerVertexRatio());
	nvDebug("  max link per vertex ratio: %f\n", this->maxLinkPerVertexRatio());
}



uint CudaSkinning::maxLinkPerBlock() const
{
	uint maxLinkNum = 0;

	const uint blockNum = blockCount();
	for (uint i = 0; i < blockNum; i++)
	{
		maxLinkNum = max(maxLinkNum, m_blockArray[i].link_num);
	}

	return maxLinkNum;
}

float CudaSkinning::linkPerVertexRatio() const
{
	return float(linkCount()) / float(vertexCount());
}

float CudaSkinning::maxLinkPerVertexRatio() const
{
	return float(maxLinkPerBlock()) / float(blockDim());
}
