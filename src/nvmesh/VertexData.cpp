// This code is in the public domain -- castanyo@yahoo.es

#include <nvmesh/VertexData.h>

using namespace nv;

/// Create a new vertex channel. Does not check for duplicates. 
VertexData::Channel * VertexData::createChannel(Semantic semantic, uint set)
{
	ChannelId id(semantic, set);
	Channel * channel = new Channel(id, m_vertexCount);
	m_channelMap.add(id, channel);
	return channel;
}

void VertexData::destroyChannel(Channel * channel)
{
	nvDebugCheck(channel != NULL);
	m_channelMap.remove(channel->id);
	delete channel;
}

const VertexData::Channel * VertexData::findChannel(Semantic semantic, uint set /*= 0*/) const
{
	Channel * channel = NULL;
	m_channelMap.get(ChannelId(semantic, set), &channel);
	return channel;
}

VertexData::Channel * VertexData::findChannel(Semantic semantic, uint set /*= 0*/)
{
	Channel * channel = NULL;
	m_channelMap.get(ChannelId(semantic, set), &channel);
	return channel;
}


// Unify channels under a single index array.
void VertexData::unifyChannels()
{
	// Array of vertex hashes
	// Paste code from weld routine.

}

/// Weld vertices that have the same values.
void VertexData::weldVertices()
{
}

/// Unweld vertices.
void VertexData::unweldVertices()
{
}


/// Compare vertices, comparing all the channels.
bool VertexData::isVertexEqual(const uint v0, const uint v1) const
{
	foreach(c, m_channelMap) {
		const Channel * channel = m_channelMap[c].value;
		nvDebugCheck(channel != NULL);

		if (channel->indexArray[v0] != channel->indexArray[v1]) {
			return false;
		}
	}

	return true;
}


/// Compute the hash for all the vertex channels.
uint VertexData::vertexHash(uint v) const
{
	uint h = 0;
	
	// We should only use one of the channel (just position 0)
	foreach(c, m_channelMap) {
		const Channel * channel = m_channelMap[c].value;
		nvDebugCheck(channel != NULL);
		
		//h += hash<Vector4>()(channel->data[v]);
		h += hash<uint>()(channel->indexArray[v]);
	}
	
	return h;
}


