// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_VERTEXDATA_H
#define NV_MESH_VERTEXDATA_H

#include <nvcore/Containers.h>
#include <nvmath/Vector.h>

namespace nv
{

/*
Vertex data should support multiple input channels.
Vertex data has information about colocal or proximal vertices.
Welding and snapping operate on VertexData.
*/



/// Vertex data.
class VertexData
{
public:
	struct Channel;
	struct ChannelId;
	struct Vertex;

	enum Semantic {
		Semantic_Position,
		Semantic_Normal,
		Semantic_TexCoord,
		Semantic_Tangent,
		Semantic_BiTangent
	};

	enum FormatHint {
		FormatHint_Float,
		FormatHint_Float2,
		FormatHint_Float3,
		FormatHint_Float4,
		FormatHint_UByte4,
		FormatHint_Short2,
		FormatHint_Short4,
		FormatHint_Dec3N,
		FormatHint_UDec3N,
	};
	
	VertexData() : m_vertexCount(0) {}
	VertexData(uint vertexNum) : m_vertexCount(vertexNum) {}
	
	// Vertex channel methods.
	Channel * createChannel(Semantic semantic, uint set);
	void destroyChannel(Channel * channel);
	
	const Channel * findChannel(Semantic semantic, uint set = 0) const;
	Channel * findChannel(Semantic semantic, uint set = 0);
	
	const Channel * positions(uint set = 0) const { return findChannel(Semantic_Position, set); }
	const Channel * normals(uint set = 0) const { return findChannel(Semantic_Normal, set); }
	const Channel * texcoords(uint set = 0) const { return findChannel(Semantic_TexCoord, set); }


	// Unify channels under a single index array.
	void unifyChannels();
	void weldVertices();
	void unweldVertices();

	bool isVertexEqual(const uint v0, const uint v1) const;
	uint vertexHash(uint v) const;
	
	const uint vertexCount() const { return m_vertexCount; }

private:
	
	uint m_vertexCount;
	HashMap<ChannelId, Channel *> m_channelMap;

	//Array<Vertex> m_vertexArray;

};


struct VertexData::ChannelId
{
	ChannelId() {}
	ChannelId(Semantic s, uint i) : semantic(s), set(i) {}
	ChannelId(const ChannelId & c) : semantic(c.semantic), set(c.set) {}
	void operator=(const ChannelId & c) { semantic = c.semantic; set = c.set; }
	
	bool operator==(const ChannelId & id) const { return semantic == id.semantic && set == id.set; }
	bool operator!=(const ChannelId & id) const { return semantic != id.semantic || set != id.set; }
	
	Semantic semantic;
	uint set;
};

/// Vertex channel.
struct VertexData::Channel
{
	Channel() {}
	Channel(Semantic semantic, uint set) : id(semantic, set) {}
	Channel(const ChannelId & id) : id(id) {}
	Channel(Semantic semantic, uint set, uint capacity) : id(semantic, set) { data.resize(capacity); }
	Channel(const ChannelId & id, uint capacity) : id(id) { data.resize(capacity); }

	ChannelId id;
	FormatHint formatHint;
	Array<Vector4> data;

	bool unified;
	Array<uint> indexArray;
};

/// Vertex link.
struct VertexData::Vertex
{
	uint id;
//	Link proximal;
};

/*
bool operator==(const VertexData::Vertex & v0, const VertexData::Vertex & v1)
{
	const int count = v0.channelIndices.count();
	nvCheck(v1.channelIndices.count() == count);

	for(int i = 0; i < count; i++) {
		if (v0.channelIndices[i] != v0.channelIndices[i]) {
			return false;
		}
	}
	
	return true;
}
*/

} // nv namespace

#endif // NV_MESH_VERTEXDATA_H
