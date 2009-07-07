// This code is in the public domain -- castanyo@yahoo.es

#include <nvcore/Stream.h>

#include "UnrealFile.h"

using namespace nv;
using namespace unreal;

namespace nv
{
	// @@ Move to Vector.cpp
	static Stream & operator << (Stream & s, Vector3 & v)
	{
		if (s.isLoading())
		{
			float x, y, z;
			s << x << y << z;
			v = Vector3(x, y, z);
			return s;
		}
		else
		{
			float x = v.x(), y = v.y(), z = v.z();
			return s << x << y << z;
		}
	}

	// @@ Move to Quaternion.cpp
	static Stream & operator << (Stream & s, Quaternion & q)
	{
		if (s.isLoading())
		{
			float x, y, z, w;
			s << x << y << z << w;
			q = Quaternion(x, y, z, w);
			return s;
		}
		else
		{
			float x = q.x(), y = q.y(), z = q.z(), w = q.w();
			return s << x << y << z << w;
		}
	}
}


Stream & nv::unreal::operator << (Stream & s, JointPos & jp)
{
	return s << jp.orientation << jp.position << jp.length << jp.xSize << jp.ySize << jp.zSize;
}

Stream & nv::unreal::operator << (Stream & s, AnimInfoBinary & aib)
{
	s.serialize(aib.name, 64);
	s.serialize(aib.group, 64);
	
	s << aib.totalBones;

	s << aib.rootInclude << aib.keyCompressionStyle << aib.keyQuotum << aib.keyReduction;
	s << aib.trackTime << aib.animRate << aib.startBone << aib.firstRawFrame << aib.numRawFrames;
	
	return s;
}

Stream & nv::unreal::operator << (Stream & s, ChunkHeader & ch)
{
	s.serialize(ch.chunkID, 20);
	return s << ch.typeFlag << ch.dataSize << ch.dataCount;

}

Stream & nv::unreal::operator << (Stream & s, Material & m)
{
	s.serialize(m.materialName, 64);
	return s << m.textureIndex << m.polyFlags << m.auxMaterial << m.auxFlags << m.lodBias << m.lodStyle;
}

Stream & nv::unreal::operator << (Stream & s, Bone & b)
{
	s.serialize(b.name, 64);
	return s << b.flags << b.numChildren << b.parentIndex << b.bonePos;
}

Stream & nv::unreal::operator << (Stream & s, RawBoneInfluence & rbi)
{
	return s << rbi.weight << rbi.pointIndex << rbi.boneIndex;
}

Stream & nv::unreal::operator << (Stream & s, QuatAnimKey & qak)
{
	return s << qak.position << qak.orientation << qak.time;
}

Stream & nv::unreal::operator << (Stream & s, Vertex & v)
{
	uint16 pad; // Align struct elements to 4 bytes.
	return s << v.pointIndex << pad << v.u << v.v << v.matIndex << v.reserved << pad ;
}

Stream & nv::unreal::operator << (Stream & s, Point & p)
{
	return s << p.point;
}

Stream & nv::unreal::operator << (Stream & s, Triangle & t)
{
	s << t.wedgeIndex[0] << t.wedgeIndex[1] << t.wedgeIndex[2];	// @@ Can we handle arrays automatically?
	return s << t.matIndex << t.auxMatIndex << t.smoothingGroups;
}

