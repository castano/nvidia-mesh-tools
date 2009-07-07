// This code is in the public domain -- castanyo@yahoo.es

#include <nvcore/StrLib.h>
#include <nvcore/Tokenizer.h>
#include <nvcore/StdStream.h>

#include <nvmath/Vector.h>

#include <nvmesh/animation/Bone.h>

#include "MeshImportPSK.h"
#include "UnrealFile.h"

using namespace nv;

namespace 
{

	struct TmpPoint
	{
		TmpPoint() : influences(4) {}
		Array<int> influences;
	};

}


/// Import a PSK mesh file.
bool MeshImportPSK::import(Stream * stream)
{
	ProtectedStream s(stream, false);
	
	Array<unreal::Point> pointArray;
	Array<unreal::Vertex> vertexArray;
	Array<unreal::Triangle> faceArray;
	Array<unreal::Material> materialArray;
	Array<unreal::Bone> boneArray;
	Array<unreal::RawBoneInfluence> influenceArray;
	
	try
	{
		// General Header    When the Typeflag value equals decimal 1999801 or lower it denotes this version 1.0 layout of the PSK file. 
		unreal::ChunkHeader header;
		s << header;
		
		if (header.typeFlag > 1999801)
		{
			// @@ Display error. Throw exception?
			return false;
		}
		
		// Points Header    Header specifying number of points, defined in Datacount. 
		unreal::ChunkHeader pointsHeader;
		s << pointsHeader;
		
		const uint pointCount = pointsHeader.dataCount;
		pointArray.resize(pointCount);
		
		// Points Data    Array with VPoint . 		
		for (uint p = 0; p < pointCount; p++)
		{
			s << pointArray[p];
		}
		
		// Wedges Header    Header specifying amount of wedges. 
		unreal::ChunkHeader wedgesHeader;
		s << wedgesHeader;
		
		const uint vertexCount = wedgesHeader.dataCount;
		vertexArray.resize(vertexCount);
		
		// Wedges Data    VVertex wedges array (a wedge consists of a UV pair with an index into the 3d points array.)		
		for (uint v = 0; v < vertexCount; v++)
		{
			s << vertexArray[v];
		}
		
		// Faces Header    Header specifying amount of faces. 
		unreal::ChunkHeader facesHeader;
		s << facesHeader;
		
		const uint faceCount = facesHeader.dataCount;
		faceArray.resize(faceCount);
		
		// Faces Data    VTriangle faces array. 
		for (uint f = 0; f < faceCount; f++)
		{
			s << faceArray[f];
		}
		
		// Materials Header    Header specifying amount of materials. 
		unreal::ChunkHeader materialsHeader;
		s << materialsHeader;
		
		const uint materialCount = materialsHeader.dataCount;
		materialArray.resize(materialCount);
		
		// Materials Data    VMaterial Materials array. 
		for (uint m = 0; m < materialCount; m++)
		{
			s << materialArray[m];
		}
		
		// Bones Header    Header specifying amount of bones. 
		unreal::ChunkHeader bonesHeader;
		s << bonesHeader;
		
		const uint boneCount = bonesHeader.dataCount;
		boneArray.resize(boneCount);
		
		// Bones Data    VBone bones array. 
		for (uint b = 0; b < boneCount; b++)
		{
			s << boneArray[b];
		}
		
		// Influences Header    Header specifying amount of bone influences.
		unreal::ChunkHeader influencesHeader;
		s << influencesHeader;
		
		const uint influenceCount = influencesHeader.dataCount;
		influenceArray.resize(influenceCount);
		
		// Influences Data    VRawBoneInfluence array of Influences.		
		for (uint i = 0; i < influenceCount; i++)
		{
			s << influenceArray[i];
		}
	}
	catch (...)
	{
		// @@ Display error.
		
		return false;
	}

	// Print info.
	nvDebug("point count = %d", pointArray.count());
	nvDebug("vertex count = %d", vertexArray.count());
	nvDebug("bone count = %d", boneArray.count());
	nvDebug("influence count = %d", influenceArray.count());
	nvDebug("face count = %d", faceArray.count());
	nvDebug("material count = %d", materialArray.count());


	// Build skeleton.
	const uint boneCount = boneArray.count();
	for (uint b = 0; b < boneCount; b++)
	{
		String name = boneArray[b].name;
		int parent = boneArray[b].parentIndex;
		
		Quaternion q = boneArray[b].bonePos.orientation;
		Vector3 v = boneArray[b].bonePos.position;

		// @@ This is weird and is not right. The legs of the characters are far to the front.

		// According to Ratcliff:
		v = Vector3(v.x(), -v.y(), v.z());

		if (b > 0)
		{
			// According to Ratcliff:
			q = Quaternion(-q.x(), q.y(), -q.z(), -q.w());
		}
		else
		{
			// According to Ratcliff:
			q = Quaternion(q.x(), -q.y(), q.z(), -q.w());

			q = q * axisAngle(Vector3(0, 0, 1), PI);
		}

		if (parent > (int)b)
		{
			nvDebugBreak();
			// This happens in:
			// - EgyptMaleA.psk
			if (parent == 2 && b == 1) parent = 0;
			// - EgyptFemaleB.psk
			if (parent == 22 && b == 21) parent = 2;
		}

		m_skeleton.addRelativeBone(name, parent, q, v);
	}


	// Build array of points, with list of influences.
	const uint pointCount = pointArray.count();
	Array<TmpPoint> tmpArray(pointCount);
	tmpArray.resize(pointCount);

	const uint influenceCount = influenceArray.count();
	for (uint i = 0; i < influenceCount; i++)
	{
		tmpArray[influenceArray[i].pointIndex].influences.append(i);
	}

	// Add influences in vertex order.
	m_skeleton.linkArray().reserve(influenceCount);
	m_skeleton.vertexArray().resize(vertexArray.count());
	
	const uint vertexCount = vertexArray.count();
	for (uint v = 0; v < vertexCount; v++)
	{
		const int pointIdx = vertexArray[v].pointIndex;

		const Vector3 point = pointArray[pointIdx].point;
		m_skeleton.vertexAt(v).pos = point;
		m_skeleton.vertexAt(v).link_first = m_skeleton.linkArray().count();

		const uint vertexInfluenceCount = tmpArray[pointIdx].influences.count();
		for (uint i = 0; i < vertexInfluenceCount; i++)
		{
			const int vertexInfluenceIdx = tmpArray[pointIdx].influences[i];

			Skeleton::Link link;
			link.bone = influenceArray[vertexInfluenceIdx].boneIndex;
			link.weight = influenceArray[vertexInfluenceIdx].weight;

			nvDebugCheck(link.bone < boneCount);

			// Transform point to bone space.
			const Matrix objectToBone = m_skeleton.boneAt(link.bone).inverse;
			link.offset = transformPoint(objectToBone, point);

			m_skeleton.linkArray().append(link);
		}

		m_skeleton.vertexAt(v).link_num = m_skeleton.linkArray().count() - m_skeleton.vertexAt(v).link_first;
	}

	// @@ Normalize weights? -> Add method to Skeleton class.

	tmpArray.clear();
	tmpArray.tighten();


	// Reset mesh.
	m_builder.reset();
	
	// Add positions.
	m_builder.hintPositionCount(pointCount);
	for (uint p = 0; p < pointCount; p++)
	{
		m_builder.addPosition(pointArray[p].point);
	}

	// Add texcoords.
	m_builder.hintTexCoordCount(vertexCount);
	for (uint v = 0; v < vertexCount; v++)
	{
		m_builder.addTexCoord(Vector2(vertexArray[v].u, vertexArray[v].v));
	}

	// Add faces.
	const uint faceCount = faceArray.count();
	m_builder.hintTriangleCount(faceCount);
	//m_builder.hintMaterialCount(materialArray.count());
	
	for (uint f = 0; f < faceCount; f++)
	{
		m_builder.beginPolygon();

		for (int i = 0; i < 3; i++)
		{
			m_builder.addVertex(vertexArray[faceArray[f].wedgeIndex[i]].pointIndex, NIL, faceArray[f].wedgeIndex[i]);
		}

		m_builder.endPolygon();
	}

	// @@ Add materials!
//	Array<unreal::Material> materialArray;
	
	m_builder.done();
	
	return true;
}



