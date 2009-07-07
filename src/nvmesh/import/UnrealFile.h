// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_UNREALFILE_H
#define NV_MESH_IMPORT_UNREALFILE_H

// References:
//
// - Binary format specifications for skeletal and vertex animation source files
//   http://udn.epicgames.com/Two/BinaryFormatSpecifications.html
// 
// - C++ data structures
//   http://udn.epicgames.com/Two/rsrc/Two/BinaryFormatSpecifications/UnrealAnimDataStructs.h

#include <nvmath/Vector.h>
#include <nvmath/Quaternion.h>

namespace nv
{
	class Stream;
	

	namespace unreal
	{
		
		// A bone: an orientation, and a position, all relative to their parent.
		struct JointPos
		{
			Quaternion orientation;
			Vector3    position;     //  
			
			float      length;       //  For collision testing / debugging drawing.  (unused)
			float      xSize;
			float      ySize;
			float      zSize;
		};
		
		
		// Binary animation info format - used to organize raw animation keys into FAnimSeqs on rebuild
		// Similar to MotionChunkDigestInfo..
		struct AnimInfoBinary
		{
			char  name[64];     // Animation's name
			char  group[64];    // Animation's group name	
			
			int   totalBones;           // TotalBones * NumRawFrames is number of animation keys to digest.
			
			int   rootInclude;          // 0 none 1 included 	(unused)
			int   keyCompressionStyle;  // Reserved: variants in tradeoffs for compression.
			int   keyQuotum;            // Max key quotum for compression	
			float keyReduction;       // desired 
			float trackTime;          // explicit - can be overridden by the animation rate
			float animRate;           // frames per second.
			int   startBone;            // - Reserved: for partial animations (unused)
			int   firstRawFrame;        //
			int   numRawFrames;         // NumRawFrames and AnimRate dictate tracktime...
		};
		
		
		// File header structure. 
		struct ChunkHeader
		{
			char chunkID[20];  // String ID of up to 19 chars (usually zero-terminated)
			int  typeFlag;      // Flags/reserved
			int  dataSize;      // Size per struct following;
			int  dataCount;     // Number of structs/
		};
		
		// Raw data material.
		struct Material
		{
			char materialName[64];
			int  textureIndex;  // Texture index ('multiskin index')
			uint polyFlags;     // ALL poly's with THIS material will have this flag.
			int  auxMaterial;   // Reserved: index into another material, eg. detailtexture/shininess/whatever.
			uint auxFlags;      // Reserved: auxiliary flags 
			int  lodBias;       // Material-specific lod bias (unused)
			int  lodStyle;      // Material-specific lod style (unused)
		};
		
		
		// Raw data bone.
		struct Bone
		{
			char     name[64];     //
			uint     flags;        // Reserved.
			int      numChildren;  // Children  (not used.)
			int      parentIndex;  // 0/NULL if this is the root bone.  
			JointPos bonePos;      // Reference position.
		};
		
		
		// Raw data bone influence.
		struct RawBoneInfluence // Just weight, vertex, and Bone, sorted later.
		{
			float weight;
			int   pointIndex;
			int   boneIndex;
		};
		
		// An animation key.
		struct QuatAnimKey
		{
			Vector3    position;     // Relative to parent.
			Quaternion orientation;  // Relative to parent.
			float      time;         // The duration until the next key (end key wraps to first...)
		};
		
		// Vertex with texturing info, akin to Hoppe's 'Wedge' concept - import only.
		struct Vertex
		{
			uint16 pointIndex;	 // Index into the 3d point table.
			float  u, v;         // Texture U, V coordinates.
			uint8  matIndex;     // At runtime, this one will be implied by the face that's pointing to us.
			uint8  reserved;     // Top secret.
		};
		
		// Points: regular FVectors 
		struct Point
		{	
			Vector3  point; 
		};
		
		// Textured triangle.
		struct Triangle
		{
			uint16  wedgeIndex[3];	 // Point to three vertices in the vertex list.
			uint8   matIndex;	     // Materials can be anything.
			uint8   auxMatIndex;     // Second material (unused).
			uint32  smoothingGroups; // 32-bit flag for smoothing groups.
		};
		
		Stream & operator << (Stream &, JointPos & );
		Stream & operator << (Stream &, AnimInfoBinary & );
		Stream & operator << (Stream &, ChunkHeader & );
		Stream & operator << (Stream &, Material & );
		Stream & operator << (Stream &, Bone & );
		Stream & operator << (Stream &, RawBoneInfluence & );
		Stream & operator << (Stream &, QuatAnimKey & );
		Stream & operator << (Stream &, Vertex & );
		Stream & operator << (Stream &, Point & );
		Stream & operator << (Stream &, Triangle & );
		
		
	} // unreal namespace

} // nv namespace
	
#endif // NV_MESH_IMPORT_UNREALFILE_H
