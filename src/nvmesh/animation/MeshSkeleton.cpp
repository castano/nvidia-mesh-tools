// This code is in the public domain -- castanyo@yahoo.es

#include "MeshSkeleton.h"
#include "Bone.h"

using namespace nv;



uint Skeleton::addAbsoluteBone(String name, int parent, Quaternion::Arg q, Vector3::Arg v)
{
	const int boneIdx = bone_array.count();

	// First bone must be the root.
	if (boneIdx == 0)
	{
		nvCheck(parent == 0 || parent == -1);
		parent = -1;
	}

	nvCheck(parent < boneIdx);

	Bone bone;
	bone.name = name;
	bone.parent = parent;
	bone.rotation = q;
	bone.offset = v;

	bone.inverse = isometryInverse(boneMatrix(bone.rotation, bone.offset));
	//bone.inverse = boneMatrix(conjugate(bone.rotation), -bone.offset);

	bone_array.append(bone);

	return boneIdx;
}

uint Skeleton::addRelativeBone(String name, int parent, Quaternion::Arg q, Vector3::Arg v)
{
	const int boneIdx = bone_array.count();

	// First bone must be the root.
	if (boneIdx == 0)
	{
		nvCheck(parent == 0 || parent == -1);
		parent = -1;
	}

	nvCheck(parent < boneIdx);

	Bone bone;
	bone.name = name;
	bone.parent = parent;

	// @@ Transform q,v to parent space.
	if (parent == -1)
	{
		bone.rotation = normalize(q);
		bone.offset = v;
	}
	else
	{
		const Bone & parentBone = bone_array[parent];
		bone.rotation = normalize(parentBone.rotation * q);
		bone.offset = parentBone.offset + transform(parentBone.rotation, v);
	}

	bone.inverse = isometryInverse(boneMatrix(bone.rotation, bone.offset));
	//bone.inverse = boneMatrix(conjugate(bone.rotation), -bone.offset);

	bone_array.append(bone);

	return boneIdx;
}


/// Reset the skeleton.
void Skeleton::reset()
{
	bone_array.clear();
	link_array.clear();
	vertex_array.clear();
}


/// Set default pose.
void Skeleton::setDefaultPose()
{
	Array<Matrix> matrix_array;
	buildDefaultMatrixArray(matrix_array);
	
	setPose(matrix_array);
}


/// Set the given pose.
void Skeleton::setPose(const Array<Matrix> & matrix_array)
{
	nvCheck(matrix_array.count() == bone_array.count());
	
	const uint vertex_num = vertex_array.count();

	for(uint v = 0; v < vertex_num; v++)
	{
		const Vertex & vertex = vertex_array[v];	
		
		Vector3 pos(zero);
		
		// @@ It's probably faster to blend the matrices and transform point only once.
		for(uint l = 0; l < vertex.link_num; l++)
		{
			const Link & link = link_array[vertex.link_first + l];
			
			pos += transformPoint(matrix_array[link.bone], link.offset) * link.weight;
		}
		
		// Set position.
		vertex_array[v].pos = pos;
	}
}


void Skeleton::buildDefaultMatrixArray(Array<Matrix> & matrix_array) const
{
	const int bone_num = bone_array.count();
	matrix_array.resize(bone_num);

	nvCheck( bone_array[0].parent == -1 );
	matrix_array[0] = boneMatrix(bone_array[0].rotation, bone_array[0].offset);

	// Convert bones to matrices.
	for(int b = 1; b < bone_num; b++)
	{
		const Bone & bone = bone_array[b];
		nvDebugCheck( bone.parent >= 0 );
		nvDebugCheck( bone.parent < b );
		
		matrix_array[b] = boneMatrix(bone.rotation, bone.offset);
	}
}

/// 
void Skeleton::buildSkinnedVertices(Array<SkinVertex> & skinVertexArray) const
{
	Array<Matrix> matrix_array;
	buildDefaultMatrixArray(matrix_array);
	
	const uint vertex_num = vertex_array.count();
	skinVertexArray.resize(vertex_num);
	
	for(uint v = 0; v < vertex_num; v++)
	{
		const Vertex & vertex = vertex_array[v];
		
		SkinVertex & skinVertex = skinVertexArray[v];
		
		// Reset vertex
		for (int i = 0; i < 4; i++)
		{
			skinVertex.bone[i] = -1;
			skinVertex.weight[i] = 0.0f;
		}
		skinVertex.pos = Vector3(zero);
		
		for(uint l = 0; l < vertex.link_num; l++)
		{
			const Link & link = link_array[vertex.link_first + l];
			
			if (l < 4)
			{
				skinVertex.weight[l] = link.weight;
				skinVertex.bone[l] = link.bone;
			}
			else
			{
				// Find link to replace.
				for (int i = 0; i < 4; i++)
				{
					if (link.weight > skinVertex.weight[i])
					{
						skinVertex.weight[i] = link.weight;
						skinVertex.bone[l] = link.bone;
					}
				}
			}
			
			skinVertex.pos += transformPoint(matrix_array[link.bone], link.offset) * link.weight;
			
			// Bubble sort weights.
			if (skinVertex.weight[1] > skinVertex.weight[0])
			{
				swap(skinVertex.weight[1], skinVertex.weight[0]);
				swap(skinVertex.bone[1], skinVertex.bone[0]);
			}
			if (skinVertex.weight[2] > skinVertex.weight[1])
			{
				swap(skinVertex.weight[2], skinVertex.weight[1]);
				swap(skinVertex.bone[2], skinVertex.bone[1]);
			}
			if (skinVertex.weight[3] > skinVertex.weight[2])
			{
				swap(skinVertex.weight[3], skinVertex.weight[2]);
				swap(skinVertex.bone[3], skinVertex.bone[2]);
			}
			if (skinVertex.weight[1] > skinVertex.weight[0])
			{
				swap(skinVertex.weight[1], skinVertex.weight[0]);
				swap(skinVertex.bone[1], skinVertex.bone[0]);
			}
			if (skinVertex.weight[2] > skinVertex.weight[1])
			{
				swap(skinVertex.weight[2], skinVertex.weight[1]);
				swap(skinVertex.bone[2], skinVertex.bone[1]);
			}
			if (skinVertex.weight[1] > skinVertex.weight[0])
			{
				swap(skinVertex.weight[1], skinVertex.weight[0]);
				swap(skinVertex.bone[1], skinVertex.bone[0]);
			}
			
			nvDebugCheck(skinVertex.weight[0] >= skinVertex.weight[1]);
			nvDebugCheck(skinVertex.weight[1] >= skinVertex.weight[2]);
			nvDebugCheck(skinVertex.weight[2] >= skinVertex.weight[3]);
			
			nvDebugCheck(skinVertex.bone[0] != -1);
			// Optimize constant access?
			//if (skinVertex.bone[1] == -1) skinVertex.bone[1] = skinVertex.bone[0];
			//if (skinVertex.bone[2] == -1) skinVertex.bone[2] = skinVertex.bone[1];
			//if (skinVertex.bone[3] == -1) skinVertex.bone[3] = skinVertex.bone[2];
			if (skinVertex.bone[1] == -1) skinVertex.bone[1] = 0;
			if (skinVertex.bone[2] == -1) skinVertex.bone[2] = 0;
			if (skinVertex.bone[3] == -1) skinVertex.bone[3] = 0;

			// Normalize weights.
			const float sum = skinVertex.weight[0] + skinVertex.weight[1] + skinVertex.weight[2] + skinVertex.weight[3];
			skinVertex.weight[0] /= sum;
			skinVertex.weight[1] /= sum;
			skinVertex.weight[2] /= sum;
			skinVertex.weight[3] /= sum;
			
			nvDebugCheck(equal(skinVertex.weight[0] + skinVertex.weight[1] + skinVertex.weight[2] + skinVertex.weight[3], 1.0f));
		}
	}
	
	// @@ Count the average of the max influences per warp.
	
}


