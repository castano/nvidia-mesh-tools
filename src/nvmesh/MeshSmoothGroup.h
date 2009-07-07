// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_SMOOTHGROUP_H
#define NV_MESH_SMOOTHGROUP_H

namespace nv
{

/// Smooth group info can be represented in different ways.
/// - MAX style smooth groups.
/// - OBJ style smooth groups.
/// - Creases.
/// - Tagged meshes.
class MeshSmoothGroup
{
public:
	MeshSmoothGroup(uint faceNum, uint edgeNum);
	~MeshSmoothGroup();

	uint group(uint face);
	bool crease(uint edge);

public:
	bool useMaxStyleGroups;
};

} // nv namespace

#endif // NV_MESH_SMOOTHGROUP_H
