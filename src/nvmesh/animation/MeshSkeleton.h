// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_SKELETON_H
#define NV_MESH_SKELETON_H

#include <nvcore/StrLib.h>

#include <nvmath/Matrix.h>
#include <nvmath/Quaternion.h>

#include <nvmesh/nvmesh.h>

namespace nv
{

	class BaseMesh;

	class Skeleton
	{
	public:

		struct Bone
		{
			String name;
			int parent;
			Quaternion rotation;
			Vector3 offset;
			Matrix inverse;
		};

		struct Link
		{
			uint bone;
			float weight;
			Vector3 offset;
		};

		struct Vertex
		{
			uint link_first;
			uint link_num;
			Vector3 pos;
		};

		struct SkinVertex
		{
			int bone[4];
			float weight[4];
			Vector3 pos;
		};
		

	public:

		Skeleton() {}

		uint boneCount() const;
		const Bone & boneAt(uint id) const;
		Bone & boneAt(uint id);

		const Array<Bone> & boneArray() const;
		Array<Bone> & boneArray();
		
		const Link & linkAt(uint id) const;
		Link & linkAt(uint id);

		const Array<Link> & linkArray() const;
		Array<Link> & linkArray();
		
		const Vertex & vertexAt(uint id) const;
		Vertex & vertexAt(uint id);

		const Array<Vertex> & vertexArray() const;
		Array<Vertex> & vertexArray();


		uint addAbsoluteBone(String name, int parent, Quaternion::Arg q, Vector3::Arg v);
		uint addRelativeBone(String name, int parent, Quaternion::Arg q, Vector3::Arg v);

		void reset();
		void setDefaultPose();

		void setPose(const Array<Matrix> & matrix_array);

		void buildDefaultMatrixArray(Array<Matrix> & matrix_array) const;
		void buildSkinnedVertices(Array<SkinVertex> & array) const;

		// @@ Add method to normalize weights.


	private:

		/// Array of bones.
		Array<Bone> bone_array;

		/// Array of links.
		Array<Link> link_array;

		/// Array of vertices.
		Array<Vertex> vertex_array;
		
	};

	
	/// Get bone count.
	inline uint Skeleton::boneCount() const {
		return bone_array.size();
	}

	/// Get a bone given its index.
	inline const Skeleton::Bone & Skeleton::boneAt(uint id) const {
		nvDebugCheck(id < bone_array.size());
		return bone_array[id];
	}
	inline Skeleton::Bone & Skeleton::boneAt(uint id) {
		nvDebugCheck(id < bone_array.size());
		return bone_array[id];
	}

	/// Get a link given its index.
	inline const Skeleton::Link & Skeleton::linkAt(uint id) const {
		nvDebugCheck(id < link_array.size());
		return link_array[id];
	}
	inline Skeleton::Link & Skeleton::linkAt(uint id) {
		nvDebugCheck(id < link_array.size());
		return link_array[id];
	}

	/// Get a vertex given its index.
	inline const Skeleton::Vertex & Skeleton::vertexAt(uint id) const {
		nvDebugCheck(id < vertex_array.size());
		return vertex_array[id];
	}
	inline Skeleton::Vertex & Skeleton::vertexAt(uint id) {
		nvDebugCheck(id < vertex_array.size());
		return vertex_array[id];
	}


	inline const Array<Skeleton::Bone> & Skeleton::boneArray() const {
		return bone_array;
	}

	inline const Array<Skeleton::Link> & Skeleton::linkArray() const {
		return link_array;
	}

	inline const Array<Skeleton::Vertex> & Skeleton::vertexArray() const {
		return vertex_array;
	}


	inline Array<Skeleton::Bone> & Skeleton::boneArray() {
		return bone_array;
	}

	inline Array<Skeleton::Link> & Skeleton::linkArray() {
		return link_array;
	}

	inline Array<Skeleton::Vertex> & Skeleton::vertexArray() {
		return vertex_array;
	}


} // nv namespace


#endif // NV_MESH_SKELETON_H
