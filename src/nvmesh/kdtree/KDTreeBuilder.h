// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#ifndef NV_MESH_KDTREEBUILDER_H
#define NV_MESH_KDTREEBUILDER_H

#include <nvcore/Containers.h>
#include <nvcore/Radix.h>

#include <nvmath/Triangle.h>

#include <nvmesh/kdtree/KDTree.h>

namespace nv
{
	class FaceLinkList;
	
	class KDTreeBuilder
	{
		NV_FORBID_COPY(KDTreeBuilder);
		NV_FORBID_HEAPALLOC();
	public:

		struct Options {
			Options() { setDefaults(); }
			NVMESH_API void setDefaults();
			
			uint leaf_face_count;
			uint max_level;
			float node_traversal_cost;
			float face_intersection_cost;
			float max_split_faces;
		};
		
		NVMESH_API KDTreeBuilder() {}
		
		// Create a kdtree.
		NVMESH_API KDTree * createTree(const Array<uint> & indices, const Array<Vector3> & vertices);

		// Options:
		Options & options() { return m_options; }

	private:
		
		struct Node;
		
		// Build methods.
		void computeBounds(const Array<Vector3> & vertices);
		void buildFaces(const Array<uint> & indices, const Array<Vector3> & vertices);

		void initStack();
		
		void buildEvents();
		
		void buildTree_r(uint level, uint idx, FaceLinkList  & head, const Box & bounds);
		bool selectSplit(uint level, const FaceLinkList & faceList, const Box & bounds, Node * best);

		bool selectSplitAccurate(uint level, const FaceLinkList & faceList, const Box & bounds, Node * node);
		bool selectSplitFast(uint level, const FaceLinkList & faceList, const Box & bounds, Node * node);

		uint addLeafFaces(const FaceLinkList & faceList);

		void splitFaces(const Node & node, const Box & leftBounds, const Box & rightBounds, FaceLinkList & head, FaceLinkList & leftLinkList, FaceLinkList & rightLinkList);
		
		void splitBounds(const Node & node, const Box & bounds, Box * box1, Box * box2);
		
		
		// Output methods.
		uint treeSize() const;	
		void writeNodes(KDTree * tree) const;
		void writeLeafFaces(KDTree * tree) const;
		void writeFaces(const Array<uint> & indices, const Array<Vector3> & vertices, KDTree * tree) const;

	private:

		Options m_options;

		Box m_bounds;
		
		// By now, only need the bounds during construction.
		struct Face {
			Triangle winding;
			Box bounds;
		};

		Array<Face> m_face_array;

		// @@ TODO: Make node more compact to save mem during construction.
		struct Node {
			int axis;		// 0 = x, 1 = y, 2 = z, -1 = leaf.
			float split;	// split distance.
			uint son[2];
		};
		
		Array<Node> m_node_array;
		Array<uint> m_leaf_face_array;

		// Statistics:
		struct Stats {
			Stats() : leaf_count(0), leaf_face_max(0), empty_leaf_count(0) { }
			uint leaf_count;
			uint leaf_face_max;
			uint empty_leaf_count;
		} stats;

		// Temporary members:
		Array<float> m_plane_array;
		RadixSort m_radix;

		uint64 time_buildTree;
		uint64 time_selectSplit;
		uint64 time_selectSplitAccurate;
		uint64 time_initPlanes;
		uint64 time_sortPlanes;
		uint64 time_sweepAndPrune;
		uint64 time_splitFaces;
		uint64 time_splitBounds;
	};

} // nv namespace

#endif // NV_MESH_KDTREEBUILDER_H
