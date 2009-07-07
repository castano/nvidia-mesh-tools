// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#ifndef NV_MESH_KDTREE_H
#define NV_MESH_KDTREE_H

#include <nvcore/Ptr.h>

#include <nvmath/Box.h>

#include <nvmesh/nvmesh.h>
#include <nvmesh/raytracing/FaceBuffer.h>


namespace nv
{

/// KDTree acceleration structure.
class KDTree// : public RefCounted
{
	friend class KDTreeBuilder;
public:

	struct Cache;
	struct Query;

	/// Ctor.
	KDTree() :
		m_nodeCount(0),
		m_nodeBuffer(NULL),
		m_leafFaceCount(0),
		m_leafFaceBuffer(NULL),
		m_faceBuffer(NULL)
	{
		// Enabled or not, reset the stats.
		m_node_test_count = 0;
		m_leaf_test_count = 0;
		m_face_test_count = 0;
	}

	/** @name Tree queries. */
	//@{
		NVMESH_API bool testRay(const Ray & ray, Query * query = NULL) const;
		NVMESH_API void testRay(const Ray & ray, Hit * hit, Query * query = NULL) const;

		//NVMESH_API bool testRay(const Cache & cache, const Ray & ray) const;
		NVMESH_API void testRay(const Cache & cache, const Ray & ray, Hit * hit, Query * query = NULL) const;

		// Debug method.
		void testAllFaces(const Ray & ray, Hit * hit) const;

		NVMESH_API Cache * createCache() const;
		NVMESH_API void deleteCache(Cache * ) const;

		NVMESH_API Query * createQuery() const;
		NVMESH_API void deleteQuery(Query * ) const;

		NVMESH_API void initSphericalCache(const Ray & ray, Cache * cache) const;
		NVMESH_API void initHemiSphericalCache(const Ray & ray, Cache * cache) const;

		Box bounds() const { return m_bounds; }
	//@}

	/** @name Tree allocation. */
    //@{
		NVMESH_API void allocateNodes(uint nodeCount);
		NVMESH_API void allocateLeafFaces(uint leafFaceCount);
		NVMESH_API void allocateFaces(uint faceCount);
		NVMESH_API void free();
	//@}


private:

	struct Node;
	struct Leaf;
	struct Face;

	// Node tests.
	void testRayIterative(const Ray & ray, Hit * hit, Query * query) const;
	bool testRayIterative(const Ray & ray, Query * query) const;

	// Face tests.
	bool testFace(const uint f, const Ray & ray, Hit * hit) const;
	bool testFaceDoubleSided(const uint f, const Ray & ray, Hit * hit) const;

	/// Returns true if the tree is built, returns false otherwise.
	bool isReady() const {
		return m_nodeBuffer != NULL;
	}


	/// Get node pointer for the given offset.
	const Node * nodeAt(const uint offset) const {
		nvDebugCheck(offset / sizeof(Node) < m_nodeCount);
		return reinterpret_cast<Node *>(m_nodeBuffer + offset);
	}

	/// Get left side child node.
	const Node * leftSideSon(const Node * node) const {
		return nodeAt(node->offset());
	}

	/// Get right side child node.
	const Node * rightSideSon(const Node * node) const {
		return nodeAt(node->offset()) + 1;
	}

	/// Get leaf face pointer for the given offset.
	const uint * leafFaces(const uint offset) const {
		nvDebugCheck(offset / sizeof(uint) < m_leafFaceCount);
		return reinterpret_cast<uint *>(m_leafFaceBuffer + offset);
	}

	friend Stream & operator<< (Stream & s, KDTree & tree);


private:

	friend Stream & operator<< (Stream & s, Node & node);

	/// Get node buffer.
	Node * nodeBuffer() {
		return reinterpret_cast<Node *>(m_nodeBuffer);
	}

	/// Get leaf face buffer.
	uint * leafFaceBuffer() {
		return reinterpret_cast<uint *>(m_leafFaceBuffer);
	}

	/// Stack item.
	struct Item {
		Item() {}
		Item(const Item & i) : node(i.node), t_near(i.t_near), t_far(i.t_far) { }
		Item(const Node * b, float n, float f) : node(b), t_near(n), t_far(f) { }

		const Node * node;
		float t_near;
		float t_far;
	};

	/// Inner node of the tree. This is public, because is used in the cache.
	struct Node
	{
		bool isLeaf() const {
			return (flags & uint(1<<31)) != 0;	// flags < 0
		}
		
		uint dimension() const {
			return flags & 0x3;
		}

		uint offset() const {
			return flags & 0x7FFFFFFC;
		}

		uint32 flags;
		// bits 0..1 : splitting dimension
		// bits 2..30 : offset bits
		// bit 31 (sign) : flag whether node is a leaf.
		
		float distance;
	};

	/// Leaf of the tree.
	struct Leaf {
		
		uint32 flags;
		// bits 0..1 : face types. (not used)
		// bits 2..30 : offset to first son.
		// bit 31 (sign) : flag whether node is a leaf.
		
		uint32 face_count;
	};


	Box m_bounds;

	uint m_nodeCount;
	uint8 * m_nodeBuffer;

	uint m_leafFaceCount;
	uint8 * m_leafFaceBuffer;

	FaceBuffer * m_faceBuffer;

	// Traversal state. Changes during the queries. @@ Move this to a query object. KD-Tree should be reentrant.
	//mutable uint m_last;
	//mutable Item m_stack[32];

public:
	// Statistics.
	mutable uint m_node_test_count;
	mutable uint m_leaf_test_count;
	mutable uint m_face_test_count;
};

} // nv namespace

#endif // NV_MESH_KDTREE_H
