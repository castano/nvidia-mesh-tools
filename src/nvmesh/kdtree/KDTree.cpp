// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

/* This is partially based on Wald thesis:
 * http://www.mpi-sb.mpg.de/~wald/PhD/
 *
 * And Jacco Bikker raytracer:
 * http://www.flipcode.com/articles/article_raytrace07.shtml
 */

#include <nvcore/Prefetch.h>

#include <nvmesh/kdtree/KDTree.h>

#if NV_CC_GNUC
#include <stdint.h> // uintptr_t
#endif


#define COMPUTE_STATS	0
#define CLIP_TO_BOX		0

using namespace nv;


/// Cache a branch of the tree.
struct KDTree::Cache {
	uint node_count;
	Node node_array[48];	// max tree depth
};

/// Local information for each ray test
struct KDTree::Query {
	uint last;
	Item stack[32];
};

// @@ Make sure allocations are aligned.
void KDTree::allocateNodes(uint nodeCount)
{
	m_nodeCount = nodeCount;
	m_nodeBuffer = (uint8 *)new Node[nodeCount];
}

void KDTree::allocateLeafFaces(uint leafFaceCount)
{
	m_leafFaceCount = leafFaceCount;
	m_leafFaceBuffer = (uint8 *)new uint[leafFaceCount];
}


/// Free the tree memory.
void KDTree::free()
{
	m_nodeCount = 0;
	delete [] (Node *) m_nodeBuffer;
	m_nodeBuffer = NULL;

	m_leafFaceCount = 0;
	delete [] (uint *) m_leafFaceBuffer;
	m_leafFaceBuffer = NULL;

	delete m_faceBuffer;
	m_faceBuffer = NULL;
}


/// Test ray, boolean query only.
bool KDTree::testRay(const Ray & ray, Query * query) const
{
	nvDebugCheck(isReady());

	// Allocate a new query if required
	const bool isQueryInternal = query == NULL;
	if (isQueryInternal) {
		query = createQuery();
	}
	nvDebugCheck(query != NULL);

#if COMPUTE_STATS
	// Reset stats.
	m_face_test_count = 0;
	m_node_test_count = 0;
	m_leaf_test_count = 0;
#endif

	// Recurse down the tree.
	float t_near, t_far;
	t_near = 0; // NV_EPSILON;
	t_far = ray.maxt;

#if CLIP_TO_BOX
	// Clip segment against box.
	if( !m_bounds.clipSegment(ray.dir, ray.origin, &t_near, &t_far) ) {
		// Clipped away.
		return false;
	}
#endif

	// Init the stack.
	query->last = 0;
	query->stack[query->last].node   = nodeAt(0);
	query->stack[query->last].t_near = t_near;
	query->stack[query->last].t_far  = t_far;
	++(query->last);

	const bool result = testRayIterative(ray, query);

	if (isQueryInternal) {
		deleteQuery(query);
	}

	return result;
}


/*/// Test ray, get distance to first hit.
bool KDTree::testRay(const Ray & ray, float * t) const
{
	nvDebugCheck(t != NULL);
	nvDebugCheck(isReady());
	
	Hit hit;
	bool result = testRay(ray, &hit);
	*t = hit.t;
	
	return result;
}
*/

/// Test ray, get first hit.
void KDTree::testRay(const Ray & ray, Hit * hit, Query * query) const
{
	nvDebugCheck(hit != NULL);
	nvDebugCheck(isReady());

	// Allocate a new query if required
	const bool isQueryInternal = query == NULL;
	if (isQueryInternal) {
		query = createQuery();
	}
	nvDebugCheck(query != NULL);

	// Init result.
	hit->face = ~0;
	hit->t = ray.maxt;

#if COMPUTE_STATS
	// Reset stats.
	m_face_test_count = 0;
	m_node_test_count = 0;
	m_leaf_test_count = 0;
#endif

	// Recurse down the tree.
	float t_near, t_far;
	t_near = 0; // NV_EPSILON;
	t_far = ray.maxt;

#if CLIP_TO_BOX
	// Clip segment against box.
	if( !m_bounds.ClipSegment(ray.dir, ray.origin, &t_near, &t_far) ) {
		// Clipped away.
		return false;
	}
#endif

	// Init the stack.
	query->last = 0;
	query->stack[query->last].node   = nodeAt(0);
	query->stack[query->last].t_near = t_near;
	query->stack[query->last].t_far  = t_far;
	++(query->last);

	testRayIterative(ray, hit, query);

	if (isQueryInternal) {
		deleteQuery(query);
	}
}


/// Test array agains the tree using the ray cache.
void KDTree::testRay(const Cache & cache, const Ray & ray, Hit * hit, Query * query) const
{
	nvDebugCheck(hit != NULL);
	nvDebugCheck(isReady());

	// Allocate a new query if required
	const bool isQueryInternal = query == NULL;
	if (isQueryInternal) {
		query = createQuery();
	}
	nvDebugCheck(query != NULL);

	// Init result.
	hit->face = ~0;
	hit->t = ray.maxt;

#if COMPUTE_STATS
	// Reset stats.
	m_face_test_count = 0;
	m_node_test_count = 0;
	m_leaf_test_count = 0;
#endif

	// Recurse down the tree.
	float t_near, t_far;
	t_near = 0; // NV_EPSILON;
	t_far = ray.maxt;

#if CLIP_TO_BOX
	// Clip segment against box.
	if( !m_bounds.ClipSegment(ray.dir, ray.origin, &t_near, &t_far) ) {
		// Clipped away.
		return false;
	}
#endif

	nvCheck(cache.node_count < 48);

	// Init the stack.
	query->last = 0;

	if (cache.node_count == 0) {
		// Cache is empty.
		return;
	}

	// Clip the ray by the cached branch.
	for(uint n = 0; n < cache.node_count - 1; n++)
	{
		// @todo Unroll and prefetch this loop.
		//nvPrefetch(&cache.node_array[n+2]);

		const Node & node = cache.node_array[n];
		const uint axis = node.dimension();

		const float d = node.distance * ray.idir.component(axis);

		if( d > t_near && d < t_far ) {
			const uint offset = node.offset();			

			// We have to visit the opposite node later, push it on the stack.
			const Node * far_node = nodeAt(offset);
			if( ray.dir.component(axis) >= 0 ) {
				far_node++;
			}

			query->stack[query->last].node   = far_node;
			query->stack[query->last].t_near = t_near;
			query->stack[query->last].t_far  = t_far;
			++(query->last);

			t_far = d;
		}
	}

	// Push leaf on the stack.
	const Node & node = cache.node_array[cache.node_count - 1];
	nvDebugCheck(node.isLeaf());

	query->stack[query->last].node   = &node;
	query->stack[query->last].t_near = t_near;
	query->stack[query->last].t_far  = t_far;
	++(query->last);

	// Start normal traversal.
	testRayIterative(ray, hit, query);

	if (isQueryInternal) {
		deleteQuery(query);
	}
}


/// Iterative traversal.
void KDTree::testRayIterative(const Ray & ray, Hit * hit, Query * query) const
{
	nvDebugCheck(query != NULL);

	do
	{
		// Pop an element from the top of the stack.
		--(query->last);
		const Node * node = query->stack[query->last].node;
		float t_near = query->stack[query->last].t_near;
		float t_far  = query->stack[query->last].t_far;
		
		while( !node->isLeaf() )
		{
#if COMPUTE_STATS
			m_node_test_count++;
#endif
			const uint axis = node->dimension();
			const uint offset = node->offset();

			// Compute ray intersection.
			const float d = (node->distance - ray.origin.component(axis)) * ray.idir.component(axis);

			const Node * near_node = nodeAt(offset);
			const Node * far_node = near_node + 1;
			
			if( ray.dir.component(axis) < 0 ) {
				near_node++;
				far_node--;
			}

			// Traverse node.
			if( d <= t_near ) {
				node = far_node;
			}
			else if( d >= t_far ) {
				node = near_node;
			}
			else {
				query->stack[query->last].node = far_node;
				query->stack[query->last].t_near = d;
				query->stack[query->last].t_far = t_far;
				++(query->last);
				nvDebugCheck(query->last < 32);

				t_far = d;
				node = near_node;
			}

			nvPrefetch(m_nodeBuffer + node->offset());
		};

		// Test leaf.
		const Leaf * leaf = reinterpret_cast<const Leaf *>( node );

#if COMPUTE_STATS
		m_leaf_test_count++;
		m_face_test_count += leaf->face_count;
#endif

		if (leaf->face_count != 0)
		{
			m_faceBuffer->testRay(leafFaces(node->offset()), leaf->face_count, ray, hit);
			
			// If the hit is inside the current leaf we are done.
			if (hit->t < t_far)
			{
				return;
			}
		}
	} while( query->last != 0 );

	// No intersection.
	return;
}


/** Iterative traversal. Does not compute the hit point, but returns as soon as an 
 * intersection is found. This is a little bit faster than the other method and is
 * useful for shadow rays.
 */
bool KDTree::testRayIterative(const Ray & ray, Query * query) const
{
	do {
		// Pop an element from the top of the stack.
		--(query->last);
		const Node * node = query->stack[query->last].node;
		float t_near = query->stack[query->last].t_near;
		float t_far  = query->stack[query->last].t_far;

		while( !node->isLeaf() ) {
#if COMPUTE_STATS
			m_node_test_count++;
#endif
			const uint axis = node->dimension();
			const uint offset = node->offset();	

			// Compute ray intersection.
			const float d = (node->distance - ray.origin.component(axis)) * ray.idir.component(axis);

			const Node * near_node = nodeAt(offset);
			const Node * far_node = nodeAt(offset) + 1;
			if( ray.dir.component(axis) < 0 ) {
				//nvSwap(far_node, near_node);
				near_node++;
				far_node--;
			}

			// Traverse node.
			if( d <= t_near ) {
				node = far_node;
			}
			else if( d >= t_far ) {
				node = near_node;
			}
			else {
				query->stack[query->last].node = far_node;
				query->stack[query->last].t_near = d;
				query->stack[query->last].t_far = t_far;
				++(query->last);

				nvDebugCheck(query->last < 32);

				t_far = d;
				node = near_node;
			}
		};


		// Test leaf, return as soon as we find an intersection.
		const Leaf * leaf = reinterpret_cast<const Leaf *>( node );

#if COMPUTE_STATS
		m_leaf_test_count++;
		m_face_test_count += leaf->face_count;
#endif

		if( leaf->face_count != 0 )
		{
			if (m_faceBuffer->testRay(leafFaces(node->offset()), leaf->face_count, ray))
			{
				return true;
			}
		}
	} while( query->last != 0 );

	// No intersection.
	return false;
}


/// Test ray against all the faces, get first hit.
void KDTree::testAllFaces(const Ray & ray, Hit * hit) const
{
	// Init result.
	hit->face = ~0;
	hit->t = ray.maxt;
	
	const uint faceCount = m_faceBuffer->faceCount();
	
	m_faceBuffer->testRay(0U, faceCount, ray, hit);
}


KDTree::Cache * KDTree::createCache() const
{
	// @@ Allocate nodes according to max depth of the tree.
	return new Cache();
}

void KDTree::deleteCache(Cache * cache) const
{
	delete cache;
}

KDTree::Query * KDTree::createQuery() const
{
	// @@ Allocates a stack to be used during transversal
	Query * query = new Query();
	query->last = 0;
	query->stack[0].node   = NULL;
	query->stack[0].t_near = 0.0f;
	query->stack[0].t_far  = 0.0f;
	return query;
}

void KDTree::deleteQuery(Query * query) const
{
	delete query;
}


/** Create a spherical ray cache. 
 * The ray origin is the center of the sphere, and maxt is its radius.
 */
void KDTree::initSphericalCache(const Ray & ray, Cache * cache) const
{
	nvDebugCheck(isReady());
	nvDebugCheck(cache != NULL);

	cache->node_count = 0;

	const Node * node = nodeAt(0);

	while( !node->isLeaf() )
	{
		const uint axis = node->dimension();

		// Compute distance to ray origin.
		const float d = (node->distance - ray.origin.component(axis));

		// Determine if the cache volume intersects with the node.
		if( fabs(d) < ray.maxt ) {
			// Push node.
			cache->node_array[cache->node_count].flags = node->flags;
			cache->node_array[cache->node_count].distance = d; // store relative distance.
			cache->node_count++;
		}

		// Choose next node.
		if( d > 0 ) {
			node = leftSideSon(node);
		}
		else {
			node = rightSideSon(node);
		}
	}

	// Push leaf.
	cache->node_array[cache->node_count] = *node;
	cache->node_count++;
	nvDebugCheck(cache->node_count < 48); // @@ Use max tree depth.
}


/** Create a hemispherical ray cache. 
 * The ray origin is the center of the hemisphere, maxt is its radius and dir
 * is the axis of the hemisphere.
 */
void KDTree::initHemiSphericalCache(const Ray & ray, Cache * cache) const
{
	nvDebugCheck(isReady());
	nvDebugCheck(cache != NULL);
	
	cache->node_count = 0;

	const Node * node = nodeAt(0);

	while( !node->isLeaf() ) {

		const float distance = node->distance;
		const uint axis = node->dimension();

		// Compute distance to ray origin.
		const float d = (distance - ray.origin.component(axis));

		// Determine if the cache volume intersects with the node.
		if( fabs(d) < ray.maxt ) {

			// @@ TODO 

			// if( d * (dir dot plane) < 0 ) -> if( d * dir[axis] < 0 )
			// or if( sqrt(1 - SQ(ray.dir.z)) < ray.maxt )

			// Push node.
			cache->node_array[cache->node_count].flags = node->flags;
			cache->node_array[cache->node_count].distance = d;
			cache->node_count++;
			nvDebugCheck(cache->node_count < 32);
		}

		// Choose next node.
		if( d < 0 ) {
			node = leftSideSon(node);
		}
		else {
			node = rightSideSon(node);
		}
	}

	// Push leaf.
	cache->node_array[cache->node_count] = *node;
	cache->node_count++;
	nvDebugCheck(cache->node_count < 48); // @@ Use max tree depth.
}


namespace nv
{
	Stream & operator<< (Stream & s, KDTree::Node & node)
	{
		return s << node.flags << node.distance;
	}

	Stream & operator<< (Stream & s, KDTree & tree)
	{
		if (s.isLoading())
		{
			tree.free();
		}

		s << tree.m_bounds;

		// Serialize nodes.
		if (s.isLoading())
		{
			uint nodeCount;
			s << nodeCount;
			tree.allocateNodes(nodeCount);
		}
		else
		{
			s << tree.m_nodeCount;
		}

		for (uint i = 0; i < tree.m_nodeCount; i++)
		{
			s << tree.nodeBuffer()[i];
		}

		// Serialize leaf faces.
		if (s.isLoading())
		{
			uint leafFaceCount;
			s << leafFaceCount;
			tree.allocateLeafFaces(leafFaceCount);
		}
		else
		{
			s << tree.m_leafFaceCount;
		}

		for (uint i = 0; i < tree.m_leafFaceCount; i++)
		{
			s << tree.leafFaceBuffer()[i];
		}

		// Serialize faces.
		if (s.isLoading())
		{
			tree.m_faceBuffer = new WaldFaceBuffer();
		}

		s << *(WaldFaceBuffer *)tree.m_faceBuffer;

		return s;
	}

}

