// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#include <algorithm>

#include <nvcore/Memory.h>
#include <nvcore/CpuInfo.h>

#include <nvmath/Polygon.h>

#include "KDTreeBuilder.h"

using namespace nv;

namespace {

	// Align offset to cache lines (32 or 64 bytes).
	static uint AlignOffset(const uint offset) {
		return (offset + 31) & ~31;
	}
	
} // namespace

namespace nv
{
	
	class FaceLink : public Link
	{
	public:
		FaceLink(uint f) : face(f) {}
		
		uint face;
	};
	
	class FaceLinkList
	{
	public:
		FaceLinkList() : m_head(NULL), m_count(0) {}
	
		void append(uint f)
		{
			attach(new FaceLink(f));
		}
		
		void attach(FaceLink * link)
		{
			nvDebugCheck(link != NULL);
			if (m_head == NULL) 
			{
				m_head = link;
			}
			else
			{
				m_head->addBefore(link);
			}
			m_count++;
			s_totalCount++;
		}
		
		FaceLink * detach()
		{
			FaceLink * link = m_head;
			m_head = (FaceLink *)m_head->next;
			link->remove();
			m_count--;
			s_totalCount--;
			
			if (m_head == link)
			{
				nvDebugCheck(m_count == 0);
				m_head = NULL;
			}

			return link;
		}

		void advance(FaceLink * & link) const
		{
			link = (link->next == m_head) ? NULL : (FaceLink *)link->next;
		}
		
		void deleteAll()
		{
			while (!isEmpty())
			{
				FaceLink * link = detach();
				delete link;
			}
		}
		
		FaceLink * head() const { return m_head; }
		uint count() const { return m_count; }
		
		bool isEmpty() const { return m_count == 0; }
		
		static uint totalCount() { return s_totalCount; }

	private:
		FaceLink * m_head;
		uint m_count;
	
		static uint s_totalCount;
	};
	
	/*static*/ uint FaceLinkList::s_totalCount = 0;
	
} // nv namespace


/// Set default kdtree build options.
void KDTreeBuilder::Options::setDefaults()
{
	leaf_face_count = 4;
	max_level = 32;
	node_traversal_cost = 1.0f;
	face_intersection_cost = 8.0f;
	max_split_faces = 0.5f;
}


/// Create a kd tree for the given faces and vertices.
KDTree * KDTreeBuilder::createTree(const Array<uint> & indices, const Array<Vector3> & vertices)
{
	KDTree * tree = new KDTree();

	nvDebug("--- Building KDTree.\n");
	
	nvDebug("---   triangle num = %u\n", indices.count() / 3);
	nvDebug("---   vertex num = %u\n", vertices.count());

	
	// Reset time counters.
	time_selectSplit = 0;
	time_selectSplitAccurate = 0;
	time_initPlanes = 0;
	time_sortPlanes = 0;
	time_sweepAndPrune = 0;
	time_splitFaces = 0;
	time_splitBounds = 0;

	// Compute bounding box.
	computeBounds(vertices);
	
	// Build faces.
	buildFaces(indices, vertices);

	uint64 t = rdtsc();
	
	const uint face_count = m_face_array.size();

	// Reserve enough memory to avoid multiple allocations.
	m_node_array.reserve(face_count * 6);
	m_leaf_face_array.reserve(face_count * 5);
	
	// Create root node pair.
	m_node_array.resize(2);    // Allocate 2 for alignment.

	FaceLinkList faceList;

	for (uint f = 0; f < face_count; f++)
	{
		faceList.append(f);
	}

	buildTree_r(0, 0, faceList, m_bounds);

	time_buildTree = rdtsc() - t;
	
	// This is true for any binary tree.
	nvCheck(stats.leaf_count == (m_node_array.size()+1)/2);

	// Display tree info.
	nvDebug("---   node num = %d\n", m_node_array.size());
	nvDebug("---   leaf num = %d\n", stats.leaf_count);
	nvDebug("---   leaf face num = %d\n", m_leaf_face_array.size());
	nvDebug("---   max leaf faces = %d\n", stats.leaf_face_max);
	nvDebug("---   empty leaf num = %d\n", stats.empty_leaf_count);
	nvDebug("---   split faces = %d\n", m_leaf_face_array.size() - m_face_array.size());
	nvDebug("---   average splits = %.2f\n", float(m_leaf_face_array.size() - m_face_array.size()) / (m_node_array.size() - stats.leaf_count));
	nvDebug("---   average leaf faces = %.2f\n", float(m_leaf_face_array.size()) / stats.leaf_count);

	nvDebug("###   time_buildTree %"POSH_I64_PRINTF_PREFIX"u\n", time_buildTree / (1024 * 1024));
	
	nvDebug("###   time_selectSplit %"POSH_I64_PRINTF_PREFIX"u\n", time_selectSplit / (1024 * 1024));
	nvDebug("###   time_selectSplitAccurate %"POSH_I64_PRINTF_PREFIX"u\n", time_selectSplitAccurate / (1024 * 1024));
	nvDebug("###   time_initPlanes %"POSH_I64_PRINTF_PREFIX"u\n", time_initPlanes / (1024 * 1024));
	nvDebug("###   time_sortPlanes %"POSH_I64_PRINTF_PREFIX"u\n", time_sortPlanes / (1024 * 1024));
	nvDebug("###   time_sweepAndPrune %"POSH_I64_PRINTF_PREFIX"u\n", time_sweepAndPrune / (1024 * 1024));
	nvDebug("###   time_splitFaces %"POSH_I64_PRINTF_PREFIX"u\n", time_splitFaces / (1024 * 1024));
	nvDebug("###   time_splitBounds %"POSH_I64_PRINTF_PREFIX"u\n", time_splitBounds / (1024 * 1024));


	m_face_array.clear();
	m_face_array.shrink();

	tree->m_bounds = m_bounds;

	// Output the tree.
	writeNodes(tree);
	m_node_array.clear();
	m_node_array.shrink();

	writeLeafFaces(tree);
	m_leaf_face_array.clear();
	m_leaf_face_array.shrink();

	tree->m_faceBuffer = new WaldFaceBuffer();
	tree->m_faceBuffer->build(indices, vertices);

	return tree;
}


/// Compute tree bounds.
void KDTreeBuilder::computeBounds(const Array<Vector3> & vertices)
{
	m_bounds.clearBounds();
	
	foreach(v, vertices) {
		m_bounds.addPointToBounds(vertices[v]);
	}
}


/// Build face array.
void KDTreeBuilder::buildFaces(const Array<uint> & indices, const Array<Vector3> & vertices)
{
	const uint index_count = indices.size();
	nvCheck(index_count % 3 == 0);
	
	const uint vertex_count = vertices.size();
	const uint face_count = index_count / 3;
	
	// Allocate tree faces.
	m_face_array.resize(face_count);

	for(uint f = 0; f < face_count; f++)
	{
		Face & face = m_face_array[f];		
		
		const uint i0 = indices[3*f + 0];
		const uint i1 = indices[3*f + 1];
		const uint i2 = indices[3*f + 2];
		
		nvCheck(i0 < vertex_count);
		nvCheck(i1 < vertex_count);
		nvCheck(i2 < vertex_count);
		
		face.winding = Triangle(vertices[i0], vertices[i1], vertices[i2]);
		face.bounds = face.winding.bounds();

		nvDebugCheck( overlap(face.winding, face.bounds) );
	}
}


/// Build tree recursively.
void KDTreeBuilder::buildTree_r(uint level, uint idx, FaceLinkList & faceList, const Box & bounds)
{
	nvCheck(idx < m_node_array.size());

	Node * node = &m_node_array[idx];
	
	if (!selectSplit(level, faceList, bounds, node))
	{
		//nvDebug("level: %d -> (%d)\n", level, node->face_array.size());

		stats.leaf_count++;

		// Set leaf flag.
		node->axis = -1;

		if (faceList.count() == 0) {
			// Empty leaf.
			node->son[0] = m_leaf_face_array.count();
			node->son[1] = 0;

			stats.empty_leaf_count++;
		}
		else {
			// Add faces to leaf face array, store first index and count in son[0] and son[1].
			node->son[0] = addLeafFaces(faceList);
			node->son[1] = faceList.count();

			faceList.deleteAll();
			
			stats.leaf_face_max = max(stats.leaf_face_max, node->son[1]);
		}
	}
	else {
		// Make sure the node is a node.
		nvDebugCheck(node->axis == 0 || node->axis == 1 || node->axis == 2);

		//nvDebug("level: %d -> %d:%f (%d)\n", level, node->axis, node->split, node->face_array.size());
		
		// Allocate two new nodes.
		const uint offset = m_node_array.size();
		node->son[0] = offset;
		node->son[1] = offset + 1;
		m_node_array.resize(offset + 2);
		
		node = NULL; // Do not use 'node' after resize!

		Box leftBounds;
		Box rightBounds;

		// Split the bounds by the node plane.
		splitBounds(m_node_array[idx], bounds, &leftBounds, &rightBounds);

		FaceLinkList leftFaceList;
		FaceLinkList rightFaceList;

		// Split the face array by the node plane.
		splitFaces(m_node_array[idx], leftBounds, rightBounds, faceList, leftFaceList, rightFaceList);

		// Recurse.
		buildTree_r(level + 1, offset, leftFaceList, leftBounds);
		buildTree_r(level + 1, offset + 1, rightFaceList, rightBounds);
	}

	// Print progress
	static int i = 0;
	if (i++ % 100 == 0) printf("\r%u", FaceLinkList::totalCount());
}


bool KDTreeBuilder::selectSplit(uint level, const FaceLinkList & faceList, const Box & bounds, Node * node)
{
	uint64 t = rdtsc();

	nvDebugCheck(node != NULL);

	const uint face_count = faceList.count();
	if (face_count <= m_options.leaf_face_count) {	// < 3
		return false;
	}
	
	if (level >= m_options.max_level) {
		//nvDebug("### Max level hit with %d faces.\n", face_count);
		return false;
	}

	bool b = selectSplitAccurate(level, faceList, bounds, node);
	//bool b = selectSplitFast(level, faceList, bounds, node);

	time_selectSplit += rdtsc() - t;

	return b;
}


/** Select best split using sweep and prune. 
 * Returns the best split in best->axis and best->split.
 * Our heuristic works as follows:
 * - Minimize surface area of the child nodes.
 * - Minimize number of split faces.
 * - Favors empty nodes when its volume is big enough.
 * The algorithm is o(N * log^2(N))
 */
bool KDTreeBuilder::selectSplitAccurate(uint level, const FaceLinkList & faceList, const Box & bounds, Node * node)
{
	uint64 t = rdtsc();

	const uint face_count = faceList.count();
	const uint max_count = uint((1.0f + m_options.max_split_faces) * face_count);
	const float total_area = bounds.area();
	const float total_area_rcp = 1.0f / total_area;
	
	// Initialize best value to the value of doing no split.
	float best_value = face_count * m_options.face_intersection_cost;

	m_plane_array.resize(2 * face_count);
	
	// Evaluate each axis.
	for (int axis = 0; axis < 3; axis++)
	{
		const float length = 2 * bounds.extents(axis);
		const float abs_min = bounds.minCorner().component(axis) + length * 0.01f;
		const float abs_max = bounds.maxCorner().component(axis) - length * 0.01f;
		
		
		if (length < NV_EPSILON) {
			// Skip this axis.
			continue;
		}
		nvDebugCheck(abs_min < abs_max);
		
		uint64 t1 = rdtsc();
		
		// Init planes.
		uint i = 0;
		for (FaceLink * link = faceList.head(); link != NULL; faceList.advance(link), i++)
		{
			const Box & faceBounds = m_face_array[link->face].bounds;
			
			m_plane_array[2 * i + 0] = faceBounds.minCorner().component(axis);
			m_plane_array[2 * i + 1] = faceBounds.maxCorner().component(axis);
		}

		time_initPlanes += rdtsc() - t1;

		uint64 t2 = rdtsc();

		// Sort the intervals.
		const uint32 * ranks = m_radix.sort(m_plane_array.buffer(), m_plane_array.count()).ranks();

		time_sortPlanes += rdtsc() - t2;

		// Make sure sort is correct!
		for (uint i = 0; i < 2*face_count-1; i++)
		{
			nvCheck(m_plane_array[ranks[i]] <= m_plane_array[ranks[i+1]]);
		}


		uint64 t3 = rdtsc();

		uint left_count = 0;
		uint right_count = face_count;

		Box left_box(bounds);
		Box right_box(bounds);

		nvDebugCheck(2 * face_count == m_plane_array.size());

		// Sweep and prune.
		for (uint i = 0; i < 2 * face_count; /*i++*/)
		{
			uint enter = 0;
			uint leave = 0;
			uint planar = 0;
			
			const float split_plane = m_plane_array[ranks[i]];
			
			// For each plane at this location.
			do {
				if (m_plane_array[ranks[i]] == m_plane_array[ranks[i]^1])
				{
					planar++;
				}
				else
				{
					// Count sides.
					if ((ranks[i] & 1) == 0) {
						enter++;
					}
					else {
						leave++;
					}
				}
				
				// Increment at least once.
				i++;
				if (i == 2 * face_count) {
					break;
				}
				
			} while (m_plane_array[ranks[i]] == split_plane);
			
			
			right_count -= leave;
			right_count += planar;
			left_count += planar;
			
			// Make sure that plane is within required bounds.
			if (split_plane > abs_min && split_plane < abs_max)
			{
				// Only consider plane if it does not split too many faces.
				if (left_count + right_count < max_count)
				{
					// Evaluate SAH
					left_box.maxCorner().setComponent(axis, split_plane);
					right_box.minCorner().setComponent(axis, split_plane);
					
					float left_area = left_box.area();
					float right_area = right_box.area();
					
					float value = m_options.node_traversal_cost + m_options.face_intersection_cost * (left_area * left_count + right_area * right_count) * total_area_rcp;
					
					if (value < best_value)
					{
						node->axis = axis;
						node->split = split_plane;
						node->son[0] = left_count;
						node->son[1] = right_count;
						best_value = value;
					}
				}
			}
			
			left_count += enter;
			right_count -= planar;
			left_count -= planar;
		}
		
		time_sweepAndPrune += rdtsc() - t3;
	}

	//uint best_splits = node->son[0] + node->son[1] - face_count;

	time_selectSplitAccurate += rdtsc() - t;

	// Evaluate (greedy) termination criterion.
	return best_value < (face_count * m_options.face_intersection_cost);
}


bool KDTreeBuilder::selectSplitFast(uint level, const FaceLinkList & faceList, const Box & bounds, Node * node)
{
	const Vector3 center = bounds.center();
	const Vector3 extents = bounds.extents();

	// Pick longest axis.
	uint axis = 0;
	if (extents.y() > extents.x()) axis = 1;
	if (extents.z() > extents.y()) axis = 2;

	float split_plane = center.component(axis);

	uint left_count = 0;
	uint right_count = 0;

	for (FaceLink * link = faceList.head(); link != NULL; faceList.advance(link))
	{
		const Box & faceBounds = m_face_array[link->face].bounds;

		if (faceBounds.minCorner().component(axis) <= split_plane) left_count++;
		if (faceBounds.maxCorner().component(axis) >= split_plane) right_count++;
	}

	node->axis = axis;
	node->split = split_plane;
	node->son[0] = left_count;
	node->son[1] = right_count;

	return true;
}


/// Add the given leaf face array to the leaf face list.
uint KDTreeBuilder::addLeafFaces(const FaceLinkList & faceList)
{
	const uint offset = m_leaf_face_array.count();
	const uint faceCount = faceList.count();

	if( faceCount == 0 ) {
		return 0;
	}
/*
	// Make sure face indices are in the right order.
	for (uint i = 1; i < faceCount; i++) {
		nvDebugCheck(face_array[i-1] < face_array[i]);
	}
*/
	for (FaceLink * link = faceList.head(); link != NULL; faceList.advance(link))
	{
		m_leaf_face_array.append(link->face);
	}
	
	// Return offset of these leaf faces.
	return offset;
}


/// Split the face links by the node plane.
void KDTreeBuilder::splitFaces(const Node & node, const Box & left_box, const Box & right_box, FaceLinkList & faceList, FaceLinkList & leftList, FaceLinkList & rightList)
{
	uint64 t = rdtsc();

	uint left_side = 0;
	uint right_side = 0;

	// Traverse face links.
	while (!faceList.isEmpty())
	{
		const uint face = faceList.head()->face;

		const Triangle & face_winding = m_face_array[face].winding;
		const Box & face_bounds = m_face_array[face].bounds;

		uint side = 0;
		
		if (node.split >= face_bounds.maxCorner().component(node.axis))
		{
			side |= 1; left_side++;
		}
		else if (node.split <= face_bounds.minCorner().component(node.axis))
		{
			side |= 2; right_side++;
		}
		else
		{
			//side = 3;
			if (overlapNoBounds(face_winding, left_box)) {
				side |= 1; left_side++;
			}
			if (overlapNoBounds(face_winding, right_box)) {
				side |= 2; right_side++;
			}
		}

		FaceLink * link = faceList.detach();
		nvDebugCheck(link != NULL);
	
		if (side == 1)
		{
			leftList.attach(link);
		}
		else if (side == 2)
		{
			rightList.attach(link);
		}
		else if (side == 3)
		{
			// @@ Clip triangle, create two new faces, etc...

			leftList.attach(link);
			rightList.append(link->face);
		}
	}

	time_splitFaces += rdtsc() - t;
}


/// Split the node bounds by its splitting plane.
void KDTreeBuilder::splitBounds(const Node & node, const Box & bounds, Box * box1, Box * box2)
{
	uint64 t = rdtsc();

	nvDebugCheck(box1 != NULL);
	nvDebugCheck(box2 != NULL);

	const float length = 2 * bounds.extents(node.axis);
	const float bias = length * 0.01f;

	// Back side node.
	*box1 = bounds;
	box1->maxCorner().setComponent(node.axis, node.split + bias);

	// Front side node.
	*box2 = bounds;
	box2->minCorner().setComponent(node.axis, node.split - bias);

	time_splitBounds += rdtsc() - t;
}


/// Get final tree size.
uint KDTreeBuilder::treeSize() const
{
	const uint node_count = m_node_array.size();
	const uint leaf_face_count = m_leaf_face_array.size();
	const uint face_count = m_face_array.size();

	const uint node_size = node_count * sizeof(KDTree::Node);
	const uint leaf_face_size = leaf_face_count * sizeof(uint);
//	const uint face_size = face_count * sizeof(KDTree::Face) + AlignOffset(node_size + leaf_face_size);
	
	const uint total_size = node_size + leaf_face_size; // + face_size;
	
	nvDebug("---   tree size = %dB = %dK.\n", total_size, total_size/1024 );
	nvDebug("---   nodes: %d%%, leaf faces: %d%%\n", (100*node_size)/total_size, (100*leaf_face_size)/total_size);
	//nvDebug("---   nodes: %d%%, leaf faces: %d%%, faces: %d%%\n", (100*node_size)/total_size, (100*leaf_face_size)/total_size, (100*face_size)/total_size);
	
	return total_size;
}


/// Write nodes to the tree memory.
void KDTreeBuilder::writeNodes(KDTree * tree) const
{
	nvCheck(tree != NULL);

	nvStaticCheck(sizeof(KDTree::Node) == 8);
	nvStaticCheck(sizeof(KDTree::Node) == sizeof(KDTree::Leaf));

	const uint nodeCount = m_node_array.count();

	tree->allocateNodes(nodeCount);

	KDTree::Node * node_ptr = tree->nodeBuffer();

	// Write compressed nodes out.
	for(uint n = 0; n < nodeCount; n++)
	{
		const Node & node = m_node_array[n];
		
		if( node.axis == -1 ) {
			// Leaf
			KDTree::Leaf * leaf_ptr = reinterpret_cast<KDTree::Leaf *>( node_ptr );

			leaf_ptr[n].flags = uint(1<<31);

			const uint offset = sizeof(uint32) * node.son[0];
			leaf_ptr[n].flags |= offset;

			leaf_ptr[n].face_count = node.son[1];
		}
		else {
			node_ptr[n].flags = node.axis;
			
			const uint offset = sizeof(KDTree::Node) * node.son[0];
			node_ptr[n].flags |= offset;
			
			node_ptr[n].distance = node.split;	
		}
	}
}


/// Write leaf faces to the tree memory.
void KDTreeBuilder::writeLeafFaces(KDTree * tree) const
{
	nvCheck(tree != NULL);

	const uint leafFaceCount = m_leaf_face_array.count();
	
	tree->allocateLeafFaces(leafFaceCount);

	uint * leaf_face_ptr = tree->leafFaceBuffer();

	// Write leaf faces.
	for (uint l = 0; l < leafFaceCount; l++) {
		leaf_face_ptr[l] = m_leaf_face_array[l];
	}
}

