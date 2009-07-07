// This code is in the public domain -- castanyo@yahoo.es

#include <nvmath/Fitting.h>

#include "HalfEdgeFace.h"
#include "HalfEdgeVertex.h"


using namespace nv;
using namespace HalfEdge;

/// Get face area.
float Face::area() const
{
	float area = 0;
	const Vector3 & v0 = m_edge->from()->pos();
	
	for(ConstEdgeIterator it(edges(m_edge->next())); it.current() != m_edge; it.advance())
	{
		const Edge * edge = it.current();
		
		const Vector3 & v1 = edge->from()->pos();
		const Vector3 & v2 = edge->to()->pos(); 
		
		area += length(cross(v1-v0, v2-v0));
	}
	
	return area * 0.5f;
}


/// Get face normal.
Vector3 Face::normal() const
{
	// Get face points eliminating duplicates.
	Array<Vector3> points(4);
	
	points.append(m_edge->prev()->from()->pos());
	
	for(ConstEdgeIterator it(edges()); !it.isDone(); it.advance())
	{
		const Edge * edge = it.current();
		nvDebugCheck(edge != NULL);
		
		const Vector3 & p = edge->from()->pos();
		if (points.back() != p)
		{
			points.append(edge->from()->pos());
		}
	}
	
	points.popBack();
	
	if (points.count() < 3)
	{
		// Invalid normal.
		return Vector3(zero);
	}
	else
	{
		// Compute regular normal.
		Vector3 normal = normalizeSafe(cross(points[1] - points[0], points[2] - points[0]), Vector3(zero), 0.0f);
		
#pragma message(NV_FILE_LINE "TODO: make sure these three points are not colinear")

		if (points.count() > 3)
		{
			// Compute best fitting plane to the points.
			Plane plane = Fit::bestPlane(points.count(), points.buffer());
			
			// Adjust normal orientation.
			if (dot(normal, plane.vector()) > 0) {
				normal = plane.vector();
			}
			else {
				normal = -plane.vector();
			}
		}
		
		nvDebugCheck(isNormalized(normal));
		return normal;
	}
}

Vector3 Face::trg_centroid() const
{
	nvCheck(m_edge->next()->next()->next() == m_edge); // check that this is triangle
	return (m_edge->from()->pos() + m_edge->next()->from()->pos() + m_edge->next()->next()->from()->pos()) / 3;
}
