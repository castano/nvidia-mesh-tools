// This code is in the public domain -- castanyo@yahoo.es

#include "HalfEdge.h"
#include "HalfEdgeVertex.h"

using namespace nv;

Vector3 HalfEdge::Edge::midPoint() const
{
	return (to()->pos() + from()->pos()) * 0.5f;
}

float HalfEdge::Edge::length() const
{
	return ::length(to()->pos() - from()->pos()); 
}

