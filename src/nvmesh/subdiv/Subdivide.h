// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_SUBDIVIDE_H
#define NV_MESH_SUBDIVIDE_H

#include <nvmath/Vector.h>

namespace nv
{
	namespace HalfEdge { class Mesh; }
	class FloatImage;

	namespace Subdivide
	{

		HalfEdge::Mesh * doCatmullClarkSplit(const HalfEdge::Mesh * mesh);
		void displace(HalfEdge::Mesh * mesh, const FloatImage * img, int channel=0, Vector2::Arg offset=Vector2(zero));

		HalfEdge::Mesh * doQuadTriangleSplit(const HalfEdge::Mesh * mesh);

		HalfEdge::Mesh * reduceToQuadOnly(const HalfEdge::Mesh * mesh);

		void projectToCatmullClarkLimitSurface(HalfEdge::Mesh * mesh);
		void projectToLoopLimitSurface(HalfEdge::Mesh * mesh);

	} // Subdivide namespace

} // nv namespace

#endif // NV_MESH_SUBDIVIDE_H
