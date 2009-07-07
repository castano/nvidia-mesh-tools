// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_BOUNDARY_H
#define NV_MESH_BOUNDARY_H

namespace nv
{

	// Subdivision surface boundary types.
	enum BoundaryMode
	{
		BoundaryMode_None,
		BoundaryMode_Spline
		// @@ Add spline mode with sharp corners.
	};

} // namespace nv

#endif // NV_MESH_BOUNDARY_H

