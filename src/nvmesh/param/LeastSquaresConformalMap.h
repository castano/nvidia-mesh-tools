// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_LEASTSQUARESCONFORMALMAP_H
#define NV_MESH_LEASTSQUARESCONFORMALMAP_H

namespace nv
{
	namespace HalfEdge { class Mesh; }
	
	void computeLeastSquaresConformalMap(HalfEdge::Mesh * mesh);
	
} // nv namespace

#endif // NV_MESH_LEASTSQUARESCONFORMALMAP_H
