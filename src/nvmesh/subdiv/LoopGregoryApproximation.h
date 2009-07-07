// Copyright NVIDIA Corporation 2007 -- Denis Kovacs <den.kovacs@gmail.com>

#ifndef NV_MESH_LOOPGREGORYAPPROXIMATION_H
#define NV_MESH_LOOPGREGORYAPPROXIMATION_H

#include <nvcore/Containers.h>
#include <nvmath/Vector.h>
#include <nvmesh/nvmesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/subdiv/Boundary.h>

namespace nv
{
	namespace HalfEdge { class Mesh; }
	class CubicGregoryMesh;
	class ACCStencils;

	/// Charles Loop subdivision surface Gregory approximation.
	class LoopGregoryApproximation
	{
	public:

		enum Mode
		{
			QuadTriangle,
			QuadsOnly,
			TrianglesOnly,
		};

		LoopGregoryApproximation(const HalfEdge::Mesh * mesh, BoundaryMode boundaryMode = BoundaryMode_None);
		~LoopGregoryApproximation();

		CubicGregoryMesh * buildMesh(Mode mode = QuadTriangle);

	private:

		HalfEdge::Mesh * m_mesh;
		BoundaryMode m_boundaryMode;
	};

} // nv namespace

#endif // NV_MESH_LOOPGREGORYAPPROXIMATION_H
