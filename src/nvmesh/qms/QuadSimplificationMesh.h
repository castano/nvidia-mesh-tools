// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_QUADSIMPLIFICATIONMESH_H
#define NV_MESH_QUADSIMPLIFICATIONMESH_H

#include <nvcore/Containers.h>

#include <nvmath/Vector.h>
#include <nvmesh/nvmesh.h>
#include <nvmesh/BaseMesh.h>

namespace nv
{
	namespace HalfEdge { class Mesh; class Edge; class Face; class Vertex; }

	namespace QuadMeshSimplification
	{
		// Mesh simplification:
		void simplify(HalfEdge::Mesh * mesh, uint maxQuad, float maxQEM, float v = 0.9f, float q = 0.05f, float d = 0.05f);

		// Mesh smoothing:
		void smooth(HalfEdge::Mesh * mesh);

		// Simplification operators:
		void polyChordCollapse(HalfEdge::Edge * edge);
		void quadCollapse(HalfEdge::Face * face);
		void doubletCollapse(HalfEdge::Vertex * vertex);

		// Helpers:
		bool isDoublet(const HalfEdge::Face * face0, const HalfEdge::Face * face1);

		struct SimplificationState;

		// Functions for interactive simplification:
		SimplificationState * simplifyStart(HalfEdge::Mesh * mesh, uint maxQuad, float maxQEM, float v = 0.9f, float q = 0.05f, float d = 0.05f);
		bool simplifyStep(SimplificationState * state);
		void simplifyEnd(SimplificationState * state);
	}

} // nv namespace

#endif // NV_MESH_QUADSIMPLIFICATIONMESH_H
