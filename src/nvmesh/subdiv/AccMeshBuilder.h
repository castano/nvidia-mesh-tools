// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_ACCMESHBUILDER_H
#define NV_MESH_ACCMESHBUILDER_H

#include <nvcore/Containers.h>

#include <nvmath/Basis.h>

#include <nvmesh/nvmesh.h>
#include <nvmesh/subdiv/Boundary.h>

namespace nv
{
	namespace HalfEdge { class Mesh; class Face; class Vertex; }
	class AccMesh;
	class FaceOneRing;
    class AccPatch;
	class TexCoordPatch;
	class BezierAccPatch;
	class GregoryAccPatch;
    class PolarPatch;
	class PmQuadAccPatch;
	class PmRegularAccPatch;
	class PmTriangleAccPatch;
	class StencilMask;
	class Basis;

	enum BuildFlags
	{
		BuildFlags_BuildBezierPatches = 0x0001,
		BuildFlags_BuildGregoryPatches = 0x0002,
        BuildFlags_BuildPmPatches = 0x0004,
        BuildFlags_BuildPolarPatches = 0x0008,
		BuildFlags_BuildTexCoordPatches = 0x0010,
		BuildFlags_BuildTangentSpacePatches = 0x0020,
		BuildFlags_BuildAll = BuildFlags_BuildBezierPatches | BuildFlags_BuildGregoryPatches | BuildFlags_BuildPmPatches | BuildFlags_BuildTexCoordPatches,

		BuildFlags_GenerateRegularPatches = 0x0100,
        BuildFlags_GenerateQuadPatches = 0x0200,
        BuildFlags_GenerateTriangularPatches = 0x0400,
        BuildFlags_GenerateAll = BuildFlags_GenerateRegularPatches | BuildFlags_GenerateQuadPatches | BuildFlags_GenerateTriangularPatches,

        BuildFlags_Default = BuildFlags_BuildGregoryPatches | BuildFlags_GenerateAll,
	};


	/// ACC mesh builder that uses Loop's approximation.
	class AccMeshBuilder
	{
	public:
		AccMeshBuilder(const HalfEdge::Mesh * mesh, BoundaryMode bm = BoundaryMode_Spline);
		~AccMeshBuilder();

		AccMesh * buildAccMesh(uint flags = BuildFlags_Default) const;

	private:

		// Bezier patch
		BezierAccPatch * buildBezierAccPatch(const FaceOneRing * face) const;
		
		void computePositionStencil(BezierAccPatch * patch) const;
		void computeCornerPositionStencil(BezierAccPatch * patch) const;
		void computeEdgePositionStencil(BezierAccPatch * patch) const;
		void computeInteriorPositionStencil(BezierAccPatch * patch) const;
		
		void computeTangentStencil(BezierAccPatch * patch) const;
		void computeCornerTangentStencil(BezierAccPatch * patch) const;
		void computeEdgeTangentStencil(BezierAccPatch * patch) const;
		void computeInteriorTangentStencil(BezierAccPatch * patch) const;

		// Gregory patch
		GregoryAccPatch * buildGregoryAccPatch(const FaceOneRing * face) const;
		
		void computeCornerStencil(GregoryAccPatch * patch) const;
		void computeEdgeStencil(GregoryAccPatch * patch) const;
		void computeInteriorStencil(GregoryAccPatch * patch) const;

        // Polar patch
		PolarPatch * buildPolarPatch(const FaceOneRing * face) const;

		// Pm quad patch
		PmQuadAccPatch * buildPmQuadAccPatch(const FaceOneRing * face) const;
		void computeCoefficients(PmQuadAccPatch * patch) const;
		
		// Pm regular patch
	    PmRegularAccPatch * buildPmRegularAccPatch(const FaceOneRing * face) const;
		void convertToBicubic(PmRegularAccPatch * patch) const;

		// Pm triangle patch
	    PmTriangleAccPatch * buildPmTriangleAccPatch(const FaceOneRing * face) const;
		void convertToP3patch(PmTriangleAccPatch * patch) const;

		// TexCoord patch
		TexCoordPatch * buildTexCoordPatch(const FaceOneRing * face) const;
		void computeTexcoords(TexCoordPatch * patch) const;

		// TangentSpace patch
		void computeChartTangents(BezierAccPatch * patch, const Array<Basis> & basisArray) const;


		// Shared methods
		void computeGlobalTangentSpace(Array<Basis> & basisArray) const;
		void computeBoundaryTangentStencils(const AccPatch * patch, const HalfEdge::Vertex * vertex, StencilMask & r0, StencilMask & r1) const;

	private:

		const HalfEdge::Mesh * m_mesh;
		BoundaryMode m_boundaryMode;

	};

} // nv namespace

#endif // NV_MESH_ACCMESHBUILDER_H
