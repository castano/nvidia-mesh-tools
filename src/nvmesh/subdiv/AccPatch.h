// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_ACCPATCH_H
#define NV_MESH_ACCPATCH_H

#include <nvcore/Containers.h>
#include <nvcore/RefCounted.h>

#include <nvmath/Vector.h>

#include <nvmesh/nvmesh.h>
#include <nvmesh/subdiv/StencilMask.h>

namespace nv
{
	namespace HalfEdge { class Vertex; class Edge; class Face; }
	class Basis;
	class AccMesh;
	struct FaceTopology;


	enum PatchType
	{
		PatchType_Bezier,
		PatchType_Polar,
		PatchType_Gregory,
		PatchType_PmQuad,
		PatchType_PmTriangle,
        PatchType_PmRegular,
		PatchType_TexCoord,
		PatchType_TangentSpace,
	};

	class FaceOneRing
	{
		NV_FORBID_COPY(FaceOneRing);
	public:

		FaceOneRing(const HalfEdge::Face * face, const HalfEdge::Edge * edge, FaceTopology * faceTopology);

		const HalfEdge::Face * face() const { return m_face; }
		const HalfEdge::Edge * firstEdge() const { return m_firstEdge; }
		FaceTopology * faceTopology() const { return m_faceTopology; }

		uint vertexCount() const;
        const HalfEdge::Vertex * vertexAt(uint i) const;
		uint vertexIndexAt(uint i) const;
		uint vertexIndex(const HalfEdge::Vertex * vertex) const;

		Vector3 evalutePositionStencil(const StencilMask & mask) const;
        Vector3 evaluteNormalStencil(const StencilMask & mask) const;

		//bool isRegular() const;
		bool isTriangle() const;
		bool isQuad() const;
		uint edgeCount() const;

		static bool isRegular(const HalfEdge::Face * face);
		static bool isPolarTriangle(const HalfEdge::Face * face);

		static bool isTriangle(const HalfEdge::Face * face);
		static bool isQuad(const HalfEdge::Face * face);
		static bool isPent(const HalfEdge::Face * face);

		static bool isBoundary(const HalfEdge::Face * face);

	protected:

		void initVertices();
		void appendVertex(const HalfEdge::Vertex * vertex);
		bool hasVertex(const HalfEdge::Vertex * vertex) const;

	protected:

		const HalfEdge::Face * const m_face;
		const HalfEdge::Edge * const m_firstEdge;
		FaceTopology * m_faceTopology;

		uint m_edgeCount;

		Array<const HalfEdge::Vertex *> m_vertexArray;
		Array<uint> m_vertexIndexArray;

	};

	class AccPatch
	{
		NV_FORBID_COPY(AccPatch);
	public:

		AccPatch(const FaceOneRing * face);
		virtual ~AccPatch();

		const HalfEdge::Face * face() const { return m_face->face(); }
		const FaceOneRing * faceOneRing() const { return m_face; }
		
		uint vertexCount() const { return m_face->vertexCount(); }
        const HalfEdge::Vertex * vertexAt(uint i) const { return m_face->vertexAt(i); }
		uint vertexIndex(const HalfEdge::Vertex * vertex) const { return m_face->vertexIndex(vertex); }

		//bool isRegular() const { return m_face->isRegular(); }
		bool isTriangle() const { return m_face->isTriangle(); }
		bool isQuad() const { return m_face->isQuad(); }

		virtual PatchType patchType() const = 0;

	protected:

		const FaceOneRing * const m_face;

	};

	
	class SurfacePatch : public AccPatch
	{
		NV_FORBID_COPY(SurfacePatch);
	public:

		SurfacePatch(const FaceOneRing * face);

		virtual void evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const = 0;
	};


	// 
	class TexCoordPatch : public AccPatch
	{
		NV_FORBID_COPY(TexCoordPatch);
	public:
		
		TexCoordPatch(const FaceOneRing * face);

		virtual PatchType patchType() const { return PatchType_TexCoord; }

		const Vector2 texCoord(uint v, uint i) const; // v ~ [0:4), i ~ [0:4)
		Vector2 & texCoord(uint v, uint i);

	protected:
		Vector2 m_texCoordArray[16];
	};


	// 
	class TangentSpacePatch : public AccPatch
	{
		NV_FORBID_COPY(TangentSpacePatch);
	public:
		
		TangentSpacePatch(const FaceOneRing * face);

		virtual PatchType patchType() const { return PatchType_TangentSpace; }

		Vector3 chartTangent(uint idx) const;
		Vector3 chartBitangent(uint idx) const;

		Vector3 evaluateChartTangent(float u, float v) const;
		Vector3 evaluateChartBitangent(float u, float v) const;

		virtual void evaluateChartFrame(float u, float v, Basis * frame) const;

	private:

		Vector2 m_chartTangentArray[4];
		Vector2 m_chartBitangentArray[4];

	};


	//
	class BezierAccStencil
	{
		NV_FORBID_COPY(BezierAccStencil);
	public:
		BezierAccStencil(const FaceOneRing * face);

		const StencilMask & positionStencil(uint i) const;
		StencilMask & positionStencil(uint i);
		float & positionStencil(uint i, const HalfEdge::Vertex * vertex);
		
		const StencilMask & tangentStencil(uint i) const;
		StencilMask & tangentStencil(uint i);
		float & tangentStencil(uint i, const HalfEdge::Vertex * vertex);
		
		const StencilMask & bitangentStencil(uint i) const;
		StencilMask & bitangentStencil(uint i);
		float & bitangentStencil(uint i, const HalfEdge::Vertex * vertex);

		bool operator==(const BezierAccStencil & other) const;

	private:
		const FaceOneRing * m_face;

		StencilMask m_positionStencil[16];
		StencilMask m_tangentStencil[12];
		StencilMask m_bitangentStencil[12];
	};



	// Bicubic Bezier Approximation to Catmull Clark.
	class BezierAccPatch : public SurfacePatch
	{
		NV_FORBID_COPY(BezierAccPatch);
	public:
		
		BezierAccPatch(const FaceOneRing * face);
		virtual ~BezierAccPatch();

		virtual PatchType patchType() const { return PatchType_Bezier; }

		const BezierAccStencil * stencil() const { return m_stencil; }
		BezierAccStencil * stencil() { return m_stencil; }

		void evaluateControlPoints();
	
		Vector3 position(uint idx) const;   // [0:16)
		Vector3 tangent(uint idx) const;    // [0:12)
		Vector3 bitangent(uint idx) const;  // [0:12)
        Vector3 normal(uint idx) const;     // [0:16)

		//Vector3 normal(uint idx) const;     // [0:4)

		Vector3 evaluateTangent(float u, float v) const;
		Vector3 evaluateBitangent(float u, float v) const;

		virtual void evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const;
		void evaluatePosition(float u, float v, Vector3 * pos) const;
		void evaluatePatchFrame(float u, float v, Basis * frame) const;

	private:
	
		Vector3 p(uint u, uint v) const;
		Vector3 tu(uint u, uint v) const;
		Vector3 tv(uint u, uint v) const;
	
	public:
		
		BezierAccStencil * m_stencil;
		
		Vector3 m_positionArray[16];
		Vector3 m_tangentArray[12];
		Vector3 m_bitangentArray[12];
        Vector3 m_normalArray[16]; // for subdivision shading

	};


	//
	class GregoryAccStencil
	{
		NV_FORBID_COPY(GregoryAccStencil);
	public:
		GregoryAccStencil(const FaceOneRing * m_face);

		const StencilMask & positionStencil(uint i) const;
		StencilMask & positionStencil(uint i);
		float & positionStencil(uint i, const HalfEdge::Vertex * vertex);
	
		bool operator==(const GregoryAccStencil & other) const;

	private:
		const FaceOneRing * m_face;

		StencilMask m_positionStencil[20];
	};

	// Approximation to Catmull Clark using Gregory patches.
	class  GregoryAccPatch : public SurfacePatch
	{
		NV_FORBID_COPY(GregoryAccPatch);
	public:
		
		GregoryAccPatch(const FaceOneRing * face);
		virtual ~GregoryAccPatch();

		virtual PatchType patchType() const { return PatchType_Gregory; }

		const GregoryAccStencil * stencil() const { return m_stencil; }
		GregoryAccStencil * stencil() { return m_stencil; }
		
		void evaluateControlPoints();
	
		Vector3 position(uint idx) const;   // [0:20)
		
		virtual void evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const;

	public:

		// quad:
		//  8     9     10     11
		// 12   0\1     2/3    13
		// 14   4/5     6\7    15
		// 16    17     18     19

		// triangle:
		//              6
        //              
		//        14   0/1   7
		//                          
		//    13   5/4     3\2   8
		//               
		// 12      11       10      9

		GregoryAccStencil * m_stencil;
		
		// These can be derived by the stencils:
		Vector3 m_positionArray[20];
		
	};


	// Jorg Peters polar subdivision approximation.
	class PolarPatch : public SurfacePatch
	{
		NV_FORBID_COPY(PolarPatch);
	public:

		PolarPatch(const FaceOneRing * face);

		virtual PatchType patchType() const { return PatchType_Polar; }

		virtual void evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const;

	private:

		//        0
		//     1 2 3 4
		//   5  6   7  8
		// 9   10   11  12

		Vector3 m_positionArray[13];

	};


	// Approximation to Catmull Clark using Pm patches.
	class PmQuadAccPatch : public SurfacePatch
	{
		NV_FORBID_COPY(PmQuadAccPatch);
	public:
		
		PmQuadAccPatch(const FaceOneRing * face);

		virtual PatchType patchType() const { return PatchType_PmQuad; }

		Vector3 position(uint idx) const;   // [0:24)

		virtual void evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const;

	private:
	
		Vector3 evaluateTangent(float u, float v) const;
		Vector3 evaluateBitangent(float u, float v) const;
        void evaluatePosition(float u, float v, Vector3 * pos) const;
		void evaluatePatchFrame(float u, float v, Basis * frame) const;

	public:

		Vector3 m_positionArray[24];
		
	};


	class PmTriangleAccPatch : public SurfacePatch
	{
		NV_FORBID_COPY(PmTriangleAccPatch);
	public:
		
		PmTriangleAccPatch(const FaceOneRing * face);

		virtual PatchType patchType() const { return PatchType_PmTriangle; }

		Vector3 position(uint idx) const;   // [0:19)

		virtual void evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const;

	private:
	
		Vector3 evaluateTangent(float u, float v) const;
		Vector3 evaluateBitangent(float u, float v) const;
        void evaluatePosition(float u, float v, Vector3 * pos) const;
		void evaluatePatchFrame(float u, float v, Basis * frame) const;

	public:

		Vector3 m_positionArray[19];
		
	};


    // @@ This should be merged with BezierAccPatch. 
    // However, control points are currenly computed in slightly different ways, and would result in slightly 
    // different values.
	class PmRegularAccPatch : public SurfacePatch
	{
		NV_FORBID_COPY(PmRegularAccPatch);
	public:
		
		PmRegularAccPatch(const FaceOneRing * face);

		virtual PatchType patchType() const { return PatchType_PmRegular; }

		Vector3 position(uint idx) const;   // [0:16)

		virtual void evaluateSurface(float u, float v, Vector3 * pos, Vector3 * du, Vector3 * dv) const;

	private:
	
		Vector3 evaluateTangent(float u, float v) const;
		Vector3 evaluateBitangent(float u, float v) const;
        void evaluatePosition(float u, float v, Vector3 * pos) const;
		void evaluatePatchFrame(float u, float v, Basis * frame) const;


	public:

		Vector3 m_positionArray[16];
		
	};

} // nv namespace

#endif // NV_MESH_ACCPATCH_H
