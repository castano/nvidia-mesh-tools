// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_ACCMESH_H
#define NV_MESH_ACCMESH_H

#include <nvcore/Containers.h>
#include <nvcore/RefCounted.h>

#include <nvmesh/nvmesh.h>

namespace nv
{
	namespace HalfEdge { class Vertex; class Face; class Edge; class Mesh; }
	class FaceOneRing;
	class AccPatch;
	class SurfacePatch;
	class TexCoordPatch;
	class TangentSpacePatch;
	class BezierAccPatch;
	class GregoryAccPatch;

	class BezierAccStencil;
	class GregoryAccStencil;


	NV_DECLARE_PTR(AccMesh);

	// Structure to describe the topology of a vertex. Keep it small.
	struct VertexTopology
	{
		VertexTopology();
		VertexTopology(const VertexTopology & vt);
		void operator=(const VertexTopology & vt);
		bool operator==(const VertexTopology & vt) const;
		VertexTopology(const HalfEdge::Vertex * vertex, const HalfEdge::Edge * firstEdge);

		uint valence : 4;		// 0-15
		uint isBoundary : 1;	// true/false
		uint isRegular : 1;		// true/false
		uint faceBitField;		// two bits per face (0=quad, 1=triangle, 3=boundary) 
	};

	// Structure to describe the toplogy of a face. Also stores the stencils for this topology.
	struct FaceTopology
	{
		uint index;
		uint topologyId;
		uint edgeCount;
		bool isRegular;

		BezierAccStencil * bezierAccStencil;
		GregoryAccStencil * gregoryAccStencil;
	};


	// An ACC mesh is a set of ACC patches.
	class AccMesh : public RefCounted
	{
		NV_FORBID_COPY(AccMesh);
	public:

		enum FaceType
		{
			FaceType_Undefined,
			FaceType_Regular,
			FaceType_Quad,
			FaceType_Tri,
		};

		struct Patch
		{
			Patch() : faceType(FaceType_Undefined), faceOneRing(NULL), bezierAccPatch(NULL), gregoryAccPatch(NULL), pmAccPatch(NULL), texCoordPatch(NULL), tangentSpacePatch(NULL) {}

			FaceType faceType;
			FaceOneRing * faceOneRing;

			BezierAccPatch * bezierAccPatch;
			GregoryAccPatch * gregoryAccPatch;
			SurfacePatch * pmAccPatch; // regular, quad or tri
			TexCoordPatch * texCoordPatch;
			TangentSpacePatch * tangentSpacePatch;
		};


		AccMesh();
		virtual ~AccMesh();

		// Add vertex, return topology index.
		FaceTopology * addFaceTopology(const HalfEdge::Face * face, const HalfEdge::Edge ** firstEdge);
		uint faceTopologyCount() const;
		const FaceTopology * faceTopologyAt(uint i) const;

		uint addVertexTopology(const HalfEdge::Vertex * vertex, const HalfEdge::Edge * firstEdge);
		uint vertexTopologyCount() const;

		uint addRegularPatch(FaceOneRing * face);
		uint addQuadPatch(FaceOneRing * face);
		uint addTriPatch(FaceOneRing * face);

		void setRegularPatch(uint p, AccPatch * patch);
		void setQuadPatch(uint p, AccPatch * patch);
		void setTriPatch(uint p, AccPatch * patch);

		uint regularPatchCount() const { return m_regularPatchArray.count(); }
		uint quadPatchCount() const { return m_quadPatchArray.count(); }
		uint triPatchCount() const { return m_triPatchArray.count(); }

		const Patch & regularPatchAt(uint i) const { return m_regularPatchArray[i]; }
		const Patch & quadPatchAt(uint i) const { return m_quadPatchArray[i]; }
		const Patch & triPatchAt(uint i) const { return m_triPatchArray[i]; }

		uint patchCount() const;
		const Patch & patchAt(uint i) const;

		enum Sort
		{
			Sort_ByVertexCount = 0x01,
			Sort_ByTopology = 0x02,
			Sort_ByFaceId = 0x04,
		};

		void sort(const Sort * order, uint orderCount);

		enum Granularity
		{
			Granularity_Single,
			Granularity_PerTopology,
			Granularity_PerVertexCount,
		};

		void optimize(Granularity granularity = Granularity_Single);

        HalfEdge::Mesh * tessellate(uint level);

	private:

		uint evaluateTopologyId(const HalfEdge::Face * face, const HalfEdge::Edge * edge);

		void sort(const Sort * order, uint orderCount, Array<Patch> & patchArray);
		void optimize(Granularity granularity/*= Granularity_Single*/, Array<Patch> & patchArray);

	private:

		Array<Patch> m_regularPatchArray;
		Array<Patch> m_quadPatchArray;
		Array<Patch> m_triPatchArray;

		Array<VertexTopology> m_vertexTopologyArray;
		HashMap<VertexTopology, uint> m_vertexTopologyMap;

		Array<FaceTopology *> m_faceTopologyArray;
		HashMap<uint, uint> m_faceTopologyMap;

	};

} // nv namespace

#endif // NV_MESH_ACCMESH_H
