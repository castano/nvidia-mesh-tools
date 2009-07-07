// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_BAKER_TRIBASESURFACE_H
#define NV_BAKER_TRIBASESURFACE_H

#include <nvcore/Ptr.h>

#include "BaseSurface.h"

namespace nv
{
	class Basis;
	class TriMesh;
	namespace HalfEdge { class Mesh; }
	
	enum PositionMode {
		PositionMode_Linear,
		PositionMode_PNTriangle,
		PositionMode_NMTriangle,
	};

	enum NormalMode {
		NormalMode_Linear,
		NormalMode_Quadratic,
		NormalMode_ImprovedQuadratic,
	};

	enum TangentMode {
		TangentMode_Linear,
		TangentMode_LinearOrtho,	// Recompute bitangent from normal/
	};


	enum TangentSpaceMode {
		TangentSpaceMode_MeshMender,
		TangentSpaceMode_Lengyel,
		TangentSpaceMode_Castano,
		TangentSpaceMode_UE3,
	};

	enum TangentSpaceProcessing {
		TangentSpaceProcessing_None,
		TangentSpaceProcessing_GramShmidtOrthogonalize,
		TangentSpaceProcessing_CastanoOrthogonalize,
		TangentSpaceProcessing_UE3,
	};

	
	class TriBaseSurface : public BaseSurface
	{
	public:
		TriBaseSurface(TriMesh * mesh);
		~TriBaseSurface();

		void prepare(TangentSpaceMode tangentMode);

		void setInterpolationMode(PositionMode positionMode, NormalMode normalMode);
		
		virtual uint faceCount() const;
		virtual void selectFace(uint i);
		
		virtual FaceDomain::Enum domain() const;
		virtual void textureCoordinates(Vector2 * texCoordArray) const;
		
		virtual void evaluate(float u, float v, Vector3 * pos, Basis * patchFrame, Basis * chartFrame) const;
	
	private:
		AutoPtr<TriMesh> m_mesh;
		uint m_faceIndex;

		Array<Basis> m_basisArray;

		PositionMode m_positionMode;
		NormalMode m_normalMode;

	};
	
};

#endif // NV_BAKER_TRIBASESURFACE_H
