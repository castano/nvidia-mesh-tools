// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_BAKER_ACCBASESURFACE_H
#define NV_BAKER_ACCBASESURFACE_H

#include "BaseSurface.h"

#include <nvcore/Ptr.h>


namespace nv
{
	class AccMesh;
	class SurfacePatch;
	class TexCoordPatch;
	class TangentSpacePatch;
	
	class AccBaseSurface : public BaseSurface
	{
	public:
		AccBaseSurface(const AccMesh * mesh);
		~AccBaseSurface();
		
		virtual uint faceCount() const;
		virtual void selectFace(uint i);
		
		virtual FaceDomain::Enum domain() const;
		virtual void textureCoordinates(Vector2 * texCoordArray) const;
		
		virtual void evaluate(float u, float v, Vector3 * pos, Basis * patchFrame, Basis * chartFrame) const;
	
	private:
		SmartPtr<const AccMesh> m_mesh;
		const SurfacePatch * m_surfacePatch;
		const TexCoordPatch * m_texCoordPatch;
		const TangentSpacePatch * m_tangentSpacePatch;
	};
	
};

#endif // NV_BAKER_ACCBASESURFACE_H
