// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_BAKER_BASESURFACE_H
#define NV_BAKER_BASESURFACE_H

#include <nvmath/Vector.h>

namespace nv
{
	class Basis;
	
	namespace FaceDomain
	{
		enum Enum {
			Triangle,
			Quad
		};
	}

	// BaseSurface interface.
	class BaseSurface
	{
	public:
		virtual ~BaseSurface() {};

		virtual uint faceCount() const = 0;
		virtual void selectFace(uint i) = 0;
		
		virtual FaceDomain::Enum domain() const = 0;
		virtual void textureCoordinates(Vector2 * texCoordArray) const = 0;

		virtual void evaluate(float u, float v, Vector3 * pos, Basis * patchFrame, Basis * chartFrame) const = 0;
	};
	
};

#endif // NV_BAKER_BASESURFACE_H
