// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_SAMPLERS_H
#define NV_SAMPLERS_H

#include <nvmath/Basis.h>

#include <nvimage/FloatImage.h>
#include <nvimage/HoleFilling.h>	// BitMap

#include "GeometryImage.h"

namespace nv
{
	class BaseSurface;
	class KDTree;

	/// Sampler that captures positions and normals.
	class GeometrySampler
	{
	public:
		GeometrySampler(GeometryImage & image, BitMap & imageMask);

		static void sampleTriCallback(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);
		static void sampleQuadCallback(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);

		void setCurrentFace(uint vertexCount, const Vector3 * positions, const Vector3 * normals);

		void sampleTri(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);
		void sampleQuad(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);

	private:

		GeometryImage & m_image;
		BitMap & m_imageMask;
		uint m_currentVertexCount;
		const Vector3 * m_normals;
		const Vector3 * m_positions;
		Vector3 m_midedgenormals[5];
	};


	/// Sampler that generates vector displacements from bezier patches given a geometry image.
	class DisplacementPatchSampler
	{
	public:

		DisplacementPatchSampler(const GeometryImage & detailedGeometryMap, const BitMap & detailedGeometryMask, GeometryImage & geometryMap, BitMap & geometryMask);

		static void sampleCallback(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);
		
		// setDetailedGeometryMap(const GeometryImage & detailedGeometryMap, const BitMap & detailedGeometryMask);
		// setOutputGeometryMap(GeometryImage & detailedGeometryMap, BitMap & detailedGeometryMask);
		
		void setTangentSpace(bool enabled);
		void setVectorDisplacement(bool enabled);
		void setCurrentSurface(const BaseSurface * surface);
		
		void sample(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);

	private:

		const GeometryImage & m_detailedGeometryMap;
		const BitMap & m_detailedGeometryMask;
		GeometryImage & m_geometryMap;
		BitMap & m_geometryMask;
		
		bool m_tangentSpace;
		bool m_vectorDisplacement;

		const BaseSurface * m_surface;

	};

} // nv namespace

#endif // NV_SAMPLERS_H
