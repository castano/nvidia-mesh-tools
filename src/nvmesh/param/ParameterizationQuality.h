// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_PARAMETERIZATIONQUALITY_H
#define NV_MESH_PARAMETERIZATIONQUALITY_H

#include <nvmesh/nvmesh.h>

namespace nv
{
	class Vector2;
	class Vector3;

	namespace HalfEdge
	{
		class HalfEdgeMesh;
	}
	
	// Estimate quality of existing parameterization.
	NVMESH_CLASS class ParameterizationQuality
	{
	public:
		ParameterizationQuality(const HalfEdge::Mesh * mesh);

		bool isValid() const;

		float rmsStretchMetric() const;
		float maxStretchMetric() const;
	
		float rmsConformalMetric() const;
		float maxAuthalicMetric() const;

	private:
	
		void processTriangle(Vector3 p[3], Vector2 t[3]);
	
	private:

		uint m_totalTriangleCount;
		uint m_flippedTriangleCount;
		
		float m_parametricArea;
		float m_geometricArea;
		
		float m_rmsStretchMetric;
		float m_maxStretchMetric;

		float m_rmsConformalMetric;
		float m_rmsAuthalicMetric;

	};
	
} // nv namespace

#endif // NV_MESH_PARAMETERIZATIONQUALITY_H
