// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmath/Montecarlo.h>

#include <nvmesh/raytracing/Raytracing.h>
#include <nvmesh/kdtree/KDTree.h>

#include "Samplers.h"
#include "BaseSurface.h"

using namespace nv;



GeometrySampler::GeometrySampler(GeometryImage & image, BitMap & imageMask) : m_image(image), m_imageMask(imageMask)
{
	nvCheck(image.width() == imageMask.width());
	nvCheck(image.height() == imageMask.height());
}

/*static*/ 
void GeometrySampler::sampleTriCallback(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	((GeometrySampler *)param)->sampleTri(x, y, bar, dx, dy, coverage);
}

/*static*/
void GeometrySampler::sampleQuadCallback(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	((GeometrySampler *)param)->sampleQuad(x, y, bar, dx, dy, coverage);
}


void GeometrySampler::setCurrentFace(uint vertexCount, const Vector3 * positions, const Vector3 * normals)
{
	nvDebugCheck(vertexCount<=4);

	m_positions = positions;
	m_normals = normals;

	for (uint k=0; k<vertexCount; k++) 
	{
		m_midedgenormals[k] = normalizeSafe(normals[k] + normals[(k+1)%vertexCount], Vector3(zero), 0);
	}

	if (vertexCount==4) {
		m_midedgenormals[4] = normalizeSafe(m_normals[0]+normals[1]+normals[2]+normals[3], Vector3(zero), 0);
	}
}



void GeometrySampler::sampleTri(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	nvDebugCheck(isFinite(coverage));

	Vector3 position = bar.x() * m_positions[0] + bar.y() * m_positions[1] + bar.z() * m_positions[2];

	m_image.addPixel(coverage * position, x, y, m_image.positionChannel());
	
	//Vector3 normal = normalizeSafe(cross(m_positions[1] - m_positions[0], m_positions[2] - m_positions[0]), Vector3(zero), 0.0f);
	Vector3 linearNormal = normalizeSafe(bar.x() * m_normals[0] + bar.y() * m_normals[1] + bar.z() * m_normals[2], Vector3(zero), 0.0f);
	float u=bar.x(), v=bar.y(), w=bar.z();

	Vector3 normal = normalizeSafe(	u*u*m_normals[0] + v*v*m_normals[1] + w*w*m_normals[2] + 
		2*(u*v*m_midedgenormals[0] + v*w*m_midedgenormals[1] + w*u*m_midedgenormals[2]), Vector3(zero), 0);

	float mx = max(max(u,v),w);

	m_image.addPixel(coverage * normal, x, y, m_image.normalChannel());

	/*if (m_image.occlusionChannel() != -1)
	{
		float occlusion = sampleOcclusion(position, normal, x, y);
		nvCheck(occlusion >= 0.0f && occlusion <= 1.0f);

		m_image.addPixel(coverage * occlusion, x, y, m_image.occlusionChannel());
	}*/

	m_image.addPixel(coverage, x, y, m_image.coverageChannel());
	
	m_imageMask.setBitAt(x, y);
}


static inline float triangleArea(Vector2::Arg a, Vector2::Arg b, Vector2::Arg c)
{
	Vector2 v0 = a - c;
	Vector2 v1 = b - c;

	return (v0.x() * v1.y() - v0.y() * v1.x());
}

void GeometrySampler::sampleQuad(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	nvDebugCheck(isFinite(coverage));

	const float u = bar.x(), U = 1-u;
	const float v = bar.y(), V = 1-v;

	// bilinear position interpolation
	Vector3 position = (U * m_positions[0] + u * m_positions[1]) * V + 
	                   (U * m_positions[3] + u * m_positions[2]) * v;

	m_image.addPixel(coverage * position, x, y, m_image.positionChannel());

	
	// biquadratic normal interpolation
	Vector3 normal = normalizeSafe(
		(U*U * m_normals[0] + 2*U*u*m_midedgenormals[0] + u*u * m_normals[1]) * V*V + 
		(U*U*m_midedgenormals[3] + 2*U*u*m_midedgenormals[4] + u*u*m_midedgenormals[1]) * 2*V*v +
		(U*U * m_normals[3] + 2*U*u*m_midedgenormals[2] + u*u * m_normals[2]) * v*v, Vector3(zero), 0);

	m_image.addPixel(coverage * normal, x, y, m_image.normalChannel());

	/*
	if (m_image.occlusionChannel() != -1)
	{
		// piecewise linear position interpolation
		#if 0
		Vector3 tripos;
		if (u < v)
		{
			float barx = triangleArea(Vector2(1,1), Vector2(0,1), Vector2(u,v));
			float bary = triangleArea(Vector2(0,0), Vector2(1,1), Vector2(u,v));
			float barz = triangleArea(Vector2(0,1), Vector2(0,0), Vector2(u,v));
			nvCheck(equal(1, barx+bary+barz));

			tripos = barx * m_positions[0] + bary * m_positions[1] + barz * m_positions[2];
		}
		else 
		{
			float barx = triangleArea(Vector2(1,0), Vector2(1,1), Vector2(u,v));
			float bary = triangleArea(Vector2(1,1), Vector2(0,0), Vector2(u,v));
			float barz = triangleArea(Vector2(0,0), Vector2(1,0), Vector2(u,v));
			nvCheck(equal(1, barx+bary+barz));

			tripos = barx * m_positions[0]  + bary * m_positions[3] + barz * m_positions[2];
		}
		#endif

		float occlusion = sampleOcclusion(position, normal, x, y);
		nvCheck(occlusion >= 0.0f && occlusion <= 1.0f);
	
		m_image.addPixel(coverage * occlusion, x, y, m_image.occlusionChannel());
	}
	*/

	m_image.addPixel(coverage, x, y, m_image.coverageChannel());
	
	m_imageMask.setBitAt(x, y);
}



DisplacementPatchSampler::DisplacementPatchSampler(const GeometryImage & detailedGeometryMap, const BitMap & detailedGeometryMask, GeometryImage & geometryMap, BitMap & geometryMask) :
	m_detailedGeometryMap(detailedGeometryMap), 
	m_detailedGeometryMask(detailedGeometryMask),
	m_geometryMap(geometryMap), 
	m_geometryMask(geometryMask)
{
	nvCheck(geometryMap.width() == geometryMask.width());
	nvCheck(geometryMap.height() == geometryMask.height());
}

/*static*/ void DisplacementPatchSampler::sampleCallback(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	((DisplacementPatchSampler *)param)->sample(x, y, bar, dx, dy, coverage);
}


// setDetailedGeometryMap(const GeometryImage & detailedGeometryMap, const BitMap & detailedGeometryMask);
// setOutputGeometryMap(GeometryImage & detailedGeometryMap, BitMap & detailedGeometryMask);

void DisplacementPatchSampler::setTangentSpace(bool enabled)
{
	m_tangentSpace = enabled;
}

void DisplacementPatchSampler::setVectorDisplacement(bool enabled)
{
	m_vectorDisplacement = enabled;
}

void DisplacementPatchSampler::setCurrentSurface(const BaseSurface * surface)
{
	m_surface = surface;
}


void DisplacementPatchSampler::sample(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	const float u = bar.x();
	const float v = bar.y();

	Vector3 origin;
	Basis patchFrame;
	Basis chartFrame;

	m_surface->evaluate(u, v, &origin, &patchFrame, &chartFrame);

	float normalDeviation = dot(chartFrame.normal, patchFrame.normal);
	patchFrame.normal = chartFrame.normal;

#if 0
	Vector3 tangentSpaceNormal = chartFrame.tangent;
//	tangentSpaceNormal = chartFrame.bitangent;
//	tangentSpaceNormal = chartFrame.normal;

	m_geometryMap.addPixel(coverage * tangentSpaceNormal, x, y, m_geometryMap.normalChannel());

	// Update coverage.
	m_geometryMap.addPixel(coverage, x, y, m_geometryMap.coverageChannel());
	m_geometryMask.setBitAt(x, y);
	
	return;
#endif

	// assuming coverage has been normalized before!!
	if (m_detailedGeometryMask.bitAt(x, y))
	{
		// Compute tangent space normal.
		if (m_detailedGeometryMap.normalChannel() != -1)
		{
			if (m_tangentSpace)
			{
				const Vector3 normal = 2 * m_detailedGeometryMap.normal(x, y) - 1;

				Vector3 tangentSpaceNormal = chartFrame.transformI(normal);
				tangentSpaceNormal = normalizeSafe(tangentSpaceNormal, Vector3(zero), 0);
			//	tangentSpaceNormal = Vector3(normalDeviation, normalDeviation, normalDeviation);

				m_geometryMap.addPixel(coverage * tangentSpaceNormal, x, y, m_geometryMap.normalChannel());
			}
		}

		if (m_detailedGeometryMap.positionChannel() != -1)
		{
			// Compute displacement.
			const Vector3 dest = m_detailedGeometryMap.position(x, y);
			const Vector3 displacement = dest - origin;

			if (m_vectorDisplacement)
			{
				if (m_tangentSpace)
				{
					Vector3 tangentSpaceDisplacement = chartFrame.transformI(displacement);
					m_geometryMap.addPixel(coverage * tangentSpaceDisplacement, x, y, m_geometryMap.displacementChannel());
				}
				else
				{
					m_geometryMap.addPixel(coverage * displacement, x, y, m_geometryMap.displacementChannel());
				}
			}
			else
			{
				float l = dot(chartFrame.normal, displacement);
				m_geometryMap.addPixel(coverage * l, x, y, m_geometryMap.displacementChannel());
			}
		}

		// Update coverage.
		m_geometryMap.addPixel(coverage, x, y, m_geometryMap.coverageChannel());
		m_geometryMask.setBitAt(x, y);
	}
}




// SCALAR DISPLACEMENTS		
/*
void sampleAA(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	nvDebugCheck(m_positions != NULL);
	nvDebugCheck(m_normals != NULL);

	nvDebugCheck(isFinite(coverage));

	const float u = bar.y();
	const float v = bar.x();

	Vector3 origin = ((1-u) * m_positions[0] + (u) * m_positions[1]) * (1-v) + 
					 ((1-u) * m_positions[3] + (u) * m_positions[2]) * (v);

	Vector3 normal = normalizeSafe(
		((1-u) * m_normals[0] + (u) * m_normals[1]) * (1-v) + 
		((1-u) * m_normals[3] + (u) * m_normals[2]) * (v), Vector3(zero), 0.0f);

	Vector3 dest = Vector3(
		m_geometryMap.pixel(x, y, 0),
		m_geometryMap.pixel(x, y, 1),
		m_geometryMap.pixel(x, y, 2));

	dest /= m_geometryMap.pixel(x, y, 6); // normalize just to be sure
	
	Vector3 displacement = dest - origin;
	
	// Make sure that displacement is roughly in the right direction.
//	if(!equal(dot(normal, normalize(displacement)), 1, NV_NORMAL_EPSILON))
//	{
	//	float diff = dot(normal, normalize(displacement));
	//	printf("%f\n", diff);
	//	nvDebugBreak();
//	}
	
	// @@ Project displacement onto normal?

	if (length(displacement) < m_outlierThreshold)
	{
		m_displacementMap.addPixel(coverage * length(displacement), x, y, 0);
		m_displacementMap.addPixel(coverage, x, y, 1);
		m_displacementMask.setBitAt(x, y);
	}
}
*/


// Vector Displacement from flat quads
/*		
void sampleAA(int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage)
{
	nvDebugCheck(m_positions != NULL);
	nvDebugCheck(m_normals != NULL);

	const float u = bar.y();
	const float v = bar.x();

	Vector3 origin = ((1-u) * m_positions[0] + (u) * m_positions[1]) * (1-v) + 
					 ((1-u) * m_positions[3] + (u) * m_positions[2]) * (v);

	Vector3 destiny = Vector3(
		m_geometryMap.pixel(x, y, 0),
		m_geometryMap.pixel(x, y, 1),
		m_geometryMap.pixel(x, y, 2));  // assuming subpixel-coverage normalized to 1
	
	Vector3 displacement = destiny - origin;
	
	if (length(displacement) < m_outlierThreshold)
	{
		m_displacementMap.setPixel(displacement.x(), x, y, 0);
		m_displacementMap.setPixel(displacement.y(), x, y, 1);
		m_displacementMap.setPixel(displacement.z(), x, y, 2);

		m_displacementMask.setBitAt(x, y);
	}
}
*/
