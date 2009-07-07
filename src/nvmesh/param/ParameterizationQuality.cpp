// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <float.h>

#include <nvcore/Debug.h>

#include <nvmath/Vector.h>

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>

#include "ParameterizationQuality.h"

using namespace nv;

/*
float triangleConformalEnergy(Vector3 q[3], Vector2 p[3])
{
	const Vector3 v1 = q[0];
	const Vector3 v2 = q[1];
	const Vector3 v3 = q[2];

	const Vector2 w1 = p[0];
	const Vector2 w2 = p[1];
	const Vector2 w3 = p[2];

	float x1 = v2.x() - v1.x();
	float x2 = v3.x() - v1.x();
	float y1 = v2.y() - v1.y();
	float y2 = v3.y() - v1.y();
	float z1 = v2.z() - v1.z();
	float z2 = v3.z() - v1.z();

	float s1 = w2.x() - w1.x();
	float s2 = w3.x() - w1.x();
	float t1 = w2.y() - w1.y();
	float t2 = w3.y() - w1.y();

	float r = 1.0f / (s1 * t2 - s2 * t1);
	Vector3 sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r, (t2 * z1 - t1 * z2) * r);
	Vector3 tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r, (s1 * z2 - s2 * z1) * r);

	Vector3 N = cross(v3-v1, v2-v1);

	// Rotate 90 around N.
}
*/

static float triangleConformalEnergy(Vector3 q[3], Vector2 p[3])
{
	// Using Denis formulas:
	Vector3 c0 = q[1] - q[2];
	Vector3 c1 = q[2] - q[0];
	Vector3 c2 = q[0] - q[1];

	Vector3 N = cross(-c0, c1);
	float T = length(N);	// 2T
	N = normalize(N, 0);

	float cot_alpha0 = dot(-c1, c2) / length(cross(-c1, c2));
	float cot_alpha1 = dot(-c2, c0) / length(cross(-c2, c0));
	float cot_alpha2 = dot(-c0, c1) / length(cross(-c0, c1));

	Vector3 t0 = -cot_alpha1 * c1 + cot_alpha2 * c2;
	Vector3 t1 = -cot_alpha2 * c2 + cot_alpha0 * c0;
	Vector3 t2 = -cot_alpha0 * c0 + cot_alpha1 * c1;

	nvCheck(equal(length(t0), length(c0)));
	nvCheck(equal(length(t1), length(c1)));
	nvCheck(equal(length(t2), length(c2)));
	nvCheck(equal(dot(t0, c0), 0));
	nvCheck(equal(dot(t1, c1), 0));
	nvCheck(equal(dot(t2, c2), 0));

	// Gradients
	Vector3 grad_u = 1.0f / T * (p[0].x() * t0 + p[1].x() * t1 + p[2].x() * t2);
	Vector3 grad_v = 1.0f / T * (p[0].y() * t0 + p[1].y() * t1 + p[2].y() * t2);

	// Rotated gradients
	Vector3 Jgrad_u = 1.0f / T * (p[0].x() * c0 + p[1].x() * c1 + p[2].x() * c2);
	Vector3 Jgrad_v = 1.0f / T * (p[0].y() * c0 + p[1].y() * c1 + p[2].y() * c2);

	// Using Lengyel's formulas:
	{ 
		const Vector3 v1 = q[0];
		const Vector3 v2 = q[1];
		const Vector3 v3 = q[2];

		const Vector2 w1 = p[0];
		const Vector2 w2 = p[1];
		const Vector2 w3 = p[2];

		float x1 = v2.x() - v1.x();
		float x2 = v3.x() - v1.x();
		float y1 = v2.y() - v1.y();
		float y2 = v3.y() - v1.y();
		float z1 = v2.z() - v1.z();
		float z2 = v3.z() - v1.z();

		float s1 = w2.x() - w1.x();
		float s2 = w3.x() - w1.x();
		float t1 = w2.y() - w1.y();
		float t2 = w3.y() - w1.y();

		float r = 1.0f / (s1 * t2 - s2 * t1);
		Vector3 sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r, (t2 * z1 - t1 * z2) * r);
		Vector3 tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r, (s1 * z2 - s2 * z1) * r);

		Vector3 Jsdir = cross(N, sdir);
		Vector3 Jtdir = cross(N, tdir);

		float x = 3;
	}

	// check: sdir == grad_u
	// check: tdir == grad_v

	return length(grad_u - Jgrad_v);
}


ParameterizationQuality::ParameterizationQuality(const HalfEdge::Mesh * mesh)
{
	nvCheck(mesh != NULL);

	m_totalTriangleCount = 0;
	m_flippedTriangleCount = 0;
	
	m_parametricArea = 0.0f;
	m_geometricArea = 0.0f;
	
	m_rmsStretchMetric = 0.0f;
	m_maxStretchMetric = 0.0f;

	m_rmsConformalMetric = 0.0f;
	m_rmsAuthalicMetric = 0.0f;

	const uint faceCount = mesh->faceCount();
	for (uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = mesh->faceAt(f);
		nvCheck(face != NULL);
		
		const HalfEdge::Vertex * vertex0 = NULL;

		Vector3 p[3];
		Vector2 t[3];
		
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			const HalfEdge::Edge * edge = it.current();
			nvCheck(edge != NULL);
			
			if (vertex0 == NULL)
			{
				vertex0 = edge->vertex();

				p[0] = vertex0->pos();
				t[0] = vertex0->tex();
			}
			else if (edge->to() != vertex0)
			{
				p[1] = edge->from()->pos();
				p[2] = edge->to()->pos();
				t[1] = edge->from()->tex();
				t[2] = edge->to()->tex();
				
				processTriangle(p, t);
			}
		}
	}
	
	float normFactor = sqrtf(m_parametricArea / m_geometricArea);
	
	m_rmsStretchMetric = sqrtf(m_rmsStretchMetric / m_geometricArea) * normFactor;
	m_maxStretchMetric = m_maxStretchMetric * normFactor;

	m_rmsConformalMetric = sqrtf(m_rmsConformalMetric / m_geometricArea);
	m_rmsAuthalicMetric = sqrtf(m_rmsAuthalicMetric / m_geometricArea);
}

bool ParameterizationQuality::isValid() const
{
	return m_flippedTriangleCount == 0 || m_flippedTriangleCount == m_totalTriangleCount;
}

float ParameterizationQuality::rmsStretchMetric() const
{
	return m_rmsStretchMetric;
}

float ParameterizationQuality::maxStretchMetric() const
{
	return m_maxStretchMetric;
}

float ParameterizationQuality::rmsConformalMetric() const
{
	return m_rmsConformalMetric;
}

float ParameterizationQuality::maxAuthalicMetric() const
{
	return m_rmsAuthalicMetric;
}



void ParameterizationQuality::processTriangle(Vector3 q[3], Vector2 p[3])
{
	// Evaluate texture stretch metric. See:
	// - "Texture Mapping Progressive Meshes", Sander, Snyder, Gortler & Hoppe
	// - "Mesh Parameterization: Theory and Practice", Siggraph'07 Course Notes, Hormann, Levy & Sheffer.

	float t1 = p[0].x();
	float s1 = p[0].y();
	float t2 = p[1].x();
	float s2 = p[1].y();
	float t3 = p[2].x();
	float s3 = p[2].y();
	
	float geometricArea = length(cross(q[1] - q[0], q[2] - q[0])) / 2;
	float parametricArea = ((s2 - s1)*(t3 - t1) - (s3 - s1)*(t2 - t1)) / 2;
	Vector3 Ss = (q[0] * (t2- t3) + q[1] * (t3 - t1) + q[2] * (t1 - t2)) / (2 * parametricArea);
	Vector3 St = (q[0] * (s3- s2) + q[1] * (s1 - s3) + q[2] * (s2 - s1)) / (2 * parametricArea);
	
	float a = dot(Ss, Ss); // E
	float b = dot(Ss, St); // F
	float c = dot(St, St); // G
	
	// Compute eigen-values of the first fundamental form:
	float sigma1 = sqrtf(0.5f * (a + c - sqrtf(square(a - c) + 4 * square(b)))); // gamma uppercase, min eigenvalue.
	float sigma2 = sqrtf(0.5f * (a + c + sqrtf(square(a - c) + 4 * square(b)))); // gamma lowercase, max eigenvalue.
	nvCheck(sigma2 >= sigma1);

	// isometric: sigma1 = sigma2 = 1
	// conformal: sigma1 / sigma2 = 1
	// authalic: sigma1 * sigma2 = 1

	float rmsStretch = sqrtf((a + c) * 0.5f);
	float rmsStretch2 = sqrtf((square(sigma1) + square(sigma2)) * 0.5f);
	nvCheck(equal(rmsStretch, rmsStretch2, 0.01f));

	m_totalTriangleCount++;

	if (parametricArea < 0.0f)
	{
		// Count flipped triangles.
		m_flippedTriangleCount++;

		parametricArea = fabsf(parametricArea);
	}

	m_rmsStretchMetric += square(rmsStretch) * geometricArea;
	m_maxStretchMetric = max(m_maxStretchMetric, sigma2);

	m_rmsConformalMetric += (sigma2 / sigma1) * geometricArea;
	m_rmsAuthalicMetric += (sigma1 * sigma2) * geometricArea;

	// Accumulate total areas.
	m_geometricArea += geometricArea;
	m_parametricArea += fabsf(parametricArea);
	
	triangleConformalEnergy(q, p);
}
