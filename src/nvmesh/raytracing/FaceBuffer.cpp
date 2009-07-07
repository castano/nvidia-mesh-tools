// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#include <nvcore/Prefetch.h>

#include <nvmath/TypeSerialization.h>

#include "FaceBuffer.h"

#define COMPUTE_STATS   0
#define USE_MAILBOX     1
#define BACKFACE_CULL   1

using namespace nv;


FaceBuffer::FaceBuffer() : m_faceCount(0), m_testCount(0)
{
}

uint FaceBuffer::faceCount() const
{
	return m_faceCount;
}

void FaceBuffer::resetTestCount()
{
	m_testCount = 0;
}

uint FaceBuffer::testCount() const
{
	return m_testCount;
}

bool FaceBuffer::testRay(uint first, uint faceCount, const Ray & ray) const
{
	return false;
}

bool FaceBuffer::testRay(const uint * faceIndices, uint faceCount, const Ray & ray) const
{
	return false;
}



// Christer Ericson ray triangle test.
void EricsonFaceBuffer::build(const Array<uint> & indices, const Array<Vector3> & vertices)
{
}

inline void EricsonFaceBuffer::testFace(uint f, const Ray & ray, Hit * hit) const
{
}

/*bool EricsonFaceBuffer::testRay(uint first, uint num, const Ray & ray) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	Hit hit;
	hit.face = -1;
	
	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, &hit);
	}

	return hit.face != -1;
}*/

void EricsonFaceBuffer::testRay(uint first, uint num, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, hit);
	}
}

void EricsonFaceBuffer::testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint f = 0; f < faceCount - 1; f++)
	{
		// Prefetch next face.
		nvPrefetch(m_faceArray.buffer() + faceIndices[f+1]);

		// Test face and update nearest hit.
		testFace(faceIndices[f], ray, hit);
	}

	// Test face and update nearest hit.
	testFace(faceIndices[faceCount - 1], ray, hit);
}


// Charles Bloom ray triangle test.
void BloomFaceBuffer::build(const Array<uint> & indices, const Array<Vector3> & vertices)
{
	m_faceCount = indices.count() / 3;
	m_faceArray.resize(m_faceCount);

	for (uint f = 0; f < m_faceCount; f++)
	{
		const Vector3 v0 = vertices[indices[f * 3 + 0]];
		const Vector3 v1 = vertices[indices[f * 3 + 1]];
		const Vector3 v2 = vertices[indices[f * 3 + 2]];

		Face & face = m_faceArray[f];

		const Vector3 edge1 = v1 - v0;
		const Vector3 edge2 = v2 - v1;

		const Vector3 n = cross(edge1, edge2);
		face.plane = normalize(Plane(n, dot(n, v1)));

		face.bary1 = Plane(cross(face.plane.vector(), edge1), v0);
		face.bary2 = Plane(cross(face.plane.vector(), edge2), v1);
		
		face.bary1 *= 1.0f / distance(face.bary1, v2);
		face.bary2 *= 1.0f / distance(face.bary1, v0);
		
		nvDebugCheck( equal(0.0f, distance(face.bary1, v0)) );
		nvDebugCheck( equal(0.0f, distance(face.bary1, v1)) );
		nvDebugCheck( equal(1.0f, distance(face.bary1, v2)) );
		nvDebugCheck( equal(1.0f, distance(face.bary2, v0)) );
		nvDebugCheck( equal(0.0f, distance(face.bary2, v1)) );
		nvDebugCheck( equal(0.0f, distance(face.bary2, v2)) );
	}
}

inline void BloomFaceBuffer::testFace(uint f, const Ray & ray, Hit * hit) const
{
	nvDebugCheck(hit != NULL);

	const Face & face = m_faceArray[f];

	const Vector3 fm = ray.origin;
	const Vector3 to = ray.origin + ray.dir * ray.maxt;
	
	const float dTo = distance(face.plane, to);
	if (dTo <= 0)
	{
		// segment doesn't reach back side of triangle
		return;
	}

	const float dFm = distance(face.plane, fm);
	if (dFm > 0)
	{
		// segment starts already on the back side
		return;
	}
	nvDebugCheck( dTo > 0 && dFm <= 0 );

	const float denom = dTo - dFm;
	nvDebugCheck( denom > 0 );

	const float uTimesDenom = dTo * dot(face.bary1.vector(), fm) - dFm * dot(face.bary1.vector(), to) - face.bary1.offset() * denom;
	if (uTimesDenom < 0 || uTimesDenom > denom)
	{
		// off triangle
		return;
	}

	const float vTimesDenom = dTo * dot(face.bary2.vector(), fm) - dFm * dot(face.bary2.vector(), to) - face.bary2.offset() * denom;
	if (vTimesDenom < 0 || (uTimesDenom + vTimesDenom) > denom)
	{
		// off triangle
		return;
	}

	// it's on the triangle, compute hit info:
	const float inv = 1.0f / denom;
	const float t = hit->t * - dFm * inv;

	hit->t = t;
	hit->face = f;
	hit->u = uTimesDenom * inv;
	hit->v = vTimesDenom * inv;
}

/*bool BloomFaceBuffer::testRay(uint first, uint num, const Ray & ray) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	Hit hit;
	hit.face = -1;
	
	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, &hit);
	}

	return hit.face != -1;
}*/

void BloomFaceBuffer::testRay(uint first, uint num, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, hit);
	}
}

void BloomFaceBuffer::testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint f = 0; f < faceCount - 1; f++)
	{
		// Prefetch next face.
		nvPrefetch(m_faceArray.buffer() + faceIndices[f+1]);

		// Test face and update nearest hit.
		testFace(faceIndices[f], ray, hit);
	}

	// Test face and update nearest hit.
	testFace(faceIndices[faceCount - 1], ray, hit);
}


// Ingo Wald ray triangle test.
// http://www.mpi-sb.mpg.de/~wald/PhD/
// This test has been modified to perform the face plane test in world space 
// instead of projection space. It also does mailboxing and backface culling.

void WaldFaceBuffer::build(const Array<uint> & indices, const Array<Vector3> & vertices)
{
	m_faceCount = indices.count() / 3;
	m_faceArray.resize(m_faceCount);

	for (uint f = 0; f < m_faceCount; f++)
	{
		const Vector3 v0 = vertices[indices[f * 3 + 0]];
		const Vector3 v1 = vertices[indices[f * 3 + 1]];
		const Vector3 v2 = vertices[indices[f * 3 + 2]];

		// Compute face acceleration info for Ingo Wald ray-face test.
		Vector3 c, b, n;
		c = v1 - v0;
		b = v2 - v0;
		n = normalize(cross(c, b), 0.0f);

		// Get largest axis.
		int k = 2;
		if( fabsf(n.x()) > fabsf(n.y()) ) {
			if( fabsf(n.x()) > fabsf(n.z()) ) {
				k = 0;
			}
		}
		else {
			if( fabsf(n.y()) > fabsf(n.z()) ) {
				k = 1;
			}
		}

		int u = (1 << k) & 3;
		int v = (1 << u) & 3;

		Face & face = m_faceArray[f];
		face.counter = 0;

	//	float i_n_k = 1.0f / n[k];
	//	face.n_u = n[u] * i_n_k;
	//	face.n_v = n[v] * i_n_k;
	//	face.n_d = dot(v0, n) * i_n_k;

		face.n = n;
		face.d = dot(v0, n);
		face.k = k;
		
		float inv = 1.0f / (b.component(u) * c.component(v) - b.component(v) * c.component(u));
		face.b_nu = b.component(u) * inv;
		face.b_nv = -b.component(v) * inv;
		face.b_d = v0.component(u) * face.b_nv + v0.component(v) * face.b_nu;

		face.c_nu = c.component(v) * inv;
		face.c_nv = -c.component(u) * inv;
		face.c_d = v0.component(u) * face.c_nu + v0.component(v) * face.c_nv;
	}
}

inline void WaldFaceBuffer::testFace(uint f, const Ray & ray, Hit * hit) const
{
	nvDebugCheck(hit != NULL);

	const Face & face = m_faceArray[f];

#if USE_MAILBOX
	if (ray.id == face.counter) {
		return;
	}

	face.counter = ray.id;
#endif

	const float nd = dot(ray.dir, face.n);

#if BACKFACE_CULL
	// Check for backfacing or coplanar faces.
	if (nd > -NV_EPSILON) {
		return;
	}
#endif

	const float ind = 1.0f / nd;
	const float t = (face.d - dot(ray.origin, face.n)) * ind;

	// Check for valid distance.
	if (hit->t < t || t < NV_EPSILON) {
		return;
	}

	const uint ku = (1 << face.k) & 3;
	const uint kv = (1 << ku) & 3;

	// Compute hitpoint position on uv plane.
	const float hu = (ray.origin.component(ku) + t * ray.dir.component(ku));
	const float hv = (ray.origin.component(kv) + t * ray.dir.component(kv));

	// Check first barycentric coordinate.
	const float lambda = (hv * face.b_nu + hu * face.b_nv - face.b_d);
	if (lambda < 0.0f) {
		return;
	}

	// Check second barycentric coordinate.
	const float mue = (hu * face.c_nu + hv * face.c_nv - face.c_d);
	if (mue < 0.0f) {
		return;
	}

	if (lambda + mue > 1.0f) {
		return;
	}

	hit->t = t;
	hit->face = f;
	hit->u = lambda;
	hit->v = mue;
}

/*bool WaldFaceBuffer::testRay(uint first, uint num, const Ray & ray) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	Hit hit;
	hit.face = -1;
	
	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, &hit);
	}

	return hit.face != -1;
}*/

void WaldFaceBuffer::testRay(uint first, uint num, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, hit);
	}
}

void WaldFaceBuffer::testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint f = 0; f < faceCount - 1; f++)
	{
		// Prefetch next face.
		nvPrefetch(m_faceArray.buffer() + faceIndices[f+1]);

		// Test face and update nearest hit.
		testFace(faceIndices[f], ray, hit);
	}

	// Test face and update nearest hit.
	testFace(faceIndices[faceCount - 1], ray, hit);
}


// Tomas Moeller ray triangle test.
void MollerFaceBuffer::build(const Array<uint> & indices, const Array<Vector3> & vertices)
{
	m_faceCount = indices.count() / 3;
	m_faceArray.resize(m_faceCount);

	for (uint f = 0; f < m_faceCount; f++)
	{
		const Vector3 v0 = vertices[indices[f * 3 + 0]];
		const Vector3 v1 = vertices[indices[f * 3 + 1]];
		const Vector3 v2 = vertices[indices[f * 3 + 2]];

		Face & face = m_faceArray[f];
		face.counter = 0;
		face.v0 = v0;
		face.v1 = v1;
		face.v2 = v2;
	}
}

inline void MollerFaceBuffer::testFace(uint f, const Ray & ray, Hit * hit) const
{
	nvDebugCheck(hit != NULL);

	const Face & face = m_faceArray[f];

#if USE_MAILBOX
	if (ray.id == face.counter) {
		return;
	}

	face.counter = ray.id;
#endif
	
	const Vector3 dir = ray.dir; // * ray.maxt;
	
	/* find vectors for two edges sharing vert0 */
	const Vector3 edge1 = face.v1 - face.v0;
	const Vector3 edge2 = face.v2 - face.v0;
	
	/* begin calculating determinant - also used to calculate U parameter */
	const Vector3 pvec = cross(dir, edge2);
	
	/* if determinant is near zero, ray lies in plane of triangle */
	const float det = dot(edge1, pvec);

	Vector3 qvec;
	float u, v;

  	if (det > NV_EPSILON)
  	{
		/* calculate distance from vert0 to ray origin */
		const Vector3 tvec = ray.origin - face.v0;
		
		/* calculate U parameter and test bounds */
		u = dot(tvec, pvec);
		if (u < 0.0f || u > det)
			return;
		
		/* prepare to test V parameter */
		qvec = cross(tvec, edge1);
		
		/* calculate V parameter and test bounds */
		v = dot(dir, qvec);
		if (v < 0.0f || u + v > det)
			return;
	}
#if !BACKFACE_CULL
	else if (det < -NV_EPSILON)
	{
		/* calculate distance from vert0 to ray origin */
		const Vector3 tvec = ray.origin - face.v0;
		
		/* calculate U parameter and test bounds */
		u = dot(tvec, pvec);
		if (u > 0.0f || u < det)
			return;
		
		/* prepare to test V parameter */
		qvec = cross(tvec, edge1);
		
		/* calculate V parameter and test bounds */
		v = dot(dir, qvec);
		if (v > 0.0f || u + v < det)
			return;
	}
#endif
	else
	{
		return;  /* ray is parallell to the plane of the triangle */
	}

	const float inv_det = 1.f / det;
	
	const float t = dot(edge2, qvec) * inv_det;
	
	// Check for valid distance.
	if (hit->t < t || t < NV_EPSILON) {
		return;
	}
	
	// calculate t, ray intersects triangle
	hit->t = t;
	hit->face = f;
	hit->u = u * inv_det;
	hit->v = v * inv_det;
}

/*bool MollerFaceBuffer::testRay(uint first, uint num, const Ray & ray) const
{
	Hit hit;
	hit.face = -1;
	
	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, &hit);
	}

	return hit.face != -1;
}*/

void MollerFaceBuffer::testRay(uint first, uint num, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint i = 0; i < num; i++)
	{
		testFace(first + i, ray, hit);
	}
}

void MollerFaceBuffer::testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const
{
#if COMPUTE_STATS
	m_testCount += num;
#endif

	for (uint f = 0; f < faceCount - 1; f++)
	{
		// Prefetch next face.
		nvPrefetch(m_faceArray.buffer() + faceIndices[f+1]);

		// Test face and update nearest hit.
		testFace(faceIndices[f], ray, hit);
	}

	// Test face and update nearest hit.
	testFace(faceIndices[faceCount - 1], ray, hit);
}


namespace nv
{
	Stream & operator<< (Stream & s, BloomFaceBuffer::Face & face)
	{
		return s << face.plane << face.bary1 << face.bary2;
	}

	Stream & operator<< (Stream & s, BloomFaceBuffer & faceBuffer)
	{
		s << faceBuffer.m_faceCount;
		s << faceBuffer.m_testCount;
		s << faceBuffer.m_faceArray;
		return s;
	}


	Stream & operator<< (Stream & s, WaldFaceBuffer::Face & face)
	{
		s << face.counter << face.n;
		s << face.d << face.k;
		s << face.b_nu << face.b_nv << face.b_d;
		s << face.c_nu << face.c_nv << face.c_d;
		return s;
	}

	Stream & operator<< (Stream & s, WaldFaceBuffer & faceBuffer)
	{
		s << faceBuffer.m_faceCount;
		s << faceBuffer.m_testCount;
		s << faceBuffer.m_faceArray;
		return s;
	}


	Stream & operator<< (Stream & s, MollerFaceBuffer::Face & face)
	{
		return s << face.counter << face.v0 << face.pad0 << face.v1 << face.pad1 << face.v2;
	}

	Stream & operator<< (Stream & s, MollerFaceBuffer & faceBuffer)
	{
		s << faceBuffer.m_faceCount;
		s << faceBuffer.m_testCount;
		s << faceBuffer.m_faceArray;
		return s;
	}

	
}

