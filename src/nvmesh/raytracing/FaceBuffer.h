// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#ifndef NV_MESH_FACEBUFFER_H
#define NV_MESH_FACEBUFFER_H

#include <nvcore/Ptr.h>
#include <nvcore/Containers.h> // Array

#include <nvmath/Plane.h>

#include <nvmesh/raytracing/Raytracing.h>


namespace nv
{
	class Stream;
	class TriMesh;
	class QuadTriMesh;
	

	class FaceBuffer
	{
		NV_FORBID_COPY(FaceBuffer);
	public:
		
		FaceBuffer();
		
		uint faceCount() const;

		void resetTestCount();
		uint testCount() const;

		// Interface.
		virtual void build(const Array<uint> & indices, const Array<Vector3> & vertices) = 0;
		
		virtual void testRay(uint first, uint faceCount, const Ray & ray, Hit * hit) const = 0;
		virtual void testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const = 0;

		virtual bool testRay(uint first, uint faceCount, const Ray & ray) const;
		virtual bool testRay(const uint * faceIndices, uint faceCount, const Ray & ray) const;

		// @@ Add double sided face tests?
		
	protected:
		
		uint m_faceCount;
		
		mutable uint m_testCount;

	};

	
	// Christer Ericson ray triangle test.
	class EricsonFaceBuffer : public FaceBuffer
	{
	public:
		void build(const Array<uint> & indices, const Array<Vector3> & vertices);
		
		virtual void testRay(uint first, uint num, const Ray & ray, Hit * hit) const;
		virtual void testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const;

	private:
		struct Face {};

		void testFace(uint f, const Ray & ray, Hit * hit) const;

		Array<Face> m_faceArray;
	};
	
	// Charles Bloom ray triangle test.
	class BloomFaceBuffer : public FaceBuffer
	{
	public:
		void build(const Array<uint> & indices, const Array<Vector3> & vertices);

		//bool testRay(uint first, uint num, const Ray & ray) const;
		virtual void testRay(uint first, uint num, const Ray & ray, Hit * hit) const;
		virtual void testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const;

		friend Stream & operator<< (Stream & s, BloomFaceBuffer & faceBuffer);

	private:

		struct Face
		{
			Plane plane;
			Plane bary1;
			Plane bary2;
		};

		void testFace(uint f, const Ray & ray, Hit * hit) const;

		friend Stream & operator<< (Stream & s, Face & face);

		Array<Face> m_faceArray;
	};
	
	// Ingo Wald ray triangle test.
	class WaldFaceBuffer : public FaceBuffer
	{
	public:
		void build(const Array<uint> & indices, const Array<Vector3> & vertices);

		//bool testRay(uint first, uint num, const Ray & ray) const;
		virtual void testRay(uint first, uint num, const Ray & ray, Hit * hit) const;
		virtual void testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const;

		friend Stream & operator<< (Stream & s, WaldFaceBuffer & faceBuffer);

	private:

		/*// Ingo Wald face.
		struct Face
		{
			uint k;		// projection axis
			float n_u;	// normal.u / normal.k
			float n_v;	// normal.v / normal.k
			float n_d;	// distance

			float b_nu;
			float b_nv;
			float b_d;
			mutable uint counter;
			//uint pad0;
			
			float c_nu;
			float c_nv;
			float c_d;
			uint pad1;
		};*/

		// Modified Ingo Wald face.
		struct Face
		{
			mutable uint counter;
			Vector3 n;

			float d;
			uint k;		// projection axis

			float b_nu;
			float b_nv;
			float b_d;

			float c_nu;
			float c_nv;
			float c_d;
		};

		void testFace(uint f, const Ray & ray, Hit * hit) const;

		friend Stream & operator<< (Stream & s, Face & face);

		Array<Face> m_faceArray;
	};
	
	// Tomas Moeller ray triangle test.
	class MollerFaceBuffer : public FaceBuffer
	{
	public:
		void build(const Array<uint> & indices, const Array<Vector3> & vertices);

		//bool testRay(uint first, uint num, const Ray & ray) const;
		virtual void testRay(uint first, uint num, const Ray & ray, Hit * hit) const;
		virtual void testRay(const uint * faceIndices, uint faceCount, const Ray & ray, Hit * hit) const;

		friend Stream & operator<< (Stream & s, MollerFaceBuffer & faceBuffer);

	private:

		struct Face
		{
			mutable uint counter;
			Vector3 v0;
			
			uint pad0;
			Vector3 v1;

			uint pad1;
			Vector3 v2;
		};

		void testFace(uint f, const Ray & ray, Hit * hit) const;

		friend Stream & operator<< (Stream & s, Face & face);

		Array<Face> m_faceArray;
	};
	
} // nv namespace


#endif // NV_MESH_FACEBUFFER_H

