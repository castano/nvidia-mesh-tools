// Copyright NVIDIA Corporation 2008 -- Edgar Velazquez-Armendariz <edgarv@nvidia.com>

#include "Array2D.h"

#include <nvcore/Algorithms.h>

using namespace nv;

const int HilbertIterator::R2D[4][2][2] = { {{0,1},{1,0}},{{1,0},{0,1}},
							{{1,0},{0,1}},{{0,-1},{-1,0}} };

const int HilbertIterator::xbit2D[4] = {0, 0, 1, 1};
const int HilbertIterator::ybit2D[4] = {0, 1, 1, 0};


// ------ Morton Iterator ------

MortonIterator::MortonIterator(uint w, uint h) : Iterator2D(w, h),
	m_currentIndex(0, 0)
{
	nvDebugCheck(w <= 0x10000 && h <= 0x10000);
}

Index2D MortonIterator::current() const
{
	nvDebugCheck(m_counter >= 0 && m_counter <= m_width*m_height);
	return m_currentIndex;
}

bool MortonIterator::moveNext()
{
	nvDebugCheck(m_counter >= 0 && m_counter <= m_width*m_height);

	// Get ready for the next one, we might need several iterations
	// when the domain is not a power of two.
	if (m_counter < m_width*m_height) {
		for(;;) {
			++m_index;

			// Extract the interleaved locations
			uint x = 0, y = 0;
			for (int i = 0, bit = 1; i < 16; i++)
			{
				x |= (m_index & bit) >> (i);     bit <<= 1;
				y |= (m_index & bit) >> (i + 1); bit <<= 1;
			}

			// Found the valid location for the next one?
			if (x < m_width && y < m_height) {
				m_currentIndex.x = x;
				m_currentIndex.y = y;
				++m_counter;
				return true;
			}
		}
	}
	else {
		// We are done
		m_currentIndex.x = m_currentIndex.y = ~0x0;
		return false;
	}
}


// ------ Hilbert Iterator ------

HilbertIterator::HilbertIterator(uint w, uint h) :
	Iterator2D(w,h),
	m_currentIndex(0, 0),
	m_sideSize( nextPowerOfTwo(max(w,h)) )
{
	// Gets the initial level, we know our value is a power of 2
	static const uint b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 
                                 0xFF00FF00, 0xFFFF0000};;
	register uint level = (m_sideSize & b[0]) != 0;
	for (register int i = 4; i > 0; i--) // unroll for speed...
	{
	  level |= ((m_sideSize & b[i]) != 0) << i;
	}
	nvDebugCheck( (1<<level) == m_sideSize );

	m_level = level;
}

inline void HilbertIterator::matrix_times_vector2D(const int (&m)[2][2], const int (&v)[2], int (&dest)[2]) {
	dest[0] = m[0][0]*v[0] + m[0][1]*v[1];
	dest[1] = m[1][0]*v[0] + m[1][1]*v[1];
}

inline void HilbertIterator::matrix_multiply2D(const int (&op1)[2][2], const int (&op2)[2][2], int (&dest)[2][2]) {
	dest[0][0] = op1[0][0]*op2[0][0] + op1[0][1]*op2[1][0];
	dest[0][1] = op1[0][0]*op2[0][1] + op1[0][1]*op2[1][1];
	dest[1][0] = op1[1][0]*op2[0][0] + op1[1][1]*op2[1][0];
	dest[1][1] = op1[1][0]*op2[0][1] + op1[1][1]*op2[1][1];
}

inline void HilbertIterator::identityMatrix2D(int (&m)[2][2]) {
	m[0][0] = m[1][1] = 1;
	m[0][1] = m[1][0] = 0;
}

inline void HilbertIterator::matrix_copy2D(const int (&src)[2][2], int (&dest)[2][2]) {
	dest[0][0] = src[0][0];
	dest[0][1] = src[0][1];
	dest[1][0] = src[1][0];
	dest[1][1] = src[1][1];
}

inline void HilbertIterator::t_to_xy(int n, int t, int &x, int &y)
{
	int rt[2][2], rq[2][2], va[2], vb[2];

	identityMatrix2D(rt);
	x = y = 0;
	for (int k = n-1; k >= 0; --k) {
		const int j = 3 & (t >> (2*k));
		va[0] = 2*xbit2D[j] - 1;
		va[1] = 2*ybit2D[j] - 1;
		matrix_times_vector2D(rt, va, vb);
		x += ((vb[0] + 1)/2) << k;
		y += ((vb[1] + 1)/2) << k;
		matrix_copy2D(rt, rq);
		matrix_multiply2D(rq, R2D[j], rt);
	}
}


bool HilbertIterator::moveNext()
{
	nvDebugCheck(m_counter >= 0 && m_counter <= m_width*m_height);
	
	// Get ready for the next one, we might need several iterations
	// when the domain is not a power of two.
	if (m_counter < m_width*m_height) {
		for(;;) {
			int x,y;
			++m_index;

			// Get the next one
			t_to_xy(m_level, m_index, x, y);
			nvDebugCheck(x >= 0 && y >= 0);

			// Found a valid location for the next one?
			if ((uint)x < m_width && (uint)y < m_height) {
				
				m_currentIndex.x = x;
				m_currentIndex.y = y;
				++m_counter;
				return true;
			}
		}
	}
	else {
		// We are done
		m_currentIndex.x = m_currentIndex.y = ~0x0;
		return false;
	}

}

Index2D HilbertIterator::current() const
{
	nvDebugCheck(m_counter >= 0 && m_counter <= m_width*m_height);
	return m_currentIndex;
}
