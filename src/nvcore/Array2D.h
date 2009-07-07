// Copyright NVIDIA Corporation 2008 
// -- Edgar Velazquez-Armendariz <edgarv@nvidia.com>

#ifndef NV_BAKER_ARRAY2D_H
#define NV_BAKER_ARRAY2D_H

#include <nvcore/Containers.h>

namespace nv {

	/// Helper struct to provide 2D addressing using convenient idioms:
	/// (row,column), (y,x), (i,j) refer all to the same element
	struct Index2D {
		union {
			/* row, y, and i refer to the same */
			uint row;
			uint y;
			uint i;
		};
		union {
			/* column, x and j refer to the same */
			uint column;
			uint x;
			uint j;
		};

		Index2D() : row(0), column(0) {}

		Index2D(uint row, uint column) : row(row), column(column) {}

		inline uint index(uint stride) const {
			nvDebugCheck( (((uint64)row)*((uint64)stride)+((uint64)column)) 
				<= 0xFFFFFFFFLL);
			return row*stride + column;
		}
	};

	/**
	 * A simple wrapper for 2D addressing of arrays. The size of the 2D array 
	 * is fixed and it is determined upon its creation.
	 */
	template <class T>
	class Array2D {

	protected:
		const uint m_rows;
		const uint m_columns;
		Array<T> m_buffer;

		// Linearizes the row and column index so that we can actually access 
		// the internal buffer
		inline uint index(uint row, uint column) const {
			nvDebugCheck(row >= 0 && row < m_rows && column >= 0 && 
				column < m_columns);
			return row*m_columns + column;
		}

	public:

		/// Creates the 2D array allocating the required space. The total 
		/// number of elements in the array must be less than or equal 
		/// to 0xFFFFFFFF and greater than zero.
		Array2D(uint rows, uint columns) :
		  m_rows(rows), m_columns(columns)
		{
			nvDebugCheck(rows*columns > 0);
			nvDebugCheck(((uint64)rows)*((uint64)columns) <= 0xFFFFFFFFLL);
			m_buffer.resize(rows*columns);
		}

		/// Creates the 2D array allocating the required space and initializing 
		/// the elements with the given value. The total number of elements in 
		/// the array must be less than or equal to 0xFFFFFFFF and greater 
		/// than zero.
		Array2D(uint rows, uint columns, const T & elem) :
		  m_rows(rows), m_columns(columns)
		{
			nvDebugCheck(rows*columns > 0);
			nvDebugCheck(((uint64)rows)*((uint64)columns) <= 0xFFFFFFFFLL);
			m_buffer.resize(rows*columns, elem);
		}

		virtual ~Array2D() {}

		/// Mutable reference access
		T& operator() (uint row, uint column) {
			return m_buffer[index(row,column)];
		}
		T& operator[] (const Index2D & index) {
			return operator()(index.row, index.column);
		}
		T& get (uint row, uint column) {
			return operator()(row, column);
		}
		T& get (const Index2D & index) {
			return operator()(index.row, index.column);
		}

		/// Constant reference access
		const T& operator() (uint row, uint column) const {
			return m_buffer[index(row,column)];
		}
		const T& operator[] (const Index2D & index) const {
			return operator()(index.row, index.column);
		}
		const T& get (uint row, uint column) const {
			return operator()(row, column);
		}
		const T& get (const Index2D & index) const {
			return operator()(index.row, index.column);
		}

		/// Accessor methods
		uint rows() const    { return m_rows; }
		uint columns() const { return m_columns; }
		uint size() const    { return m_buffer.size(); }

	};


	/**
	 * Iterator interface for a matrix-like 2D domain. The iterators simply
	 * provide indices which in turn may be further used.
	 * The iterators use Java's iterators semantics: initially they are
	 * located right before the first element.
	 */
	struct Iterator2D
	{
		Iterator2D(uint w, uint h) : m_width(w), m_height(h), 
			m_index(-1), m_counter(0) {
			nvDebugCheck(((uint64)w)*((uint64)h) <= 0xFFFFFFFFLL);
		}

		virtual ~Iterator2D() {}

		/// next() and hasNext() are Java-style iterators

		/// Returns true if the iteration has more indices, so that
		/// next() would actually return something valid.
		virtual bool hasNext() const {
			nvDebugCheck(m_counter >= 0 && m_counter <= m_width*m_height);
			return m_counter < (m_width*m_height);
		}

		/// Returns the next index. If hasNext() was not true before
		/// calling this method its results are undefined.
		virtual Index2D next() {
			nvDebugCheck(hasNext());
			moveNext();
			return current();
		}

		/// current() and moveNext() are .NET-style iterators

		/// Returns the current index. If this is called before next()
		/// or moveNext() the result is undetermined.
		virtual Index2D current() const = 0;

		/// Advances one position. Returns false if it moves beyond
		/// the domain.
		virtual bool moveNext() {
			nvDebugCheck(m_index < (m_width*m_height));
			nvDebugCheck(m_counter <= (m_width*m_height));
			++m_index;
			++m_counter;
			nvDebugCheck((m_index < (m_width*m_height)) ==
				(m_counter <= (m_width*m_height)) );
			return m_index < (m_width*m_height);
		}

		uint width()  const { return m_width; }
		uint height() const { return m_height; }

	protected:
		const uint m_width;
		const uint m_height;
		int64 m_index;
		uint  m_counter;
	};

	/**
	 * Iterator that advances in scanline order: from left to right
	 * and from top to bottom (row major). The maximum number of
	 * elements for this iterator is 0xFFFFFFFE.
	 */
	struct ScanlineIterator : public Iterator2D 
	{
		ScanlineIterator(uint w, uint h) : Iterator2D(w,h) {
			nvDebugCheck(w*h < 0xFFFFFFFF);
		}

		virtual Index2D current() const {
			nvDebugCheck(m_index >= 0 && m_index < m_width * m_height);
			return Index2D((uint)m_index / m_width, (uint)m_index % m_width);
		}
	};

	/**
	 * Iterator that advances in column major order: from top to bottom
	 * and from left to right. The maximum number of
	 * elements for this iterator is 0xFFFFFFFE.
	 */
	struct ColumnIterator : public ScanlineIterator {

		ColumnIterator(uint w, uint h) : ScanlineIterator(w,h) {}

		virtual Index2D current() const {
			nvDebugCheck(m_index >= 0 && m_index < m_width * m_height);
			return Index2D((uint)m_index % m_height, (uint)m_index / m_height);
		}
	};

	/**
	 * Iterator which uses the Morton (Z) order. The maximum size
	 * of the domain is 64k-by-64k inclusive (0x10000-by-0x10000).
	 */
	struct MortonIterator : public Iterator2D {

		MortonIterator(uint w, uint h);

		virtual Index2D current() const;
		virtual bool moveNext();

	private:
		Index2D m_currentIndex;
	};

	/**
	 * Iterator which uses the Peano-Hilbert order. The maximum size
	 * of the domain is 64k-by-64k inclusive (0x10000-by-0x10000).
	 */
	struct HilbertIterator : public Iterator2D {

		HilbertIterator(uint w, uint h);
		
		virtual Index2D current() const;
		virtual bool moveNext();

	private:
		// The base size for constructing the object, since
		// the Hilbert fractal assumes a power of two domain.
		// This must be at most 0x10000)
		const uint m_sideSize;

		// Level: this is the log2 of the sideSize
		uint m_level;

		// The current index to return
		Index2D m_currentIndex;

		// Using method from
		// MAX, Nelson. "Visualizing Hilbert Curves", IEEE Visualization 1998
		// http://portal.acm.org/citation.cfm?id=288348
		static const int R2D[4][2][2];
		static const int xbit2D[4];
		static const int ybit2D[4];

		static void matrix_times_vector2D(const int (&m)[2][2], 
			const int (&v)[2], int (&dest)[2]);
		static void matrix_multiply2D(const int (&op1)[2][2], 
			const int (&op2)[2][2], int (&dest)[2][2]);
		static void identityMatrix2D(int (&m)[2][2]);
		static void matrix_copy2D(const int (&src)[2][2], int (&dest)[2][2]);
		static void t_to_xy(/*order*/int n, /*value*/int t, int &x, int &y);
	};


}	// nv namespace

#endif // NV_BAKER_ARRAY2D_H
