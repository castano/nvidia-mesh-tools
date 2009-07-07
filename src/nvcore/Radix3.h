///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 *	Contains source code from the article "Radix Sort Revisited".
 *	\file		IceRevisitedRadix.h
 *	\author		Pierre Terdiman
 *	\date		April, 4, 2000
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Include Guard
#ifndef NV_CORE_RADIXSORT_H
#define NV_CORE_RADIXSORT_H

#include <nvcore/nvcore.h>
#include <nvcore/Containers.h>

namespace nv
{

	class NVCORE_CLASS RadixSort3
	{
		NV_FORBID_COPY(RadixSort3);
	public:
		// Constructor/Destructor
		RadixSort3();
		~RadixSort3();

		// Sorting methods
		RadixSort3 & sort(const uint32 * input, uint nb, bool signedValues = true);
		RadixSort3 & sort(const float* input, uint nb);

		// Helpers
		RadixSort3 & sort(const Array<int> & input);
		RadixSort3 & sort(const Array<uint> & input);
		RadixSort3 & sort(const Array<float> & input);

		//! Access to results. mRanks is a list of indices in sorted order, i.e. in the order you may further process your data
		inline /*const*/ uint32 * ranks() const	{ return mRanks; }

		//! mIndices2 gets trashed on calling the sort routine, but otherwise you can recycle it the way you want.
		inline uint32 * recyclable() const { return mRanks2; }

		// Stats
		//! Returns the total number of calls to the radix sorter.
		inline uint32 totalCalls() const { return mTotalCalls;	}
		
		//! Returns the number of eraly exits due to temporal coherence.
		inline uint32 hits() const { return mNbHits; }

		bool setRankBuffers(uint32 * ranks0, uint32 * ranks1);

	private:
		
		uint32 mCurrentSize;    //!< Current size of the indices list
		uint32 * mRanks;        //!< Two lists, swapped each pass
		uint32 * mRanks2;
		
		// Stats
		uint32 mTotalCalls;		//!< Total number of calls to the sort routine
		uint32 mNbHits;         //!< Number of early exits due to coherence
		
		// Stack-radix
		bool mDeleteRanks;
		
		// Internal methods
		void checkResize(uint32 nb);
		bool resize(uint32 nb);
	};


} // nv namespace

#endif // NV_CORE_RADIXSORT_H
