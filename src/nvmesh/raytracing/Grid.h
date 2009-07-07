// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#ifndef NV_MESH_GRID_H
#define NV_MESH_GRID_H

#include <nvcore/Ptr.h>
#include <nvcore/RefCounted.h>
#include <nvcore/Containers.h>

#include <nvmath/Box.h>

#include <nvmesh/nvmesh.h>


namespace nv
{
	class TriMesh;
	class QuadTriMesh;
	class FaceBuffer;

	NV_DECLARE_PTR(Grid);
	
	// Space grid.
	class Grid : public RefCounted
	{
		NV_FORBID_COPY(Grid);
	public:

		Grid(Vector3 extents, float size);

		void addMesh(const TriMesh & mesh);

	private:
		
		FaceBuffer * m_faceArray;
		
	};
	

	class GridBuilder
	{
		NV_FORBID_COPY(GridBuilder);
	public:

		GridBuilder();
		
		Grid * build(const TriMesh & mesh, float cellSize);

	private:

		void setupBounds(const Box & bounds, float cellSize);

		uint allocateTriangle(const Triangle & tri);
		void insertTriangle(uint f, const Triangle & tri);
		
		Box cellBounds(uint x, uint y, uint z);

		Box m_gridBounds;
		float m_cellSize;
		uint m_xDim, m_yDim, m_zDim;

		Array<uint> m_cellCountArray;
		Array<uint> m_cellOffsetArray;
		Array<uint> m_triangleArray;
	};
	
} // nv namespace


#endif // NV_MESH_GRID_H
