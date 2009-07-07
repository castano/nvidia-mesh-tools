// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#include <nvmath/Triangle.h>

#include <nvmesh/TriMesh.h>
#include <nvmesh/geometry/Bounds.h>

#include "Grid.h"

// Fast voxel traversal (amanatides & woo)
// http://www.cse.yorku.ca/~amana/research/grid.pdf


using namespace nv;


GridBuilder::GridBuilder()
{
}


Grid * GridBuilder::build(const TriMesh & mesh, float cellSize)
{
	setupBounds(MeshBounds::box(&mesh), cellSize);

	//const uint zero = 0;
	const uint cellCount = m_xDim * m_yDim * m_zDim;
	m_cellCountArray.resize(cellCount, 0);
	m_cellOffsetArray.resize(cellCount);

	uint insertions = 0;

	// Compute Cell counts.
	const uint faceCount = mesh.faceCount();
	for (uint f = 0; f < faceCount; f++)
	{
		const TriMesh::Face & face = mesh.faceAt(f);
		
		const Vector3 v0 = mesh.vertexAt(face.v[0]).pos;
		const Vector3 v1 = mesh.vertexAt(face.v[1]).pos;
		const Vector3 v2 = mesh.vertexAt(face.v[2]).pos;
		
		insertions += allocateTriangle(Triangle(v0, v1, v2));
	}


	// Compute Cell offsets.
	m_cellOffsetArray[0] = 0;
	for (uint i = 1; i < cellCount; i++)
	{
		m_cellOffsetArray[i] = m_cellOffsetArray[i - 1] + m_cellCountArray[i - 1];
	}

	// Reset counts.
	foreach(i, m_cellCountArray)
	{
		m_cellCountArray[i] = 0;
	}


	// Build triangle indices.
	m_triangleArray.resize(insertions);

	for (uint f = 0; f < faceCount; f++)
	{
		const TriMesh::Face & face = mesh.faceAt(f);
		
		const Vector3 v0 = mesh.vertexAt(face.v[0]).pos;
		const Vector3 v1 = mesh.vertexAt(face.v[1]).pos;
		const Vector3 v2 = mesh.vertexAt(face.v[2]).pos;
		
		insertTriangle(f, Triangle(v0, v1, v2));
	}
	
	return NULL;
}


void GridBuilder::setupBounds(const Box & bounds, float cellSize)
{
	float xDim = ceilf(bounds.extents().x() / cellSize);
	float yDim = ceilf(bounds.extents().y() / cellSize);
	float zDim = ceilf(bounds.extents().z() / cellSize);

	m_gridBounds.setCenterExtents(bounds.center(), Vector3(xDim, yDim, zDim) * cellSize);
	m_cellSize = cellSize;
	m_xDim = (uint)xDim;
	m_yDim = (uint)yDim;
	m_zDim = (uint)zDim;
}


uint GridBuilder::allocateTriangle(const Triangle & tri)
{
	Box triBounds = tri.bounds();

	Vector3 minCorner = triBounds.minCorner() - m_gridBounds.minCorner();
	Vector3 maxCorner = triBounds.maxCorner() - m_gridBounds.minCorner();

	uint minZ = (uint)floorf(minCorner.z());
	uint maxZ = (uint)ceilf(maxCorner.z());
	uint minY = (uint)floorf(minCorner.y());
	uint maxY = (uint)ceilf(maxCorner.y());
	uint minX = (uint)floorf(minCorner.x());
	uint maxX = (uint)ceilf(maxCorner.x());

	uint insertions = 0;

	for (uint z = minZ; z <= maxZ; z++)
	{
		for (uint y = minY; y <= maxY; y++)
		{
			for (uint x = minX; x <= maxX; x++)
			{
				const uint idx = ((z * m_yDim) + y) * m_xDim + x;
				
				if (overlapNoBounds(tri, cellBounds(x, y, z)))
				{
					// Insert triangle.
					m_cellCountArray[idx]++;
					insertions++;
				}
			}
		}
	}

	return insertions;
}


void GridBuilder::insertTriangle(uint f, const Triangle & tri)
{
	Box triBounds = tri.bounds();

	Vector3 minCorner = triBounds.minCorner() - m_gridBounds.minCorner();
	Vector3 maxCorner = triBounds.maxCorner() - m_gridBounds.minCorner();

	uint minZ = (uint)floorf(minCorner.z());
	uint maxZ = (uint)ceilf(maxCorner.z());
	uint minY = (uint)floorf(minCorner.y());
	uint maxY = (uint)ceilf(maxCorner.y());
	uint minX = (uint)floorf(minCorner.x());
	uint maxX = (uint)ceilf(maxCorner.x());

	uint insertions = 0;

	for (uint z = minZ; z <= maxZ; z++)
	{
		for (uint y = minY; y <= maxY; y++)
		{
			for (uint x = minX; x <= maxX; x++)
			{
				const uint idx = ((z * m_yDim) + y) * m_xDim + x;
				
				if (overlapNoBounds(tri, cellBounds(x, y, z)))
				{
					// Insert triangle.
					m_triangleArray[m_cellOffsetArray[idx] + m_cellCountArray[idx]] = f;
					m_cellCountArray[idx]++;
				}
			}
		}
	}
}


Box GridBuilder::cellBounds(uint x, uint y, uint z)
{
	// @@ TODO
	return Box();
}


