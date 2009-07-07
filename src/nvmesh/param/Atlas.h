// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_ATLAS_H
#define NV_MESH_ATLAS_H

#include <nvcore/Containers.h>
#include <nvcore/Ptr.h>
#include <nvmesh/nvmesh.h>


namespace nv
{
	namespace HalfEdge { class Mesh; }

	class Chart;
	
	/*
	Atlas atlas(mesh);
	atlas.computeSeamlessTextureAtlas();
	mesh = atlas.result();
	*/

	/// An atlas is a set of charts.
	class Atlas
	{
	public:
		
		Atlas(const HalfEdge::Mesh * mesh);
		~Atlas();

		uint chartCount() const;

		const Chart * chartAt(uint i) const;
		Chart * chartAt(uint i);

		enum SeamlessTextureAtlas {
			Tiles,
			GroupTiles,
			AdaptiveTiles,
			Default = Tiles
		};
		
		// Compute a trivial seamless texture
		bool computeSeamlessTextureAtlas(uint w = 1024, uint h = 1024, SeamlessTextureAtlas mode=Default);
		
		// Extract the charts of the input mesh.
		void extractCharts();

		// Add chart scaling and packing.
		
		HalfEdge::Mesh * mergeCharts() const;

	private:

		bool staTiles();
		bool staGTiles();
		bool staATiles();

	private:
		
		const HalfEdge::Mesh * m_mesh;

		Array<Chart *> m_chartArray;
		
	};

	
	/// A chart is a connected set of faces with a certain topology (usually a disk).
	class Chart
	{
	public:
		
		Chart();
	
		void build(const Array<uint> & faceArray, const HalfEdge::Mesh * originalMesh);

		void closeHoles();

		bool isDisk() const;

		uint vertexCount() const { return m_mesh->vertexCount(); }
		uint faceCount() const { return m_faceArray.count(); }

		const HalfEdge::Mesh * mesh() const { return m_mesh.ptr(); }
		HalfEdge::Mesh * mesh() { return m_mesh.ptr(); }

		uint vertexIndex(uint i) const { return m_vertexIndexArray[i]; }
		
	private:
		
		// Chart mesh.
		AutoPtr<HalfEdge::Mesh> m_mesh;
		
		// List of faces of the original mesh that belong to this chart.
		Array<uint> m_faceArray;
		
		// Map vertices of the chart mesh to vertices of the original mesh.
		Array<uint> m_vertexIndexArray;
	};
	
	
} // nv namespace

#endif // NV_MESH_ATLAS_H
