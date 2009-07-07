// Copyright NVIDIA Corporation 2007 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_MAYAMESHBUILDER_H
#define NV_MESH_MAYAMESHBUILDER_H

#include <nvmesh/MeshBuilder.h>

class MDagPath;

namespace nv
{
	class TriMesh;
	class QuadTriMesh;
	class PolyMesh;
	namespace HalfEdge { class Mesh; }

	/// Build options.
	struct MayaMeshBuilderOptions
	{
		MayaMeshBuilderOptions();
		MayaMeshBuilderOptions(const MayaMeshBuilderOptions & options);
		const MayaMeshBuilderOptions & operator=(const MayaMeshBuilderOptions & options);

		bool addTexcoords;
		bool addNormals;
	};


	/// Maya mesh builder.
	class MayaMeshBuilder
	{
	public:

		MayaMeshBuilder(const MayaMeshBuilderOptions & options);

		bool addScene();
		bool addSelection();
		bool addNode(const MDagPath & dagPath);

		void done();
		void reset();

		TriMesh * buildTriMesh() const;
		QuadTriMesh * buildQuadTriMesh() const;
		//PolyMesh * buildPolyMesh() const;

		HalfEdge::Mesh * buildHalfEdgeMesh() const;

		// Expose attribute indices of the unified vertex array.
		uint vertexCount() const;
		uint positionIndex(uint vertex) const;
		uint normalIndex(uint vertex) const;
		uint texcoordIndex(uint vertex) const;

	private:

		const MayaMeshBuilderOptions m_options;

		MeshBuilder m_builder;

	};

}


#endif // NV_MESH_MAYAMESHBUILDER_H
