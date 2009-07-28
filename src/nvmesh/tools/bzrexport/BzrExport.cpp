// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmesh/halfedge/HalfEdge.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>

#include <nvmesh/import/MeshImport.h>

#include <nvmesh/MeshTopology.h>
#include <nvmesh/geometry/MeshNormals.h>
#include <nvmesh/geometry/MeshTransform.h>

#include <nvmesh/subdiv/Subdivide.h>
#include <nvmesh/subdiv/AccMeshBuilder.h>
#include <nvmesh/subdiv/AccMesh.h>
#include <nvmesh/subdiv/AccPatch.h>
#include <nvmesh/subdiv/RemapFaces.h>

#include <nvmesh/param/Atlas.h>
#include <nvmesh/param/Seams.h>

#include <nvmath/Vector.h>
#include <nvmath/TypeSerialization.h>

//#include <nvcore/FileSystem.h>
#include <nvcore/StdStream.h>
#include <nvcore/TextWriter.h>
#include <nvcore/Ptr.h>

#include <stdio.h> // printf

using namespace nv;

uint controlPointOrderRegular[16] = {
    0, 1, 4, 5,
    3, 7, 2, 6,
    15, 14, 11, 10,
    12, 8, 13, 9
};

uint controlPointOrderGregory[20] = {
    8, 9, 12, 1, 0,
    11, 13, 10, 3, 2,
    19, 18, 15, 6, 7,
    16, 14, 17, 4, 5
};

const bool s_outputControlPoints = true;
const bool s_outputStencils = false;
const bool s_outputMesh = true;
const bool s_outputTwoStageData = true;


static void OutputHeader(StdOutputStream & stream, const AccMesh * accMesh, const HalfEdge::Mesh * mesh, bool exportText)
{
	if (!exportText)
	{
		uint header = uint('B') | (uint('Z') << 8) | (uint('R') << 16) | (uint(' ') << 24);
		uint version = 0x0100; // 1.0

		stream << header << version;
	}

	//
	uint regularPatchCount = accMesh->regularPatchCount();
	uint quadPatchCount = accMesh->quadPatchCount();
	uint triPatchCount = accMesh->triPatchCount();

	TextWriter writer(&stream);

	if (exportText) {
		writer << "const int regular_patch_count = " << regularPatchCount << ";\n";
		writer << "const int quad_patch_count = " << quadPatchCount << ";\n";
		writer << "const int tri_patch_count = " << triPatchCount << ";\n";
	}
	else {
		stream << regularPatchCount;
		stream << quadPatchCount;
		stream << triPatchCount;
	}
}

static void OutputControlPoints(StdOutputStream & stream, const AccMesh * accMesh, const HalfEdge::Mesh * mesh, bool reorder, bool exportText)
{
	const uint regularPatchCount = accMesh->regularPatchCount();
	const uint quadPatchCount = accMesh->quadPatchCount();
	const uint triPatchCount = accMesh->triPatchCount();

    TextWriter writer(&stream);

	// Output regular patches.
	if (exportText) writer << "\n// Regular Patches\n\n";

	if (regularPatchCount > 0)
	{
		if (exportText) writer << "const float3 regular_bezier_control_points[regular_patch_count][16] =\n{\n";

		// Bezier control points.
		for (uint p = 0; p < regularPatchCount; p++)
		{
			int idx = p;
			if (reorder) idx = controlPointOrderRegular[p];

			//const BezierAccPatch * patch = accMesh->regularPatchAt(idx).bezierAccPatch;
			const PmRegularAccPatch * patch = static_cast<const PmRegularAccPatch *>(accMesh->regularPatchAt(p).pmAccPatch);
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 16; i++)
			{
				Vector3 v = patch->position(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float2 regular_tex_coords[regular_patch_count][16] =\n{\n";
		}

		// Texcoords.
		for (uint p = 0; p < regularPatchCount; p++)
		{
			const TexCoordPatch * patch = accMesh->regularPatchAt(p).texCoordPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 16; i++)
			{
				Vector2 v = patch->texCoord(i / 4, i % 4);
				if (exportText) writer.write("{ %f, %f }, ", v.x(), v.y());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float3 regular_normals[regular_patch_count][16] =\n{\n";
		}

		// Normals.
		for (uint p = 0; p < regularPatchCount; p++)
		{
			const BezierAccPatch * patch = accMesh->regularPatchAt(p).bezierAccPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 16; i++)
			{
				Vector3 v = patch->normal(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) writer << "};\n\n";
	}
	else
	{
		if (exportText) {
			writer << "float3 regular_bezier_control_points[1][16];\n";
			writer << "float2 regular_tex_coords[1][16];\n\n";
			writer << "float3 regular_normals[1][16];\n";
		}
	}

	// Output quad patches.
	if (exportText) writer << "\n// Quad Patches\n\n";

	if (quadPatchCount > 0)
	{
		if (exportText) writer << "const float3 quad_bezier_control_points[quad_patch_count][32] =\n{\n";

		// Bezier control points.
		for (uint p = 0; p < quadPatchCount; p++)
		{
			const BezierAccPatch * patch = accMesh->quadPatchAt(p).bezierAccPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 16; i++)
			{
				Vector3 v = patch->position(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}
			for (uint i = 0; i < 8; i++)
			{
				Vector3 v = patch->tangent((i < 4) ? i : i + 4);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}
			for (uint i = 0; i < 8; i++)
			{
				Vector3 v = patch->bitangent((i < 4) ? i : i + 4);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float3 quad_gregory_control_points[quad_patch_count][20] =\n{\n";
		}

		// Gregory position control points.
		for (uint p = 0; p < quadPatchCount; p++)
		{
			int idx = p;
			if (reorder) idx = controlPointOrderGregory[p];

			const GregoryAccPatch * patch = accMesh->quadPatchAt(idx).gregoryAccPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 20; i++)
			{
				Vector3 v = patch->position(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float3 quad_pm_control_points[quad_patch_count][24] =\n{\n";
		}

		// PmQuad position control points.
		for (uint p = 0; p < quadPatchCount; p++)
		{
			const PmQuadAccPatch * patch = static_cast<const PmQuadAccPatch *>(accMesh->quadPatchAt(p).pmAccPatch);
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 24; i++)
			{
				Vector3 v = patch->position(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float2 quad_tex_coords[quad_patch_count][16] =\n{\n";
		}

		for (uint p = 0; p < quadPatchCount; p++)
		{
			const TexCoordPatch * patch = accMesh->quadPatchAt(p).texCoordPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 16; i++)
			{
				Vector2 v = patch->texCoord(i / 4, i % 4);
				if (exportText) writer.write("{ %f, %f }, ", v.x(), v.y());
				else stream << v;
			}
			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float3 quad_normals[quad_patch_count][16] =\n{\n";
		}

		// Normals.
		for (uint p = 0; p < quadPatchCount; p++)
		{
			const BezierAccPatch * patch = accMesh->quadPatchAt(p).bezierAccPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 16; i++)
			{
				Vector3 v = patch->normal(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) writer << "};\n\n";
	}
	else
	{
		if (exportText) {
			writer << "float3 quad_bezier_control_points[1][32];\n";
			writer << "float3 quad_gregory_control_points[1][20];\n";
			writer << "float3 quad_pm_control_points[1][24];\n";
			writer << "float2 quad_tex_coords[1][16];\n\n";
			writer << "float3 quad_normals[1][16];\n";
		}
	}


	// Output triangle patches.
	if (exportText) writer << "\n// Triangle Patches\n\n";

	if (triPatchCount > 0)
	{
		if (exportText) writer << "const float3 tri_gregory_control_points[tri_patch_count][15] =\n{\n";

		// Gregory position control points.
		for (uint p = 0; p < triPatchCount; p++)
		{
			const GregoryAccPatch * patch = accMesh->triPatchAt(p).gregoryAccPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 15; i++)
			{
				Vector3 v = patch->position(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float3 tri_pm_control_points[tri_patch_count][19] =\n{\n";
		}

		// PmTriangle position control points.
		for (uint p = 0; p < triPatchCount; p++)
		{
			const PmTriangleAccPatch * patch = static_cast<const PmTriangleAccPatch *>(accMesh->triPatchAt(p).pmAccPatch);
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 19; i++)
			{
				Vector3 v = patch->position(i);
				if (exportText) writer.write("{ %f, %f, %f }, ", v.x(), v.y(), v.z());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) {
			writer << "};\n\n";
			writer << "const float2 tri_tex_coords[tri_patch_count][12] =\n{\n";
		}

		// Texcoords.
		for (uint p = 0; p < triPatchCount; p++)
		{
			const TexCoordPatch * patch = accMesh->triPatchAt(p).texCoordPatch;
			nvCheck (patch != NULL);

			if (exportText) writer << "\t{ ";

			for (uint i = 0; i < 12; i++)
			{
				Vector2 v = patch->texCoord(i / 4, i % 4);
				if (exportText) writer.write("{ %f, %f }, ", v.x(), v.y());
				else stream << v;
			}

			if (exportText) writer << "},\n";
		}

		if (exportText) writer << "};\n\n";
	}
	else
	{
		if (exportText) {
			writer << "float3 tri_gregory_control_points[1][15];\n";
			writer << "float3 tri_pm_control_points[1][19];\n";
			writer << "float2 tri_tex_coords[1][12];\n\n";
		}
	}
}

static void OutputStencils(StdOutputStream & stream, const AccMesh * accMesh, const HalfEdge::Mesh * mesh, int maxVertexCount, bool exportText)
{
	const uint regularPatchCount = accMesh->regularPatchCount();
	const uint quadPatchCount = accMesh->quadPatchCount();
	const uint triPatchCount = accMesh->triPatchCount();

	TextWriter writer(&stream);

	// Stencils are only exported as text so far.
	if (exportText)
	{
		writer << "\n// Face Topologies, stencils\n\n";

		uint faceTopologyCount = accMesh->faceTopologyCount();

		writer << "const int face_topology_count = " << faceTopologyCount << "; \n\n";

		uint primitiveSize = 0;

		const uint patchCount = accMesh->patchCount();
		for	(uint p = 0; p < patchCount; p++)
		{
			const FaceOneRing * faceOneRing = accMesh->patchAt(p).faceOneRing;

			uint patchVertexCount = faceOneRing->vertexCount();
			if (patchVertexCount > 32)
			{
				printf("Warning: patch %u has %u vertices.\n", p, patchVertexCount);
				patchVertexCount = 32;
			}

			primitiveSize = max(primitiveSize, patchVertexCount);
		}

		if (primitiveSize < 24) primitiveSize = 24;

		printf("Primitive Size = %u\n", primitiveSize);

		writer << "const int primitive_size = " << primitiveSize << ";\n\n";

		writer << "const float bezier_stencils[face_topology_count * 32 * primitive_size] =\n{\n";

		for (uint t = 0; t < faceTopologyCount; t++)
		{
			const FaceTopology * faceTopology = accMesh->faceTopologyAt(t);
			nvDebugCheck(faceTopology != NULL);

			float stencilTable[32*32];
			memset(stencilTable, 0, sizeof(stencilTable));

			if (faceTopology->bezierAccStencil != NULL)
			{
				for (uint i = 0; i < 32; i++)
				{
					const StencilMask * stencil;

					if (i < 16) stencil = & faceTopology->bezierAccStencil->positionStencil(i);
					else if (i < 24) stencil = & faceTopology->bezierAccStencil->tangentStencil((i-16 < 4) ? i-16 : i-16 + 4);
					else stencil = & faceTopology->bezierAccStencil->bitangentStencil((i-24 < 4) ? i-24 : i-24 + 4);

					for (uint v = 0; v < stencil->count(); v++)
					{
						stencilTable[32 * v + i] = (*stencil)[v];
					}
				}
			}

			writer.write("\t// %d (%.8X)\n", faceTopology->index, faceTopology->topologyId);

			for (uint v = 0; v < primitiveSize; v++)
			{
				writer << "\t";

				for (uint i = 0; i < 32; i++)
				{
					writer << stencilTable[32 * v + i] << ", ";
				}

				writer << "\n";
			}
		}

		writer << "};\n\n";

		writer << "const float gregory_stencils[face_topology_count * 20 * primitive_size] =\n{\n";

		for (uint t = 0; t < faceTopologyCount; t++)
		{
			const FaceTopology * faceTopology = accMesh->faceTopologyAt(t);
			nvDebugCheck(faceTopology != NULL);

			float stencilTable[20*32];
			memset(stencilTable, 0, sizeof(stencilTable));

			if (faceTopology->gregoryAccStencil != NULL)
			{
				for (uint i = 0; i < 20; i++)
				{
					const StencilMask & stencil = faceTopology->gregoryAccStencil->positionStencil(i);

					for (uint v = 0; v < stencil.count(); v++)
					{
						stencilTable[20 * v + i] = stencil[v];
					}
				}
			}

			writer.write("\t// %d (%.8X)\n", faceTopology->index, faceTopology->topologyId);

			for (uint v = 0; v < primitiveSize; v++)
			{
				writer << "\t";

				for (uint i = 0; i < 20; i++)
				{
					writer << stencilTable[20 * v + i] << ", ";
				}

				writer << "\n";
			}
		}

		writer << "};\n\n";

		writer << "\n// Face Indices\n\n";

		// Output regular face indices
		if (regularPatchCount > 0)
		{
			writer << "const unsigned int regular_face_vertex_count[regular_patch_count] =\n{\n";

			for (uint p = 0; p < regularPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->regularPatchAt(p).faceOneRing;

				const uint vertexCount = face->vertexCount();
				nvCheck (vertexCount <= 16);

				writer << "\t" << vertexCount << ",\n";
			}

			writer << "};\n\n";

			writer << "const unsigned int regular_face_topology_index[regular_patch_count] =\n{\n";

			for (uint p = 0; p < regularPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->regularPatchAt(p).faceOneRing;

				const uint index = face->faceTopology()->index;
				nvCheck (index < faceTopologyCount);

				writer << "\t" << index << ",\n";
			}

			writer << "};\n\n";

			writer << "const unsigned short regular_vertex_indices[regular_patch_count * 16] =\n{\n";

			for (uint p = 0; p < regularPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->regularPatchAt(p).faceOneRing;

				const uint vertexCount = face->vertexCount();
				nvCheck (vertexCount <= 16);

				writer << "\t";

				uint i = 0, idx;
				for (; i < vertexCount; i++)
				{
					uint sorted_idx = face->vertexIndexAt(i);

					idx = face->vertexAt(sorted_idx)->id();

					writer << idx << ", ";
				}
				for (; i < 16; i++)
				{
					// Repeat the last index
					writer << idx << ", ";
				}

				writer << "\n";
			}

			writer << "};\n\n";


			writer << "const unsigned int regular_stencil_indices[regular_patch_count * 16] =\n{\n";

			for (uint p = 0; p < regularPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->regularPatchAt(p).faceOneRing;

				const uint vertexCount = face->vertexCount();
				nvCheck (vertexCount <= primitiveSize);

				writer << "\t";

				uint i = 0, idx;
				for (; i < vertexCount; i++)
				{
					idx = face->vertexIndexAt(i);
					writer << idx << ", ";
				}
				for (; i < 16; i++)
				{
					// Repeat the last index
					writer << idx << ", ";
				}

				writer << "\n";
			}

			writer << "};\n\n";
		}
		else
		{
			writer << "const unsigned int regular_face_vertex_count[1] = { 0 };\n";
			writer << "const unsigned int regular_face_topology_index[1] = { 0 };\n";
			writer << "const unsigned short regular_vertex_indices[1] = { 0 };\n";
			writer << "const unsigned int regular_stencil_indices[1] = { 0 };\n\n";
		}

		// Output quad face indices
		if (quadPatchCount > 0)
		{
			writer << "const unsigned int quad_face_vertex_count[quad_patch_count] =\n{\n";

			for (uint p = 0; p < quadPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->quadPatchAt(p).faceOneRing;

				const uint vertexCount = face->vertexCount();
				nvCheck (vertexCount <= 32);

				writer << "\t" << vertexCount << ",\n";
			}

			writer << "};\n\n";

			writer << "const unsigned int quad_face_topology_index[quad_patch_count] =\n{\n";

			for (uint p = 0; p < quadPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->quadPatchAt(p).faceOneRing;

				const uint index = face->faceTopology()->index;
				nvCheck (index < faceTopologyCount);

				writer << "\t" << index << ",\n";
			}

			writer << "};\n\n";

			uint quadPatchVertexCount = maxVertexCount;
			if (maxVertexCount == 0) quadPatchVertexCount = 24;

			writer << "const unsigned short quad_vertex_indices[quad_patch_count * " << quadPatchVertexCount << "] =\n{\n";

			for (uint p = 0; p < quadPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->quadPatchAt(p).faceOneRing;

				const uint vertexCount = face->vertexCount();
				nvCheck (vertexCount <= 32);

				writer << "\t";

				uint i = 0, idx;
				for (; i < vertexCount; i++)
				{
					uint sorted_idx = face->vertexIndexAt(i);

					idx = face->vertexAt(sorted_idx)->id();

					writer << idx << ", ";
				}
				for (; i < quadPatchVertexCount; i++)
				{
					// Repeat the last index
					writer << idx << ", ";
				}

				writer << "\n";
			}

			writer << "};\n\n";


			writer << "const unsigned int quad_stencil_indices[quad_patch_count * primitive_size] =\n{\n";

			for (uint p = 0; p < quadPatchCount; p++)
			{
				const FaceOneRing * face = accMesh->quadPatchAt(p).faceOneRing;

				const uint vertexCount = face->vertexCount();
				nvCheck (vertexCount <= primitiveSize);

				writer << "\t";

				uint i = 0, idx;
				for (; i < vertexCount; i++)
				{
					idx = face->vertexIndexAt(i);
					writer << idx << ", ";
				}
				for (; i < primitiveSize; i++)
				{
					// Repeat the last index
					writer << idx << ", ";
				}

				writer << "\n";
			}

			writer << "};\n\n";
		}
		else
		{
			writer << "const unsigned int quad_face_vertex_count[1] = { 0 };\n";
			writer << "const unsigned int quad_face_topology_index[1] = { 0 };\n";
			writer << "const unsigned short quad_vertex_indices[1] = { 0 };\n";
			writer << "const unsigned int quad_stencil_indices[1] = { 0 };\n\n";
		}
	}
}

static void OutputMesh(StdOutputStream & stream, const AccMesh * accMesh, const HalfEdge::Mesh * mesh, int maxValence, bool exportText)
{
    nvCheck(maxValence <= 16);

	const uint regularPatchCount = accMesh->regularPatchCount();
	const uint quadPatchCount = accMesh->quadPatchCount();
	const uint triPatchCount = accMesh->triPatchCount();

	const uint patchCount = accMesh->patchCount();
	uint vertexCount = mesh->vertexCount();

    TextWriter writer(&stream);

    if (exportText) {
        writer << "\n// Input Mesh\n\n";
        writer << "\n// Vertices\n\n";
        writer << "const int vertex_count = " << vertexCount << ";\n\n";
    }
    else {
        stream << vertexCount;
    }


    if (exportText) {
        writer << "const float3 vertices[vertex_count] =\n{\n";
    }

    // Output mesh vertices
    for (uint v = 0; v < vertexCount; v++)
    {
        const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
        Vector3 pos = vertex->pos();
        
        if (exportText) writer.write("\t{%f, %f, %f},\n", pos.x(), pos.y(), pos.z());
        else stream << pos;

    }
	if (exportText) writer << "};\n\n";

    if (exportText) {
        writer << "const int valences[vertex_count] =\n{\n";
    }

#pragma message(NV_FILE_LINE "TODO: Move this to OutputTwoStageData")

    // Output mesh valences
    for (uint v = 0; v < vertexCount; v++)
    {
        const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
        int n = vertex->valence();
        
        if (exportText) writer.write("\t%d,\n", n);
        else stream << n;
    }
    if (exportText) writer << "};\n\n";

	// Output maxvalenve
    if (exportText) {
        writer << "const int max_valence = " << maxValence << ";\n\n";
    }
    else {
        stream << maxValence;
    }
    

    // Output face indices
    if (exportText) writer << "\n// Face Indices\n\n";

	if (regularPatchCount > 0)
    {
        if (exportText) writer << "const int regular_face_indices[regular_patch_count][4] =\n{\n";

		for (uint p = 0; p < regularPatchCount; ++p)
        {
			const HalfEdge::Face * face = accMesh->regularPatchAt(p).faceOneRing->face();

			if (exportText) writer << "\t{ ";

			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
			{
				uint id = it.current()->vertex()->id();
				if (exportText) writer << id << ", ";
				else stream << id;
			}

			if (exportText) writer << "},\n";
        }

        if (exportText) writer << "};\n\n";
    }
    else
    {
        if (exportText) writer << "const int regular_face_indices[1][4] = {{0, 0, 0, 0}};\n\n";
    }

	if (quadPatchCount > 0)
    {
        if (exportText) writer << "const int quad_face_indices[quad_patch_count][4] =\n{\n";

		for (uint p = 0; p < quadPatchCount; ++p)
		{
			const HalfEdge::Face * face = accMesh->quadPatchAt(p).faceOneRing->face();

			if (exportText) writer << "\t{ ";

			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
			{
				uint id = it.current()->vertex()->id();
				if (exportText) writer << id << ", ";
				else stream << id;
			}

			if (exportText) writer << "},\n";
        }

        if (exportText) writer << "};\n\n";
    }
    else
    {
        if (exportText) writer << "const int quad_face_indices[1][4] = {{0, 0, 0, 0}};\n\n";
    }

	if (triPatchCount > 0)
    {
        if (exportText) writer << "const int tri_face_indices[tri_patch_count][3] =\n{\n";

		for (uint p = 0; p < triPatchCount; ++p)
		{
			const HalfEdge::Face * face = accMesh->triPatchAt(p).faceOneRing->face();

			if (exportText) writer << "\t{ ";

			for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance())
			{
				uint id = it.current()->vertex()->id();
				if (exportText) writer << id << ", ";
				else stream << id;
			}

			if (exportText) writer << "},\n";
		}

        if (exportText) writer << "};\n\n";
    }
    else
    {
        if (exportText) writer << "const int tri_face_indices[1][3]={{0,0,0}};\n\n";
    }
}


static void OutputTwoStageData(StdOutputStream & stream, const AccMesh * accMesh, const HalfEdge::Mesh * mesh, int maxValence, bool exportText)
{
	nvCheck(maxValence <= 16);

	const uint regularPatchCount = accMesh->regularPatchCount();
	const uint quadPatchCount = accMesh->quadPatchCount();
	const uint triPatchCount = accMesh->triPatchCount();

	uint faceCount = mesh->faceCount();
	uint vertexCount = mesh->vertexCount();

	TextWriter writer(&stream);

	// Part4: precomputed two-stage info
	if (exportText) writer << "\n// precomputed two-stage info\n\n";

	// Output vertex nbrs
	if (exportText) {
		writer << "const int vertex_neighbors[vertex_count][3 * max_valence] =\n{\n";
	}

	for (uint v = 0; v < vertexCount; v++)
	{
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);

		if (exportText) writer << "\t{ ";

		int i = 0;
		for (HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), i++)
		{
			const HalfEdge::Edge * edge = it.current();

			uint directNeighbor = edge->to()->id();
			// direct--> diag must be in ccw order
			uint diagonalNeighbor = edge->next()->to()->id();
			uint directNeighbor_p = edge->pair()->next()->to()->id();

			if (diagonalNeighbor == directNeighbor_p) {
				// triangle
				if (exportText) writer.write("%u, %u, %u, ", directNeighbor, directNeighbor, directNeighbor_p);
				else stream << directNeighbor << directNeighbor << directNeighbor_p;
			}
			else {
				// quad
				if (exportText) writer.write("%u, %u, %u, ", directNeighbor, diagonalNeighbor, diagonalNeighbor);
				else stream << directNeighbor << diagonalNeighbor << diagonalNeighbor;
			}

		}
		uint pad = ~0;
		for (; i < maxValence; i++)
		{
			if (!exportText) stream << pad << pad << pad;
			else {writer.write("%u, %u, %u, ", 0,0,0);}
		}

		if (exportText) writer << "},\n";
	}

	if (exportText) {
		writer << "};\n\n";
	}

	// Output face offsets.
	Array<int> regularFaceOffsets(8 * regularPatchCount);
	Array<int> quadFaceOffsets(8 * quadPatchCount);
	Array<int> triFaceOffsets(6 * triPatchCount);

	const uint patchCount = accMesh->patchCount();
	for (uint p = 0; p < patchCount; p++)
	{
		const HalfEdge::Face * face = accMesh->patchAt(p).faceOneRing->face();

		Array<int> * offsets = NULL;

		AccMesh::FaceType faceType = accMesh->patchAt(p).faceType;
		if (faceType == AccMesh::FaceType_Regular) offsets = &regularFaceOffsets;
		else if (faceType == AccMesh::FaceType_Quad) offsets = &quadFaceOffsets;
		else if (faceType == AccMesh::FaceType_Tri) offsets = &triFaceOffsets;

#pragma message(NV_FILE_LINE "TODO: Use packet offsets")
		//uint packedOffsets = 0;

		// For each vertex in the quad:
		int i = 0;
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++)
		{
			const HalfEdge::Edge * edge = it.current();

			const HalfEdge::Vertex * vertex = edge->from();
			const uint valence = vertex->valence();

			// Record location of the edges of this face with respect to this vertex.
			int j = 0;
			for (HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance(), j++)
			{
				if (edge == eit.current())
				{
					int k = (j + valence - 1) % valence;
					offsets->append(j);
					offsets->append(k);
					//packedOffsets |= j << 4 * (2 * i + 0);
					//packedOffsets |= k << 4 * (2 * i + 1);
					break;
				}
			}
		}

		// offsets->append(packetOffsets);
	}

	nvCheck(regularFaceOffsets.count() == 8 * regularPatchCount);
	nvCheck(quadFaceOffsets.count() == 8 * quadPatchCount);
	nvCheck(triFaceOffsets.count() == 6 * triPatchCount);

	if (exportText) writer << "\n// Face Offsets\n\n";

	if (regularPatchCount > 0)
	{
		if (exportText) writer << "const int regular_face_offsets[8 * regular_patch_count] =\n{\n\t";

		foreach(i, regularFaceOffsets) {
			if (exportText) writer << regularFaceOffsets[i] << ", ";
			else stream << regularFaceOffsets[i];
		}

		if (exportText) writer << "\n};\n\n";
	}
	else
	{
		if (exportText) writer << "const int regular_face_offsets[1]={0};\n\n";
	}

	if (quadPatchCount > 0)
	{
		if (exportText) writer << "const int quad_face_offsets[8 * quad_patch_count] =\n{\n\t";

		foreach(i, quadFaceOffsets) {
			if (exportText) writer << quadFaceOffsets[i] << ", ";
			else stream << quadFaceOffsets[i];
		}

		if (exportText) writer << "\n};\n\n";
	}
	else
	{
		if (exportText) writer << "const int quad_face_offsets[1]={0};\n\n";
	}

	if (triPatchCount > 0)
	{
		if (exportText) writer << "const int tri_face_offsets[6 * tri_patch_count] =\n{\n\t";

		foreach(i, triFaceOffsets) {
			if (exportText) writer << triFaceOffsets[i] << ", ";
			else stream << triFaceOffsets[i];
		}

		if (exportText) writer << "\n};\n\n";
	}
	else
	{
		if (exportText) writer << "const int tri_face_offsets[1] = { 0 };\n\n";
	}
}

struct MyMessageHandler : public MessageHandler
{
	void log(const char * str, va_list arg)
	{
		va_list tmp;
		va_copy(tmp, arg);
		vprintf(str, arg);

#if _DEBUG && NV_OS_WIN32
//		static StringBuilder buffer(1024);
//		buffer.format(str, arg);
//		OutputDebugStringA(buffer.str());
#endif

		va_end(tmp);
	}
};

int main(int argc, const char * argv[])
{
	MyMessageHandler messageHandler;
	debug::setMessageHandler(&messageHandler);

	const char * fileName = NULL;

    bool exportText = false;
    bool reorder = false;
	bool quadOnly = false;
	int maxVertexCount = 0;

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-text") == 0)
        {
            exportText = true;
        }
        else if (strcmp(argv[i], "-reorder") == 0)
        {
            reorder = true;
        }
        else if (strcmp(argv[i], "-quadOnly") == 0)
        {
            quadOnly = true;
        }
        else if (strcmp(argv[i], "-maxVertexCount") == 0)
        {
            maxVertexCount = atoi(argv[++i]);
        }
		else if (argv[i][0] != '-')
        {
            fileName = argv[i];
        }
    }

    if (fileName == NULL) {
        printf("Usage: nvbzrexport filename\n");
        return EXIT_SUCCESS;
    }

    AutoPtr<MeshImport> importer(MeshImport::importer(fileName));

    if (importer == NULL) {
        printf("Error, unkown file type '%s'\n", fileName);
        return EXIT_FAILURE;
    }

    /*if (!FileSystem::fileExists(fileName)) {
		printf("Error, file '%s' does not exist.\n", fileName);
		return EXIT_FAILURE;
    }*/

	StdInputStream inputStream(fileName);
	if (inputStream.isError()) {
		printf("Error, cannot open '%s'\n", fileName);
		return EXIT_FAILURE;
	}

	importer->import(&inputStream);

	AutoPtr<HalfEdge::Mesh> halfEdgeMesh(importer->builder().buildHalfEdgeMesh());

	RemapFaces::consistentPatchOrientation(halfEdgeMesh.ptr());
	//RemapFaces::minimizeTopologyCount(halfEdgeMesh.ptr());

	MeshNormals::computeCatmullClarkNormals(halfEdgeMesh.ptr());

	// Validate input mesh.
	uint maxValence = 16;
	uint valence = 0;

	for (HalfEdge::Mesh::ConstVertexIterator v(halfEdgeMesh.ptr()); !v.isDone(); v.advance())
	{
		const uint currentValence = v.current()->valence();
		if (currentValence > maxValence)
		{
			printf("Error, vertex valence (%d) exceeds maximum (%d).\n", currentValence, maxValence);
			return EXIT_FAILURE;
		}

		if (currentValence > valence)
		   valence = currentValence;
	}
	maxValence = valence;

	for (HalfEdge::Mesh::ConstFaceIterator f(halfEdgeMesh->faces()); !f.isDone(); f.advance())
	{
		const HalfEdge::Face * face = f.current();

		const uint edgeCount = face->edgeCount();
		if (edgeCount != 3 && edgeCount != 4)
		{
			printf("Error, only triangle and quad meshes supported.\n");
			return EXIT_FAILURE;
		}
	}

	// Build ACC mesh.
	AccMeshBuilder builder(halfEdgeMesh.ptr(), BoundaryMode_Spline);

	uint flags = BuildFlags_BuildAll;
	if (quadOnly) flags |= BuildFlags_GenerateQuadPatches;
	else flags |= BuildFlags_GenerateAll;

	AutoPtr<AccMesh> accMesh(builder.buildAccMesh(flags));

//	AccMesh::Sort order[] = { AccMesh::Sort_ByTopology, AccMesh::Sort_ByVertexCount };
//	accMesh->sort(order, sizeof(order) / sizeof(order[0]));
//	accMesh->optimize(AccMesh::Granularity_Single); // Optimize without respecting the topology order or the vertex counts.

    // Create bzr file.
	Path outputFileName(fileName);
    outputFileName.stripExtension();
    outputFileName.append(".bzr"); 

    if (exportText) {
        outputFileName.append(".h");
    }

    StdOutputStream outputStream(outputFileName);

    // Bezier File Format :
    //   Header ('BZR ')            | sizeof(uint)
    //   Version (1.0)              | sizeof(uint)
//	//   Flags	                    | sizeof(uint)
	//   Regular patch count        | sizeof(uint)
	//   Quad patch count           | sizeof(uint)
	//   Triangle patch count       | sizeof(uint)

	//   Part 1.  Precomputed Control Points:
    //     Regular Patches:
    //       Bezier control points        | 16 * regularPatchCount * sizeof(float3)
    //       Texture coordinates          | 16 * regularPatchCount * sizeof(float2)
	//       Normal control points        | 16 * regularPatchCount * sizeof(float3)
    //     Quad Patches:
    //       Bezier control points        | 32 * quadPatchCount * sizeof(float3)
    //       Gregory control points       | 20 * quadPatchCount * sizeof(float3)
    //       Pm control points            | 24 * quadPatchCount * sizeof(float3)
    //       Texture coordinates          | 16 * quadPatchCount * sizeof(float2)
	//       Normal control points        | 16 * quadPatchCount * sizeof(float3)
    //     Triangle Patches:
    //       Gregory control points       | 15 * trianglePatchCount * sizeof(float3)
    //       Pm control points            | 19 * trianglePatchCount * sizeof(float3)
    //       Texture coordinates          | 12 * trianglePatchCount * sizeof(float2)

	//   Part 2. Stencils:

	//   Part 3. Input Mesh Topology:
    //     Vertex count                   | sizeof(uint)
    //     Vertices                       | vertexCount * sizeof(float3)
    //     Valences                       | vertexCount * sizeof(int)
    //     Max valence                    | sizeof(uint)
    //     Regular face indices           | 4 * regularPatchCount * sizeof(uint)
    //     Quad face indices              | 4 * irregularpatchCount * sizeof(uint)
    //     Triangle face indices          | 3 * trianglePatchCount * sizeof(uint)

	//   Part 4. 2-stage precomputed info:
    //     Vertex neighbors               | 3 * maxValence * sizeof(uint)
    //     Regular face offsets           | 8 * regularPatchCount * sizeof(uint)
    //     Quad face offsets              | 8 * irregularPatchCount * sizeof(uint)
    //     Triangle Face offsets          | 6 * trianglePatchCount * sizeof(uint)

	OutputHeader(outputStream, accMesh.ptr(), halfEdgeMesh.ptr(), exportText);

	if (s_outputControlPoints)
	{
		OutputControlPoints(outputStream, accMesh.ptr(), halfEdgeMesh.ptr(), reorder, exportText);
	}

	if (s_outputStencils)
	{
		OutputStencils(outputStream, accMesh.ptr(), halfEdgeMesh.ptr(), maxVertexCount, exportText);
	}

	if (s_outputMesh)
	{
		OutputMesh(outputStream, accMesh.ptr(), halfEdgeMesh.ptr(), maxValence, exportText);
	}

	if (s_outputTwoStageData)
	{
		OutputTwoStageData(outputStream, accMesh.ptr(), halfEdgeMesh.ptr(), maxValence, exportText);
	}

    return EXIT_SUCCESS;
}
