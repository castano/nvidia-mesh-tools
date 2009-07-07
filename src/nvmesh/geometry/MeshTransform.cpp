// This code is in the public domain -- castanyo@yahoo.es

#include <nvmath/Vector.h>
#include <nvmath/Matrix.h>
#include <nvmesh/BaseMesh.h>
#include <nvmesh/geometry/MeshTransform.h>
#include <nvmesh/geometry/Bounds.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>

using namespace nv;

void nv::MeshTransform::transform(HalfEdge::Mesh * mesh, const Matrix & matrix)
{
	nvCheck(mesh != NULL);

	Matrix normalXForm = transpose(inverse(matrix));

	const uint vertexCount = mesh->vertexCount();
	for (uint v = 0; v < vertexCount; v++)
	{
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);
		
		vertex->setPos(transformPoint(matrix, vertex->pos()));
		vertex->setNor(transformVector(normalXForm, vertex->nor()));
	}
}

float nv::MeshTransform::fitBox(HalfEdge::Mesh * mesh, const Box & box)
{
	nvCheck(mesh != NULL);

	Vector3 boxCenter = box.center();
	Vector3 boxExtents = box.extents();

	Box bounds = nv::MeshBounds::box(mesh);
	Vector3 boundsCenter = bounds.center();
	Vector3 boundsExtents = bounds.extents();

	Matrix matrix(identity);
	matrix.translate(boxCenter);

	Vector3 vecScale(boxExtents.x() / boundsExtents.x(), boxExtents.y() / boundsExtents.y(), boxExtents.z() / boundsExtents.z());
	float scale = min(min(vecScale.x(), vecScale.y()), vecScale.z());
	matrix.scale(scale);

	matrix.translate(-boundsCenter);

	transform(mesh, matrix);

	return scale;
}

void nv::MeshTransform::translate(HalfEdge::Mesh * mesh, const Vector3 & v)
{
	nvCheck(mesh != NULL);

	Matrix matrix(identity);
	matrix.translate(v);

	transform(mesh, matrix);
}

/*
void nv::MeshTransform::flipAxis(HalfEdge::Mesh * mesh, int a, int b, int c)
{
	nvCheck(mesh != NULL);

	// TBD
}
*/


void nv::MeshTransform::transform(BaseMesh * mesh, const Matrix & matrix)
{
	nvCheck(mesh != NULL);

	Matrix normalXForm = transpose(inverse(matrix));

	const uint vertexCount = mesh->vertexCount();
	for(uint v = 0; v < vertexCount; v++)
	{
		BaseMesh::Vertex & vertex = mesh->vertexAt(v);
		
		vertex.pos = transformPoint(matrix, vertex.pos);
		
		if (!equal(vertex.nor, Vector3(zero)))
		{
			vertex.nor = normalize(transformVector(normalXForm, vertex.nor));
		}
	}
}

float nv::MeshTransform::fitBox(BaseMesh * mesh, const Box & box)
{
	nvCheck(mesh != NULL);

	Vector3 boxCenter = box.center();
	Vector3 boxExtents = box.extents();

	Box bounds = MeshBounds::box(mesh);
	Vector3 boundsCenter = bounds.center();
	Vector3 boundsExtents = bounds.extents();

	Matrix matrix(identity);
	matrix.translate(boxCenter);

	Vector3 vecScale(boxExtents.x() / boundsExtents.x(), boxExtents.y() / boundsExtents.y(), boxExtents.z() / boundsExtents.z());
	float scale = min(min(vecScale.x(), vecScale.y()), vecScale.z());
	matrix.scale(scale);

	matrix.translate(-boundsCenter);

	transform(mesh, matrix);

	return scale;
}

void nv::MeshTransform::translate(BaseMesh * mesh, const Vector3 & v)
{
	nvCheck(mesh != NULL);

	Matrix matrix(identity);
	matrix.translate(v);

	transform(mesh, matrix);
}
