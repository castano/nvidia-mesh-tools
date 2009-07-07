// This code is in the public domain -- castanyo@yahoo.es

#include <nvmesh/TriMesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/param/Atlas.h>
#include <nvmesh/param/SingleFaceMap.h>
#include <nvmesh/param/LeastSquaresConformalMap.h>
//#include <nvmesh/geometry/MeshNormals.h>
//#include <nvmesh/geometry/TangentSpace.h>

#include "MeshNormals.h"
#include "TangentSpace.h"


using namespace nv;


void nv::geometry::computeMeshTangents(const TriMesh * mesh, Array<Basis> & meshBasis)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	meshBasis.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		meshBasis[v].tangent = Vector3(zero);
		meshBasis[v].bitangent = Vector3(zero);
	}

	const uint faceCount = mesh->faceCount();
    for (uint f = 0; f < faceCount; f++)
    {
		const TriMesh::Face & face = mesh->faceAt(f);

        const uint i1 = face.v[0];
        const uint i2 = face.v[1];
        const uint i3 = face.v[2];
        
        const Vector3 & v1 = mesh->vertexAt(i1).pos;
        const Vector3 & v2 = mesh->vertexAt(i2).pos;
        const Vector3 & v3 = mesh->vertexAt(i3).pos;
        
        const Vector2 & w1 = mesh->vertexAt(i1).tex;
        const Vector2 & w2 = mesh->vertexAt(i2).tex;
        const Vector2 & w3 = mesh->vertexAt(i3).tex;
        
		float x1 = v2.x() - v1.x();
		float x2 = v3.x() - v1.x();
		float y1 = v2.y() - v1.y();
		float y2 = v3.y() - v1.y();
		float z1 = v2.z() - v1.z();
		float z2 = v3.z() - v1.z();
	    
		float s1 = w2.x() - w1.x();
		float s2 = w3.x() - w1.x();
		float t1 = w2.y() - w1.y();
		float t2 = w3.y() - w1.y();

		const float denom = s1 * t2 - s2 * t1;

		if (denom != 0.0f)
		{
			const float r = 1.0f / denom;

			Vector3 sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r, (t2 * z1 - t1 * z2) * r);
			Vector3 tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r, (s1 * z2 - s2 * z1) * r);

			meshBasis[i1].tangent += sdir;
			meshBasis[i1].bitangent += tdir;

			meshBasis[i2].tangent += sdir;
			meshBasis[i2].bitangent += tdir;

			meshBasis[i3].tangent += sdir;
			meshBasis[i3].bitangent += tdir;
		}
    }

	meshBasis.resize(vertexCount);

	// Normalize, but do not orthonormalize.
	for (uint v = 0; v < vertexCount; v++)
	{
		const TriMesh::Vertex & vertex = mesh->vertexAt(v);

		meshBasis[v].normal = normalize(vertex.nor, 0.0f);
		meshBasis[v].tangent = normalize(meshBasis[v].tangent, 0.0f);
		meshBasis[v].bitangent = normalize(meshBasis[v].bitangent, 0.0f);
	}
}


struct VertexBasis
{
	VertexBasis() : 
		uu_sum(0),
		vv_sum(0),
		uv_sum(0),
		qu_sum(zero),
		qv_sum(zero) {}

	float uu_sum;
	float vv_sum;
	float uv_sum;
	Vector3 qu_sum;
	Vector3 qv_sum;

	void evaluateTangents(Vector3 & tan0, Vector3 & tan1) const
	{
		const float denom = uu_sum * vv_sum - uv_sum * uv_sum;
		nvCheck(!equal(denom, 0.0f, 0.0f));

		const float factor = 1.0f / denom;

		tan0 = (qu_sum * vv_sum - qv_sum * uv_sum) * factor;
		tan1 = (qv_sum * uu_sum - qu_sum * uv_sum) * factor;
	}
};

void nv::geometry::computeLeastSquaresMeshTangents(const TriMesh * mesh, Array<Basis> & meshBasis)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();

	Array<VertexBasis> vertexBasisArray;
	vertexBasisArray.resize(vertexCount);

	const uint faceCount = mesh->faceCount();
    for (uint f = 0; f < faceCount; f++)
    {
		const TriMesh::Face & face = mesh->faceAt(f);

		// @@ non-boundary and non-seam edges are being added twice!

		for (uint i = 0, j = 3; i < 3; j = i, i++)
		{
			const TriMesh::Vertex & from = mesh->vertexAt(face.v[j]);
			const TriMesh::Vertex & to = mesh->vertexAt(face.v[i]);

			const Vector3 q = to.pos - from.pos;
			const Vector2 v = to.tex - from.tex;

			vertexBasisArray[face.v[i]].qu_sum += q * v.x();
			vertexBasisArray[face.v[i]].qv_sum += q * v.y();

			vertexBasisArray[face.v[i]].uu_sum += v.x() * v.x();
			vertexBasisArray[face.v[i]].uv_sum += v.x() * v.y();
			vertexBasisArray[face.v[i]].vv_sum += v.y() * v.y();

			vertexBasisArray[face.v[j]].qu_sum += q * v.x();
			vertexBasisArray[face.v[j]].qv_sum += q * v.y();

			vertexBasisArray[face.v[j]].uu_sum += v.x() * v.x();
			vertexBasisArray[face.v[j]].uv_sum += v.x() * v.y();
			vertexBasisArray[face.v[j]].vv_sum += v.y() * v.y();
		}
	}

	for (uint v = 0; v < vertexCount; v++)
	{
		const TriMesh::Vertex & vertex = mesh->vertexAt(v);

		vertexBasisArray[v].evaluateTangents(meshBasis[v].tangent, meshBasis[v].bitangent);

		// Normalize basis vectors:
		meshBasis[v].normal = normalize(vertex.nor, 0.0f);
		meshBasis[v].tangent = normalize(meshBasis[v].tangent, 0.0f);
		meshBasis[v].bitangent = normalize(meshBasis[v].bitangent, 0.0f);
	}
}




// This is Lengyel's algorithm adapted to polygonal meshes. It considers the triangles formed by 
// each pair of edges.
//
void nv::geometry::computeLengyelTangentBasis(const HalfEdge::Vertex * vertex, Basis * basis)
{
	Vector3 tan0(zero);
	Vector3 tan1(zero);

	// traverse edges around vertex.
	for(HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge0 = it.current();
		nvDebugCheck(edge0 != NULL);

		if (edge0->face() == NULL) continue; // skip boundary edges.
		
		nvDebugCheck(edge0->from()->pos() == vertex->pos());
		//nvDebugCheck(edge0->from() == vertex);
		if (edge0->from()->tex() != vertex->tex()) continue;
		//if (edge0->from() != vertex) continue;

		const HalfEdge::Edge * edge1 = edge0->prev();
		nvDebugCheck(edge1 != NULL);
		//nvDebugCheck(edge1->to() == vertex);
		nvDebugCheck(edge1->to()->tex() == vertex->tex());
		//if (edge1->to() != vertex) continue;

		Vector3 v1 = vertex->pos();
		Vector3 v2 = edge0->to()->pos();
		Vector3 v3 = edge1->from()->pos();

		Vector2 w1 = vertex->tex();
		Vector2 w2 = edge0->to()->tex();
		Vector2 w3 = edge1->from()->tex();

		float x1 = v2.x() - v1.x();
		float x2 = v3.x() - v1.x();
		float y1 = v2.y() - v1.y();
		float y2 = v3.y() - v1.y();
		float z1 = v2.z() - v1.z();
		float z2 = v3.z() - v1.z();
	    
		float s1 = w2.x() - w1.x();
		float s2 = w3.x() - w1.x();
		float t1 = w2.y() - w1.y();
		float t2 = w3.y() - w1.y();

		float r = 1.0F / (s1 * t2 - s2 * t1);

		nvCheck(!equal(s1 * t2 - s2 * t1, 0.0f, 0.0f));

		//if (!equal(s1 * t2 - s2 * t1, 0.0f, 0.0f))
		{
			Vector3 sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r, (t2 * z1 - t1 * z2) * r);
			Vector3 tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r, (s1 * z2 - s2 * z1) * r);

			tan0 += sdir;
			tan1 += tdir;
		}
	}

	basis->normal = vertex->nor();
	basis->tangent = tan0;
	basis->bitangent = tan1;
}


// Solve this equation system in the least squares sense:
//
// [Q1 .. QN] = [V1 .. VN] [B1 B2]
//
// (VtV)^-1 VtQ = B
//
// Ax = b -> (VtV) B = VtQ
//
void nv::geometry::computeLeastSquaresTangentBasis(const HalfEdge::Vertex * vertex, Basis * basis)
{
	nvCheck(vertex != NULL);
	nvCheck(basis != NULL);

	float uu_sum = 0;
	float vv_sum = 0;
	float uv_sum = 0;

	Vector3 qu_sum(zero);
	Vector3 qv_sum(zero);

	int i = 0;
	for (HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), i++)
	{
		const HalfEdge::Edge * edge = it.current();
		nvDebugCheck(edge != NULL);

		nvDebugCheck(edge->from()->pos() == vertex->pos());

		const HalfEdge::Vertex * from = vertex;
		const HalfEdge::Vertex * to = edge->to();

		// Skip edges that are not part of this vertex's chart.
		if (edge->face() == NULL || edge->from()->tex() != vertex->tex())
		{
			const HalfEdge::Edge * pair = edge->pair();
			nvCheck(pair->to()->pos() == vertex->pos());

			if (pair->face() == NULL || pair->to()->tex() != vertex->tex()) {
				continue;
			}

			to = pair->from();
		}

		const Vector3 q = to->pos() - from->pos();
		const Vector2 v = to->tex() - from->tex();

		// Setup edge equations:
		qu_sum += q * v.x();
		qv_sum += q * v.y();

		uu_sum += v.x() * v.x();
		uv_sum += v.x() * v.y();
		vv_sum += v.y() * v.y();
	}

	// Solve equation system:
	float denom = uu_sum * vv_sum - uv_sum * uv_sum;
	nvCheck(!equal(denom, 0.0f, 0.0f));

	float factor = 1.0f / denom;

	Vector3 tan0 = (qu_sum * vv_sum - qv_sum * uv_sum) * factor;
	Vector3 tan1 = (qv_sum * uu_sum - qu_sum * uv_sum) * factor;

	basis->normal = vertex->nor();
	basis->tangent = tan0;
	basis->bitangent = tan1;
}




void nv::geometry::computeMeshTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	meshBasis.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		Basis basis;
		//computeLengyelTangentBasis(vertex, &basis);
		computeLeastSquaresTangentBasis(vertex, &basis);

		basis.normal = MeshNormals::computeCatmullClarkNormal(vertex);	// @@ !!!

		//basis.robustOrthonormalize();
		basis.tangent -= basis.normal * dot(basis.normal, basis.tangent);
		basis.tangent = normalize(basis.tangent);
		basis.bitangent -= basis.normal * dot(basis.normal, basis.bitangent);
		basis.bitangent = normalize(basis.bitangent);

		meshBasis[v] = basis;
	}
}


static void computeCatmullClarkBasis(const HalfEdge::Vertex * vertex, Basis * basis)
{
	Vector3 normal = MeshNormals::computeCatmullClarkNormal(vertex);

	Vector3 tan0(zero);
	Vector3 tan1(zero);

	uint colocalCount = vertex->colocalCount();
	uint valence = vertex->valence();

	// traverse edges around vertex.
	for(HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge0 = it.current();
		nvDebugCheck(edge0 != NULL);

		if (edge0->face() == NULL) continue; // skip boundary edges.
		
		nvDebugCheck(edge0->from()->pos() == vertex->pos());
		if (edge0->from()->tex() != vertex->tex()) continue;
	
		const HalfEdge::Edge * edge1 = edge0->prev();
		nvDebugCheck(edge1 != NULL);
		nvDebugCheck(edge1->to()->tex() == vertex->tex());

		Vector3 v1 = vertex->pos();
		Vector3 v2 = edge0->to()->pos();
		Vector3 v3 = edge1->from()->pos();

		Vector2 w1 = vertex->tex();
		Vector2 w2 = edge0->to()->tex();
		Vector2 w3 = edge1->from()->tex();

		Vector3 e1 = MeshNormals::computeCatmullClarkTangent(edge0);
		Vector3 e2 = MeshNormals::computeCatmullClarkTangent(edge1->pair());

		if (equal(e1, -e2))
		{
			// This happens in corners, but also on some interior vertices...
			e1 = normalize(v2 - v1, 0);
			e2 = normalize(v3 - v1, 0);
		}
		else
		{
			// @@ This should be true in the case above!! OK for corners, but not for interior!!
			nvCheck(dot(cross(e1, e2), normal) > 0.0f);
		}

		float x1 = e1.x();
		float x2 = e2.x();
		float y1 = e1.y();
		float y2 = e2.y();
		float z1 = e1.z();
		float z2 = e2.z();
		
		Vector2 st1 = w2 - w1;
		Vector2 st2 = w3 - w1;

		//st1 = normalize(st1, 0);
		//st2 = normalize(st2, 0);

		float s1 = st1.x();
		float s2 = st2.x();
		float t1 = st1.y();
		float t2 = st2.y();

		float r = 1.0F / (s1 * t2 - s2 * t1);

		nvCheck(!equal(s1 * t2 - s2 * t1, 0.0f, 0.0f));

		if (!equal(s1 * t2 - s2 * t1, 0.0f, 0.0f))
		{
			Vector3 sdir((t2 * x1 - t1 * x2) * r, (t2 * y1 - t1 * y2) * r, (t2 * z1 - t1 * z2) * r);
			Vector3 tdir((s1 * x2 - s2 * x1) * r, (s1 * y2 - s2 * y1) * r, (s1 * z2 - s2 * z1) * r);

			tan0 += sdir;
			tan1 += tdir;
		}
	}

	basis->normal = normal;
	basis->tangent = tan0;
	basis->bitangent = tan1;
}

static void computeLeastSquaresCatmullClarkBasis(const HalfEdge::Vertex * vertex, Basis * basis)
{
	nvCheck(vertex != NULL);
	nvCheck(basis != NULL);

	Vector3 normal = MeshNormals::computeCatmullClarkNormal(vertex);

	float uu_sum = 0;
	float vv_sum = 0;
	float uv_sum = 0;

	Vector3 qu_sum(zero);
	Vector3 qv_sum(zero);

	int i = 0;
	for (HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), i++)
	{
		const HalfEdge::Edge * edge = it.current();
		nvDebugCheck(edge != NULL);

		nvDebugCheck(edge->from()->pos() == vertex->pos());

		const HalfEdge::Vertex * from = vertex;
		const HalfEdge::Vertex * to = edge->to();

		if (edge->face() == NULL || edge->from()->tex() != vertex->tex())
		{
			const HalfEdge::Edge * pair = edge->pair();
			nvCheck(pair->to()->pos() == vertex->pos());

			if (pair->face() == NULL || pair->to()->tex() != vertex->tex()) {
				continue;
			}

			to = pair->from();
		}

		const Vector3 ccq = MeshNormals::computeCatmullClarkTangent(edge);
		const Vector3 q = to->pos() - from->pos();
		const Vector2 v = to->tex() - from->tex();

		// Setup edge equations:
		qu_sum += ccq * v.x();
		qv_sum += ccq * v.y();

		uu_sum += v.x() * v.x();
		uv_sum += v.x() * v.y();
		vv_sum += v.y() * v.y();
	}

	// Solve equation system:
	float denom = uu_sum * vv_sum - uv_sum * uv_sum;
	if (denom == 0)
	{
		// @@ this happens on isolated vertices.
		// @@ this happens on colocal that belonged to a face that was removed (ie. when removing triangles from quad-tri mesh).

		if (length(normal) == 0.0f)
		{
			basis->tangent = Vector3(1, 0, 0);
			basis->bitangent = Vector3(0, 1, 0);
			basis->normal = Vector3(0, 0, 1);
		}
		else
		{
			basis->buildFrameForDirection(normalize(normal));
		}
	}
	else
	{
		float factor = 1.0f / denom;

		Vector3 tan0 = (qu_sum * vv_sum - qv_sum * uv_sum) * factor;
		Vector3 tan1 = (qv_sum * uu_sum - qu_sum * uv_sum) * factor;

		basis->normal = normal;
		basis->tangent = tan0;
		basis->bitangent = tan1;
	}
}


void nv::geometry::computeCatmullClarkTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	meshBasis.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		Basis basis;
		//computeCatmullClarkBasis(vertex, &basis);
		computeLeastSquaresCatmullClarkBasis(vertex, &basis);

	//	basis.robustOrthonormalize();
		basis.tangent -= basis.normal * dot(basis.normal, basis.tangent);
		basis.tangent = normalize(basis.tangent);
		basis.bitangent -= basis.normal * dot(basis.normal, basis.bitangent);
		basis.bitangent = normalize(basis.bitangent);

		meshBasis[v] = basis;
	}

#if 0
	const uint faceCount = mesh->faceCount();
	for (uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = mesh->faceAt(f);
		nvDebugCheck(face != NULL);

		const uint edgeCount = face->edgeCount();
		nvCheck(edgeCount <= 4);

		Basis basis[4];
		float det[4];

		uint b = 0;
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), b++)
		{
			const HalfEdge::Edge * edge = it.current();

			computeCatmullClarkBasis(edge->vertex(), &basis[b]);

			det[b] = basis[b].determinant();
		}

		// @@ This is is not always true:
		/*nvCheck(det[0] != 0);
		for (uint i = 1; i < edgeCount; i++)
		{
			nvCheck(det[i] != 0);
			nvCheck((det[0] < 0) == (det[i] < 0));
		}*/
	}
#endif
}


void nv::geometry::computeConformalMeshTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	meshBasis.resize(vertexCount);

	Atlas atlas(mesh);
	atlas.extractCharts();

	const uint chartCount = atlas.chartCount();
	for (uint i = 0; i < chartCount; i++)
	{
		Chart * chart = atlas.chartAt(i);
		nvDebugCheck(chart != NULL);

		chart->closeHoles();

		if (chart->isDisk())
		{
			if (chart->faceCount() == 1)
			{
				computeSingleFaceMap(chart->mesh());
			}
			else
			{
				computeLeastSquaresConformalMap(chart->mesh());
			}
		}
		else
		{
			nvDebug("### Warning: chart '%d' is not a disk.\n", i);
		}

		Array<Basis> chartBasis;
		geometry::computeMeshTangents(chart->mesh(), chartBasis);

		// Store chart basis frames in mesh basis array.
		const uint chartVertexCount = chart->vertexCount();
		for (uint i = 0; i < chartVertexCount; i++)
		{
			meshBasis[chart->vertexIndex(i)] = chartBasis[i];
		}
	}
}

void nv::geometry::computeConformalCatmullClarkTangents(const HalfEdge::Mesh * mesh, Array<Basis> & meshBasis)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	meshBasis.resize(vertexCount);

	// Backup original texture coordinates.
	Array<Vector2> texCoordArray;
	texCoordArray.resize(vertexCount);

	for (uint v = 0; v < vertexCount; v++)
	{
		texCoordArray[v] = mesh->vertexAt(v)->tex();
	}

	// We are going to modify the original mesh, but we will restore it later.
	HalfEdge::Mesh * mutableMesh = const_cast<HalfEdge::Mesh *>(mesh);


	Atlas atlas(mesh);
	atlas.extractCharts();

	const uint chartCount = atlas.chartCount();
	for (uint i = 0; i < chartCount; i++)
	{
		Chart * chart = atlas.chartAt(i);
		nvDebugCheck(chart != NULL);

		chart->closeHoles();

		if (chart->isDisk())
		{
			if (chart->faceCount() == 1)
			{
				computeSingleFaceMap(chart->mesh());
			}
			else
			{
				computeLeastSquaresConformalMap(chart->mesh());
			}
		}
		else
		{
			nvDebug("### Warning: chart '%d' is not a disk.\n", i);
		}

		// We cannot compute the tangents per patch, we need neighborhood information.
		// So, we transfer the texture coordinates to the original mesh.

		// Store chart vertices in original mesh.
		const uint chartVertexCount = chart->vertexCount();
		for (uint i = 0; i < chartVertexCount; i++)
		{
			mutableMesh->vertexAt(chart->vertexIndex(i))->setTex(chart->mesh()->vertexAt(i)->tex());
		}
	}

	// Compute tangent space basis now.
	geometry::computeCatmullClarkTangents(mesh, meshBasis);

	// Restore original texcoords.
	for (uint v = 0; v < vertexCount; v++)
	{
		mutableMesh->vertexAt(v)->setTex(texCoordArray[v]);
	}
}


