// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmath/Basis.h>
#include <nvmesh/TriMesh.h>
#include <nvmesh/mender/NVMeshMender.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/geometry/TangentSpace.h>
#include <nvmesh/geometry/MeshNormals.h>

#include "TriBaseSurface.h"

using namespace nv;

#define DELTA   (0.00001f)

static void UE3_CreateOrthonormalBasis(Vector3 & XAxis, Vector3 & YAxis, Vector3 & ZAxis)
{
	// Project the X and Y axes onto the plane perpendicular to the Z axis.
	XAxis -= dot(XAxis, ZAxis) / dot(ZAxis, ZAxis) * ZAxis;
	YAxis -= dot(YAxis, ZAxis) / dot(ZAxis, ZAxis) * ZAxis;

	// If the X axis was parallel to the Z axis, choose a vector which is orthogonal to the Y and Z axes.
	if (length_squared(XAxis) < DELTA*DELTA)
	{
		XAxis = cross(YAxis, ZAxis);
	}

	// If the Y axis was parallel to the Z axis, choose a vector which is orthogonal to the X and Z axes.
	if(length_squared(YAxis) < DELTA*DELTA)
	{
		YAxis = cross(XAxis, ZAxis);
	}

	// Normalize the basis vectors.
	XAxis = normalize(XAxis);
	YAxis = normalize(YAxis);
	ZAxis = normalize(ZAxis);
}

static void UE3_CreateOrthonormalBasis(Basis & basis)
{
	UE3_CreateOrthonormalBasis(basis.tangent, basis.bitangent, basis.normal);
}

static void UE3_QuantizeVector(Vector3 & v)
{
	v *= 127.5f;
	v += Vector3(127.5f, 127.5f, 127.5f);

	float x = float(clamp(int(v.x()), 0, 255));
	float y = float(clamp(int(v.y()), 0, 255));
	float z = float(clamp(int(v.z()), 0, 255));

	v.set(x, y, z);
	
	v -= Vector3(127.5f, 127.5f, 127.5f);
	v /= 127.5f;
}

static void UE3_QuantizeBasis(Basis & basis)
{
	float handness = basis.handness();

	UE3_QuantizeVector(basis.normal);
	UE3_QuantizeVector(basis.tangent);

	// Reconstruct bitangent as done in the vertex shader:
	basis.bitangent = cross(basis.normal, basis.tangent) * handness;

	// @@ Does the vertex shader normalize normal, tangent, and bitangent?

}




TriBaseSurface::TriBaseSurface(TriMesh * mesh) : m_mesh(mesh), m_faceIndex(0)
{
}

TriBaseSurface::~TriBaseSurface()
{
}

void TriBaseSurface::prepare(TangentSpaceMode tangentSpaceMode)
{
	if (tangentSpaceMode == TangentSpaceMode_MeshMender)
	{
		std::vector<MeshMender::Vertex> vertices;
		
		const uint vertexCount = m_mesh->vertexCount();
		vertices.resize(vertexCount);
		
		for (uint i = 0; i < vertexCount; i++)
		{
			const TriMesh::Vertex & vertex = m_mesh->vertexAt(i);

			vertices[i].pos = vertex.pos;
			vertices[i].normal = vertex.nor;
			vertices[i].s = vertex.tex.x();
			vertices[i].t = vertex.tex.y();
			vertices[i].tangent = Vector3(zero);
			vertices[i].binormal = Vector3(zero);
		}

		std::vector<unsigned int> indices;
		
		const uint faceCount = m_mesh->faceCount();
		indices.resize(3 * faceCount);

		for (uint i = 0; i < faceCount; i++)
		{
			const TriMesh::Face & face = m_mesh->faceAt(i);

			indices[3 * i + 0] = face.v[0];
			indices[3 * i + 1] = face.v[1];
			indices[3 * i + 2] = face.v[2];
		}

		std::vector<unsigned int> xrefs;

		// Default options.
		float minNormalsCreaseCosAngle = 0.0f;
		float minTangentsCreaseCosAngle = 0.0f;
		float minBinormalsCreaseCosAngle = 0.0f;
		float weightNormalsByArea = 1.0f;

		MeshMender::NormalCalcOption computeNormals = MeshMender::CALCULATE_NORMALS;
		MeshMender::ExistingSplitOption respectExistingSplits = MeshMender::DONT_RESPECT_SPLITS;
		MeshMender::CylindricalFixOption fixCylindricalWrapping = MeshMender::DONT_FIX_CYLINDRICAL;

		MeshMender mender;
		mender.Mend(vertices, indices, xrefs,
			minNormalsCreaseCosAngle,
			minTangentsCreaseCosAngle,
			minBinormalsCreaseCosAngle,
			weightNormalsByArea,
			computeNormals,
			respectExistingSplits,
			fixCylindricalWrapping);

		// Create new triMesh.
		const uint newFaceCount = indices.size() / 3;
		const uint newVertexCount = vertices.size();

		TriMesh * triMesh = new TriMesh(newFaceCount, newVertexCount);
		triMesh->faces().resize(newFaceCount);
		triMesh->vertices().resize(newVertexCount);
		m_basisArray.resize(newVertexCount);

		for (uint i = 0; i < newVertexCount; i++)
		{
			TriMesh::Vertex & vertex = triMesh->vertexAt(i);

			vertex.pos = vertices[i].pos;
			vertex.nor = vertices[i].normal;
			vertex.tex = Vector2(vertices[i].s, vertices[i].t);

			m_basisArray[i].normal = vertices[i].normal;
			m_basisArray[i].tangent = vertices[i].tangent;
			m_basisArray[i].bitangent = vertices[i].binormal;
		}

		for (uint i = 0; i < newFaceCount; i++)
		{
			TriMesh::Face & face = triMesh->faceAt(i);

			face.v[0] = indices[3 * i + 0];
			face.v[1] = indices[3 * i + 1];
			face.v[2] = indices[3 * i + 2];
		}
	}
	else if (tangentSpaceMode == TangentSpaceMode_Lengyel)
	{
		if (!MeshNormals::hasNormals(m_mesh.ptr()))
		{
			MeshNormals::computeNormals(m_mesh.ptr(), WeightFaceArea);
		}

		geometry::computeMeshTangents(m_mesh.ptr(), m_basisArray);

		foreach(i, m_basisArray)
		{
			Basis & basis = m_basisArray[i];
		//	basis.orthonormalize();
			UE3_CreateOrthonormalBasis(basis);
		}
	}
	else if (tangentSpaceMode == TangentSpaceMode_Castano)
	{
		if (!MeshNormals::hasNormals(m_mesh.ptr()))
		{
			MeshNormals::computeNormals(m_mesh.ptr(), WeightFaceArea);
		}

		geometry::computeMeshTangents(m_mesh.ptr(), m_basisArray);

		foreach(i, m_basisArray)
		{
			Basis & basis = m_basisArray[i];
			basis.robustOrthonormalize();
		}
	}
	else
	{
		nvCheck(tangentSpaceMode == TangentSpaceMode_UE3);

		if (!MeshNormals::hasNormals(m_mesh.ptr()))
		{
			MeshNormals::computeNormals(m_mesh.ptr(), WeightFaceArea);
		}

		// Compute triangle basis.
		const uint faceCount = m_mesh->faceCount();

		Array<Basis> triangleBasis;
		triangleBasis.resize(faceCount);

		for (uint i = 0; i < faceCount; i++)
		{
			const TriMesh::Face & face = m_mesh->faceAt(i);
			const TriMesh::Vertex & vertex1 = m_mesh->vertexAt(face.v[0]);
			const TriMesh::Vertex & vertex2 = m_mesh->vertexAt(face.v[1]);
			const TriMesh::Vertex & vertex3 = m_mesh->vertexAt(face.v[2]);

			Vector3 p1 = vertex1.pos;
			Vector3 p2 = vertex2.pos;
			Vector3 p3 = vertex3.pos;

			Matrix parameterToLocal(
				Vector4(p2 - p1,  0),
				Vector4(p3 - p1,  0),
				Vector4(p1,       0),
				Vector4(0, 0, 0,  1));

			Vector2 t1 = vertex1.tex;
			Vector2 t2 = vertex2.tex;
			Vector2 t3 = vertex3.tex;

			Matrix parameterToTexture(
				Vector4(t2 - t1, 0, 0),
				Vector4(t3 - t1, 0, 0),
				Vector4(t1,      1, 0),
				Vector4(0, 0,    0, 1));

			Vector3 tangentX(zero);
			Vector3 tangentY(zero);
			Vector3 normal = normalizeSafe(cross(p2 - p1, p3 - p1), Vector3(zero), 0.0f);

			if (parameterToTexture.determinant() != 0.0f)
			{
				const Matrix textureToLocal = mul(parameterToLocal, inverse(parameterToTexture));

				tangentX = normalizeSafe(transformVector(textureToLocal, Vector3(1,0,0)), Vector3(zero), 0.0f);
				tangentY = normalizeSafe(transformVector(textureToLocal, Vector3(0,1,0)), Vector3(zero), 0.0f);
			}

			tangentX -= normal * dot(tangentX, normal);
			tangentY -= normal * dot(tangentY, normal);

			triangleBasis[i].tangent = normalizeSafe(tangentX, Vector3(zero), 0.0f);
			triangleBasis[i].bitangent = normalizeSafe(tangentY, Vector3(zero), 0.0f);
			triangleBasis[i].normal = normal;
		}

		const uint vertexCount = m_mesh->vertexCount();
		m_basisArray.resize(vertexCount);

#if 0
		for (uint v = 0; v < vertexCount; v++)
		{
			m_basisArray[v].tangent = Vector3(zero);
			m_basisArray[v].bitangent = Vector3(zero);
			m_basisArray[v].normal = Vector3(zero);
		}


		// @@ Average triangle basis?

	    for (uint f = 0; f < faceCount; f++)
		{
			const TriMesh::Face & face = mesh->faceAt(f);

			for (int i = 0; i < 3; i++)
			{
				const uint idx = face.v[i];

				m_basisArray[idx].tangent += triangleBasis[f].tangent;
				m_basisArray[idx].bitangent += triangleBasis[f].bitangent;
				m_basisArray[idx].normal += triangleBasis[f].normal;
			}
		}
#else
		for (uint f = 0; f < faceCount; f++)
		{
			const TriMesh::Face & face = m_mesh->faceAt(f);

			Vector3 vertexTangentX[3];
			Vector3 vertexTangentY[3];
			Vector3 vertexTangentZ[3];

			for (uint i = 0; i < 3; i++)
			{
				vertexTangentX[i] = Vector3(zero);
				vertexTangentY[i] = Vector3(zero);
				vertexTangentZ[i] = Vector3(zero);
			}

			Vector3 triangleNormal = triangleBasis[f].normal;
			float determinant = triangleBasis[f].determinant();

			for (uint g = 0; g < faceCount; g++)
			{
				const TriMesh::Face & otherFace = m_mesh->faceAt(g);

				Vector3 otherTriangleNormal = triangleBasis[g].normal;
				float otherDeterminant = triangleBasis[g].determinant();

				for (uint i = 0; i < 3; i++)
				{
					for(uint j = 0; j < 3; j++)
					{
						if (equal(m_mesh->vertexAt(face.v[i]).pos, m_mesh->vertexAt(otherFace.v[j]).pos, 0.0f))
						{
							if (determinant * otherDeterminant > 0.0f && 
								equal(m_mesh->vertexAt(face.v[i]).tex, m_mesh->vertexAt(otherFace.v[j]).tex, 0.0f))
							{
								vertexTangentX[i] += triangleBasis[g].tangent;
								vertexTangentY[i] += triangleBasis[g].bitangent;
							}

							// @@ IC: Ignore "smoothing" edges. Always add normals together.
							// Only contribute 'normal' if the vertices are truly one and the same to obey hard "smoothing" edges baked into 
							// the mesh by vertex duplication.. - Erik
							{
								vertexTangentZ[i] += otherTriangleNormal;
							}
						}
					}
				}
			}

			// ...

			for (uint i = 0; i < 3; i++)
			{
				Vector3 tangentX = normalize(vertexTangentX[i], 0.0f);
				Vector3 tangentY = normalize(vertexTangentY[i], 0.0f);
				Vector3 tangentZ = normalize(vertexTangentZ[i], 0.0f);

				tangentY -= tangentX * dot(tangentX, tangentY);
				tangentY = normalize(tangentY);

				tangentX -= tangentZ * dot(tangentZ, tangentX);
				tangentY -= tangentZ * dot(tangentZ, tangentY);

				tangentX = normalize(tangentX);
				tangentY = normalize(tangentY);

				m_basisArray[face.v[i]].tangent = tangentX;
				m_basisArray[face.v[i]].bitangent = -tangentY;	// @@ For some strange reason UE3 negates the bitangent!
				m_basisArray[face.v[i]].normal = tangentZ;
			}
		}
#endif

	}
}

void TriBaseSurface::setInterpolationMode(PositionMode positionMode, NormalMode normalMode)
{
	m_positionMode = positionMode;
	m_normalMode = normalMode;
}


uint TriBaseSurface::faceCount() const
{
	nvCheck(m_mesh != NULL);
	return m_mesh->faceCount();
}

void TriBaseSurface::selectFace(uint i)
{
	nvDebugCheck(m_mesh != NULL);
	m_faceIndex = i;

	// @@ Compute PN Triangles and MN Triangles control points.

}

FaceDomain::Enum TriBaseSurface::domain() const
{
	return FaceDomain::Triangle;
}

void TriBaseSurface::textureCoordinates(Vector2 * texCoordArray) const
{
	nvDebugCheck(m_mesh != NULL);

	const TriMesh::Face & face = m_mesh->faceAt(m_faceIndex);

	texCoordArray[0] = m_mesh->vertexAt(face.v[0]).tex;
	texCoordArray[1] = m_mesh->vertexAt(face.v[1]).tex;
	texCoordArray[2] = m_mesh->vertexAt(face.v[2]).tex;
}


void TriBaseSurface::evaluate(float u, float v, Vector3 * pos, Basis * patchFrame, Basis * chartFrame) const
{
	const float w = 1 - u - v;

	const TriMesh::Face & face = m_mesh->faceAt(m_faceIndex);

	uint i0 = face.v[0];
	uint i1 = face.v[1];
	uint i2 = face.v[2];

	Vector3 pos0 = m_mesh->vertexAt(face.v[0]).pos;
	Vector3 pos1 = m_mesh->vertexAt(face.v[1]).pos;
	Vector3 pos2 = m_mesh->vertexAt(face.v[2]).pos;

	Vector3 nor0 = m_basisArray[i0].normal;
	Vector3 nor1 = m_basisArray[i1].normal;
	Vector3 nor2 = m_basisArray[i2].normal;

	//if (m_positionMode == PositionMode_Linear)
	{
		*pos = pos0 * u + pos1 * v + pos2 * w;
	}

	//if (m_normalMode == NormalMode_Linear)
	{
		patchFrame->normal = normalizeSafe(nor0 * u + nor1 * v + nor2 * w, Vector3(zero), 0.0f);
	}

	Vector3 tan0 = m_basisArray[i0].tangent;
	Vector3 tan1 = m_basisArray[i1].tangent;
	Vector3 tan2 = m_basisArray[i2].tangent;

	chartFrame->tangent = normalizeSafe(tan0 * u + tan1 * v + tan2 * w, Vector3(zero), 0.0f);

	Vector3 bit0 = m_basisArray[i0].bitangent;
	Vector3 bit1 = m_basisArray[i1].bitangent;
	Vector3 bit2 = m_basisArray[i2].bitangent;

	chartFrame->bitangent = normalizeSafe(bit0 * u + bit1 * v + bit2 * w, Vector3(zero), 0.0f);

	chartFrame->normal = normalizeSafe(cross(chartFrame->tangent, chartFrame->bitangent), Vector3(zero), 0.0f);

	// Align chart frame to surface.
	if (true)
	{
		chartFrame->normal = patchFrame->normal;

		// Project bitangent to surface normal.
		chartFrame->tangent -= patchFrame->normal * dot(patchFrame->normal, chartFrame->tangent);
		chartFrame->tangent = normalize(chartFrame->tangent);

		// Project bitangent to surface normal.
		chartFrame->bitangent -= patchFrame->normal * dot(patchFrame->normal, chartFrame->bitangent);
		chartFrame->bitangent = normalize(chartFrame->bitangent);	
	}
}

