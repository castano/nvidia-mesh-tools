// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include "AccBaseSurface.h"

#include <nvmesh/subdiv/AccMesh.h>
#include <nvmesh/subdiv/AccPatch.h>

#include <nvmath/Basis.h>

using namespace nv;


AccBaseSurface::AccBaseSurface(const AccMesh * mesh) : m_mesh(mesh), m_surfacePatch(NULL), m_texCoordPatch(NULL)
{
	nvDebugCheck(mesh != NULL);
}

AccBaseSurface::~AccBaseSurface()
{
}

uint AccBaseSurface::faceCount() const
{
    return m_mesh->regularPatchCount() + m_mesh->quadPatchCount() + m_mesh->triPatchCount();
}

void AccBaseSurface::selectFace(uint i)
{
    AccMesh::Patch patch = m_mesh->patchAt(i);

    if (patch.bezierAccPatch != NULL) m_surfacePatch = patch.bezierAccPatch;
    else if (patch.gregoryAccPatch != NULL) m_surfacePatch = patch.gregoryAccPatch;
    //else if (patch.pmAccPatch != NULL) m_surfacePatch = patch.pmAccPatch;

	m_texCoordPatch = patch.texCoordPatch;
	m_tangentSpacePatch = patch.tangentSpacePatch;
}

FaceDomain::Enum AccBaseSurface::domain() const
{
	if (m_surfacePatch->isTriangle())
	{
		return FaceDomain::Triangle;
	}
	
	return FaceDomain::Quad;
}

void AccBaseSurface::textureCoordinates(Vector2 * texCoordArray) const
{
	nvDebugCheck(texCoordArray);
	
	const uint vertexCount = m_texCoordPatch->isTriangle() ? 3 : 4;

	for (uint i = 0; i < vertexCount; i++)
	{
		// Copy interior texture coordinates.
		texCoordArray[i] = m_texCoordPatch->texCoord(i, 0);
	}
}

void AccBaseSurface::evaluate(float u, float v, Vector3 * pos, Basis * patchFrame, Basis * chartFrame) const
{
	u = clamp(u, 0.0f, 1.0f);
	v = clamp(v, 0.0f, 1.0f);

	//m_patch->evaluatePosition(u, v, pos);
	//m_patch->evaluatePatchFrame(u, v, patchFrame);

	m_surfacePatch->evaluateSurface(u, v, pos, &patchFrame->tangent, &patchFrame->bitangent);
	patchFrame->normal = normalizeSafe(cross(patchFrame->tangent, patchFrame->bitangent), Vector3(zero), 0);

    // @@ Enable tangent space patch generation!
	/*m_tangentSpacePatch->evaluateChartFrame(u, v, chartFrame);

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
	}*/

	/*
	if (!basis.isValid())
	{
		bool b = basis.isValid();
		
		nvDebug("[%f] %d %d %d\n", basis.determinant(), b, x, y);
		nvDebug("%f %f %f\n", basis.normal.x(), basis.normal.y(), basis.normal.z());
		nvDebug("%f %f %f\n", basis.tangent.x(), basis.tangent.y(), basis.tangent.z());
		nvDebug("%f %f %f\n", basis.bitangent.x(), basis.bitangent.y(), basis.bitangent.z());
	}
	*/
}


