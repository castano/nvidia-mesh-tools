// Copyright NVIDIA Corporation 2006 -- Denis Kovacs <dkovacs@nvidia.com>

#include <nvcore/StdStream.h>

#include <nvmath/Box.h>

#include <nvimage/ImageIO.h>

#include <nvmesh/MeshTopology.h>
#include <nvmesh/TriMesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>
#include <nvmesh/raster/Raster.h>

#include <nvmesh/import/MeshImport.h>

#include <nvmesh/subdiv/Subdivide.h>
#include <nvmesh/subdiv/AccMeshBuilder.h>
#include <nvmesh/subdiv/AccMesh.h>
#include <nvmesh/subdiv/RemapFaces.h>

#include <nvmesh/geometry/MeshTransform.h>

#include "BaseMeshPass.h"
#include "AccBaseSurface.h"
#include "TriBaseSurface.h"
#include "Samplers.h"


using namespace nv;


BaseMeshPass::BaseMeshPass() : 
	m_mesh(NULL),
	m_geometryMap(0),
	m_geometryMask(0),
	m_meshFileName(""),
	m_subdivisionLevels(0),
	m_baseMeshType(0),
	m_displacements3D(false),
	m_tangentSpace(false)
{
}

BaseMeshPass::~BaseMeshPass()
{
}

struct BaseMeshPass::Options : CmdOptions
{
private:
	TCLAP::ValuesConstraint<int> baseMeshTypeConstr;
	const static int baseMeshTypeValues[3];
	static std::vector<int> baseMeshTypeValuesV;

public:
	TCLAP::ValueArg<std::string> meshFileNameArg;
	TCLAP::ValueArg<int>  subdivisionLevelsArg;
	TCLAP::SwitchArg displacements3DArg;
	TCLAP::SwitchArg tangentSpaceArg;
	TCLAP::ValueArg<int> baseMeshTypeArg;

	// The constructor initializes the options
	Options() :
		baseMeshTypeConstr(baseMeshTypeValuesV),
		meshFileNameArg("", "lomesh", "Low resolution mesh.", true, /*default*/ "", "filename"),
		subdivisionLevelsArg("", "losub", "Low resolution subdivision levels (defaults to 0).", false, 0, "levels"),
		displacements3DArg("", "disp3D", "Use 3D displacements.", /*default*/ false), 
		tangentSpaceArg("", "tangentSpace", "Use tangent space.", /*default*/ false),
		baseMeshTypeArg("", "lomode", "Base mesh mode:\n\t0\tBezier ACC Surface (default)\n\t1\tGregory ACC Surface\n\t2\tTriangular Mesh.", false, /*defautl*/ 0, &baseMeshTypeConstr) {}

	void add(TCLAP::CmdLine &cmd) {
		cmd.add(baseMeshTypeArg);
		cmd.add(tangentSpaceArg);
		cmd.add(displacements3DArg);
		cmd.add(subdivisionLevelsArg);
		cmd.add(meshFileNameArg);
	}
};

// @@ This is a bit of a hack.
const int BaseMeshPass::Options::baseMeshTypeValues[3] = {0,1,2};
std::vector<int> BaseMeshPass::Options::baseMeshTypeValuesV(baseMeshTypeValues, baseMeshTypeValues+3);

CmdOptions * BaseMeshPass::getCmdOptions() {
	return new Options();
}

void BaseMeshPass::setCmdOptions(CmdOptions * opt) {
	nvDebugCheck(opt != NULL);
	Options & options = * dynamic_cast<Options*>(opt);

	nvDebugCheck( options.meshFileNameArg.isSet() );
	m_meshFileName = options.meshFileNameArg.getValue();

	m_subdivisionLevels = options.subdivisionLevelsArg.getValue();
	m_displacements3D   = options.displacements3DArg.getValue();
	m_tangentSpace      = options.tangentSpaceArg.getValue();
	m_baseMeshType      = options.baseMeshTypeArg.getValue();
}


bool BaseMeshPass::loadMesh()
{
	freeMesh();

	if (m_meshFileName.empty()) {
		printf("No base mesh given. Skipping displacement map.\n");
		return false;
	}

	nvDebug("Importing '%s'\n", m_meshFileName.c_str());
	AutoPtr<MeshImport> importer(MeshImport::importer(m_meshFileName.c_str()));

	if (importer == NULL) {
		nvDebug("Error, unkown file type '%s'\n", m_meshFileName.c_str());
		return false;
	}

	StdInputStream stream(m_meshFileName.c_str());
	if (stream.isError()) {
		nvDebug("Error, cannot open '%s'\n", m_meshFileName.c_str());
		return false;
	}

	importer->import(&stream);

	if (m_baseMeshType == 0 || m_baseMeshType == 1)
	{
		AutoPtr<HalfEdge::Mesh> mesh( importer->builder().buildHalfEdgeMesh() );

		if (mesh == NULL) {
			nvDebug("Failure to build HalfEdge Mesh\n");
			return false;
		}

        // @@ WTF! Do not scale the input mesh, and if we do we must scale the detailed mesh too.
		//MeshTransform::fitBox(mesh.ptr(), Box(Vector3(-1,-1,-1), Vector3(1,1,1)));

		for (int i = 0; i < m_subdivisionLevels; i++)
		{
			mesh = Subdivide::doCatmullClarkSplit(mesh.ptr());
		}

		// Make sure the orientation of the patches 
		// is the same as in the render.
		RemapFaces::minimizeTopologyCount(mesh.ptr());

        uint buildFlags = BuildFlags_GenerateAll | BuildFlags_BuildTexCoordPatches;
		if (m_baseMeshType == 1)
		{
            buildFlags |= BuildFlags_BuildGregoryPatches;
		}

		AccMeshBuilder builder(mesh.ptr(), BoundaryMode_Spline);
		m_mesh = new AccBaseSurface(builder.buildAccMesh(buildFlags));
	}
	else if (m_baseMeshType == 2)
	{
		TriMesh * mesh = importer->builder().buildTriMesh();

		//MeshTransform::fitBox(mesh, Box(Vector3(-1,-1,-1), Vector3(1,1,1)));
		MeshTransform::transform(mesh, 
			Matrix(Vector4(1, 0, 0, 0), 
			Vector4(0, -1, 0, 0), 
			Vector4(0, 0, -1, 0), 
			Vector4(0, 0, 0, 1)));

		AutoPtr<TriBaseSurface> baseSurface(new TriBaseSurface(mesh));
		
		//baseSurface->prepare(TangentSpaceMode_MeshMender);
		//baseSurface->prepare(TangentSpaceMode_Lengyel);
		//baseSurface->prepare(TangentSpaceMode_Castano);
		baseSurface->prepare(TangentSpaceMode_UE3);

		m_mesh = baseSurface.release();
	}

	return true;
}

void BaseMeshPass::freeMesh()
{
	m_mesh = NULL;
}


void BaseMeshPass::initGeometryMap(Vector2::Arg extents)
{
	uint flags = GeometryImage::DisplacementFlag;

	if (m_displacements3D) {
		nvDebug(" Creating 3D displacements.\n");
		flags |= GeometryImage::VectorDisplacementFlag;
	}
	if (m_tangentSpace) {
		flags |= GeometryImage::NormalFlag;	// For tangent space normals!
	}


	m_geometryMap = new GeometryImage(flags, extents);
	m_geometryMask = 
		new BitMap(m_geometryMap->width(), m_geometryMap->height());
	m_geometryMask->clearAll();
}


void BaseMeshPass::rasterizeMesh(GeometryImage * detailedGeometryMap, 
  BitMap * detailedGeometryMask)
{
	nvCheck(m_mesh != NULL);
	nvDebugCheck(detailedGeometryMap != NULL);
	nvDebugCheck(detailedGeometryMask != NULL);

	Vector2 extents(float(m_geometryMap->width()), 
		float(m_geometryMap->height()));

	printf("--- Rasterizing base mesh\n");

	// Rasterize lores mesh.
	DisplacementPatchSampler displacementSampler(*detailedGeometryMap, 
		*detailedGeometryMask, *m_geometryMap, *m_geometryMask);

	Raster::SamplingCallback cb = DisplacementPatchSampler::sampleCallback;
	displacementSampler.setTangentSpace(m_tangentSpace);
	displacementSampler.setVectorDisplacement(m_displacements3D);
	displacementSampler.setCurrentSurface(m_mesh.ptr());

	const uint faceCount = m_mesh->faceCount();
	for (uint f = 0, percentage = 0; f < faceCount; f++)
	{
		if (percentage != (100 * f) / faceCount)
		{
			percentage = (100 * f) / faceCount;
			printf("\r%u%%", percentage);
			fflush(stdout);
		}

		m_mesh->selectFace(f);

		Vector2 texcoords[4];
		m_mesh->textureCoordinates(texcoords);
		for (uint k = 0; k < 4; k++) {
			texcoords[k] *= extents;
		}

		if (m_mesh->domain() == FaceDomain::Triangle)
		{
			Raster::drawTriangle(RASTER_ANTIALIAS, 
				extents, texcoords, cb, &displacementSampler);
		}
		else
		{
			nvDebugCheck(m_mesh->domain() == FaceDomain::Quad);
			Raster::drawQuad(RASTER_ANTIALIAS, 
				extents, texcoords, cb, &displacementSampler);
		}
	}
	printf("\r");
}


void BaseMeshPass::applyFilters()
{
	nvDebug("--- Applying filters.\n");

	m_geometryMap->fillSubpixels(m_geometryMask.ptr(), 0.001f);

	fillQuadraticExtrapolate(4, m_geometryMap->img(), m_geometryMask.ptr());
	//fillVoronoi(m_geometryMap->img(), m_geometryMask.ptr());
	//fillPullPush(m_geometryMap->img(), m_geometryMask.ptr());
//	fillPullPushLinear(m_geometryMap->img(), m_geometryMask.ptr());

	if (!m_displacements3D)
	{
		m_geometryMap->img()->scaleBias(m_geometryMap->displacementChannel(), 1, 1.0f, 0.5f);
	}

	if (m_tangentSpace)
	{
		m_geometryMap->img()->normalize(m_geometryMap->normalChannel());
		m_geometryMap->img()->packNormals(m_geometryMap->normalChannel());
	}
}


bool BaseMeshPass::saveDisplacementMap(const char * fileName) const
{
	nvDebug("--- Saving displacement map.\n");

	int dspChannels = 1;
	if (m_displacements3D) dspChannels = 3;

	if (!ImageIO::saveFloat(fileName, m_geometryMap->img(), m_geometryMap->displacementChannel(), dspChannels))
	{
		printf("Error, cannot save '%s'\n", fileName);
		return false;
	}

	return true;
}

bool BaseMeshPass::saveNormalMap(const char * fileName) const
{
	nvDebug("--- Saving normal map.\n");

	if (!ImageIO::saveFloat(fileName, m_geometryMap->img(), m_geometryMap->normalChannel(), 3))
	{
		printf("Error, cannot save '%s'\n", fileName);
		return false;
	}

	return true;
}

void BaseMeshPass::saveMaps(const char * name) const
{
	Path fileName;
	
	if (m_tangentSpace)
	{
		// Tangent space normal map.
		fileName = name;
		fileName.stripExtension();
		fileName.append("_tnmap.tga");

		saveNormalMap(fileName.str());
	}


	if (m_displacements3D)
	{
		if (m_tangentSpace)
		{
			// Vector displacement in tangent space.
			fileName = name;
			fileName.stripExtension();
			fileName.append("_tvmap.tif");

			saveDisplacementMap(fileName.str());
		}
		else
		{
			// Vector displacement in object space.
			fileName = name;
			fileName.stripExtension();
			fileName.append("_vmap.tif");

			saveDisplacementMap(fileName.str());
		}
	}
	else
	{
		// 1D displacement only.
		fileName = name;
		fileName.stripExtension();
		fileName.append("_dmap.tif");

		saveDisplacementMap(fileName.str());
	}
}

