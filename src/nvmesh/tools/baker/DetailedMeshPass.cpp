// Copyright NVIDIA Corporation 2006 -- Denis Kovacs <dkovacs@nvidia.com>

#include <nvcore/Ptr.h>
#include <nvcore/StdStream.h>

#include <nvmath/Montecarlo.h>

#include <nvimage/Image.h>
#include <nvimage/ImageIO.h>
#include <nvimage/Filter.h>

#include <nvmesh/TriMesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>
#include <nvmesh/raster/Raster.h>

#include <nvmesh/import/MeshImport.h>
#include <nvmesh/geometry/MeshNormals.h>
#include <nvmesh/geometry/Bounds.h>
#include <nvmesh/subdiv/Subdivide.h>
#include <nvmesh/weld/VertexWeld.h>
#include <nvmesh/kdtree/MeshKDTree.h>
#include <nvmesh/kdtree/KDTree.h>

#include "DetailedMeshPass.h"
#include "Samplers.h"
#include "TiledTask.h"

#include <string>

using namespace nv;

DetailedMeshPass::DetailedMeshPass() :
	m_mesh(NULL),
	m_geometryMap(NULL),
	m_geometryMask(NULL),
	m_meshFileName(""),
	m_subdivisionLevels(0),
	m_displacementFileName(""),
	m_nmapFileName(""),
	m_superSampling(false),
	m_outputPositionMap(false),
	m_outputOcclusionMap(false),
	m_outputBentNormalMap(false),
	m_rebuildCache(false),
	m_maxDistRatio(0.25f),
	m_useIlmOcclusion(false)
{
}

DetailedMeshPass::~DetailedMeshPass()
{
}

struct DetailedMeshPass::Options : CmdOptions {

	TCLAP::ValueArg<std::string> meshFileNameArg;
	TCLAP::ValueArg<int> subdivisionLevelsArg;
	TCLAP::ValueArg<std::string> displacementFileNameArg;
	TCLAP::ValueArg<std::string> nmapFileNameArg;
	TCLAP::SwitchArg superSamplingArg;
	TCLAP::SwitchArg outputPositionMapArg;
	TCLAP::SwitchArg outputOcclusionMapArg;
	TCLAP::SwitchArg outputBentNormalMapArg;
	TCLAP::ValueArg<int> raysArg;
	TCLAP::ValueArg<float> maxDistRatioArg;
	TCLAP::SwitchArg useIlmOcclusionArg;
	TCLAP::ValueArg<std::string> cacheDirArg;
	TCLAP::SwitchArg rebuildCacheArg;

	// The constructor initializes the options
	Options() :
	meshFileNameArg("", "himesh", 
		"High resolution mesh.", true, /*default*/ "", "filename"),
	subdivisionLevelsArg("", "hisub",
		"High resolution subdivision levels (default 0).", false, 0, "int"),
	displacementFileNameArg("",
		"hidisp", "Displacement file.", false, /*default*/ "", "filename"),
	nmapFileNameArg("", "hinmap",
		"Normal Map file", false, /*default*/ "", "filename"),
	superSamplingArg("", "hisupersampling",
		"Enable supersampling."),
	outputPositionMapArg("", "hipos",
		"Output position map."),
	outputOcclusionMapArg("", "occlusion",
		"Output occlusion map."),
	outputBentNormalMapArg("", "bentnormal",
		"Output bent normal map."),
	raysArg("", "rays", "Number of rays for ambient occlusion (default 64).", 
		false, 64, "int"),
	maxDistRatioArg("", "maxdistratio",
		"Ratio of the b.b. extents to which the max ray distance is set.",
		false, 0.25f, "ratio"),
	useIlmOcclusionArg("", "ilmocclusion",
		"Computes the occlusion as present by ILM in SIGGRAPH 2002, "
		"using binary instead of attenuated visibility."),
	cacheDirArg("", "hicachedir", 
		"Directory for saving and reading the caches.", 
		false, /*default*/ "", "dir"),
	rebuildCacheArg("", "rebuildcache",
		"Rebuild the mesh cache.")  {}

	void add(TCLAP::CmdLine &cmd) {
		cmd.add(rebuildCacheArg);
		cmd.add(cacheDirArg);
		cmd.add(useIlmOcclusionArg);
		cmd.add(maxDistRatioArg);
		cmd.add(raysArg);
		cmd.add(outputBentNormalMapArg);
		cmd.add(outputOcclusionMapArg);
		cmd.add(outputPositionMapArg);
		cmd.add(superSamplingArg);
		cmd.add(nmapFileNameArg);
		cmd.add(displacementFileNameArg);
		cmd.add(subdivisionLevelsArg);
		cmd.add(meshFileNameArg);
	}
};



CmdOptions * DetailedMeshPass::getCmdOptions() {
	return new Options();
}

void DetailedMeshPass::setCmdOptions(CmdOptions * opt) {
	nvDebugCheck(opt != NULL);
	Options & options = * dynamic_cast<Options*>(opt);

	nvDebugCheck( options.meshFileNameArg.isSet() );
	m_meshFileName = options.meshFileNameArg.getValue();
	
	if (options.displacementFileNameArg.isSet()) {
		m_displacementFileName = options.displacementFileNameArg.getValue();
	}
	if (options.nmapFileNameArg.isSet()) {
		m_nmapFileName = options.nmapFileNameArg.getValue();
	}

	m_subdivisionLevels   = options.subdivisionLevelsArg.getValue();
	m_superSampling       = options.superSamplingArg.getValue();
	m_outputPositionMap   = options.outputPositionMapArg.getValue();
	m_outputOcclusionMap  = options.outputOcclusionMapArg.getValue();
	m_outputBentNormalMap = options.outputBentNormalMapArg.getValue();
	m_rebuildCache        = options.rebuildCacheArg.getValue();
	m_numRays             = options.raysArg.getValue();
	m_maxDistRatio        = options.maxDistRatioArg.getValue();
	m_useIlmOcclusion     = options.useIlmOcclusionArg.getValue();

	// Sets the base path for the caches (path without the extension)
	const std::string & cacheDir = options.cacheDirArg.getValue();
	if ( !cacheDir.empty() ) {
		m_cacheBasePath = cacheDir.c_str();
		if (cacheDir[cacheDir.length()-1] != '/' &&
			cacheDir[cacheDir.length()-1] != '\\') {
			m_cacheBasePath.append("/");
		}
		m_cacheBasePath.translatePath();
		m_cacheBasePath.append( Path::fileName(m_meshFileName.c_str()) );
		m_cacheBasePath.stripExtension();
	}
	else {
		m_cacheBasePath = m_meshFileName.c_str();
		m_cacheBasePath.translatePath();
		m_cacheBasePath.stripExtension();
	}

}

bool DetailedMeshPass::hasValidMesh() const
{
	return !m_meshFileName.empty();
}


bool DetailedMeshPass::loadCachedMesh()
{
	nvDebugCheck( m_cacheBasePath.length() > 0 );
	
	Path cacheFileName(m_cacheBasePath);
	cacheFileName.append(".mesh_cache");

	StdInputStream stream(cacheFileName);
	if (!stream.isError())
	{
		m_mesh = new QuadTriMesh();
		stream << *m_mesh;
		return true;
	}

	return false;
}

void DetailedMeshPass::saveCachedMesh() const
{
	nvCheck(m_mesh != NULL);
	nvDebugCheck( m_cacheBasePath.length() > 0 );
	
	Path cacheFileName(m_cacheBasePath);
	cacheFileName.append(".mesh_cache");

	StdOutputStream stream(cacheFileName);
	if (!stream.isError())
	{
		stream << *m_mesh;
	}
}

bool DetailedMeshPass::loadCachedTree()
{
	nvDebugCheck( m_cacheBasePath.length() > 0 );
	
	Path cacheFileName(m_cacheBasePath);
	cacheFileName.append(".tree_cache");

	StdInputStream stream(cacheFileName);
	if (!stream.isError())
	{
		m_kdTree = new KDTree();
		stream << *m_kdTree;
		return true;
	}

	return false;
}

void DetailedMeshPass::saveCachedTree() const
{
	nvCheck(m_kdTree != NULL);
	nvDebugCheck( m_cacheBasePath.length() > 0 );
	
	Path cacheFileName(m_cacheBasePath);
	cacheFileName.append(".tree_cache");

	StdOutputStream stream(cacheFileName);
	if (!stream.isError())
	{
		stream << *m_kdTree;
	}
}



bool DetailedMeshPass::loadMesh()
{
	if (m_meshFileName.empty()) return false;

	if (m_rebuildCache || !loadCachedMesh())
	{
		nvDebug("Importing '%s'\n", m_meshFileName.c_str());
		AutoPtr<MeshImport> 
			importer(MeshImport::importer(m_meshFileName.c_str()));

		if (importer == NULL) {
			nvDebug("Error, unkown file type '%s'\n", 
				m_meshFileName.c_str());
			return false;
		}

		StdInputStream stream(m_meshFileName.c_str());
		if (stream.isError()) {
			nvDebug("Error, cannot open '%s'\n", 
				m_meshFileName.c_str());
			return false;
		}

		nvDebug("Loading mesh.\n");
		importer->import(&stream);

		m_mesh = importer->builder().buildQuadTriMesh();

		// Print mesh info:
		printf("--- Loaded mesh with %d faces and %d vertices.\n", 
			m_mesh->faceCount(), m_mesh->vertexCount());

        Box box = MeshBounds::box(m_mesh.ptr());
        printf("---  bounds: (%f %f %f)-(%f %f %f)\n", 
            box.minCorner().x(), box.minCorner().y(), box.minCorner().z(), 
            box.maxCorner().x(), box.maxCorner().y(), box.maxCorner().z());
		fflush(stdout);


		// if no subdivision:
		if (m_subdivisionLevels > 0)
		{
            printf("--- Subdividing mesh:\n");

			AutoPtr<HalfEdge::Mesh> mesh(importer->builder().buildHalfEdgeMesh());

			/*
			FloatImage * dispMap = NULL;
			if (m_displacementName != NULL) {
				dispMap = ImageIO::loadFloatTIFF(m_displacementName);
			}
			
			if (dispMap != NULL) {
				// Quantize texture coordinates.
				adjustSeamTexCoords(mesh, dispMap->width(), dispMap->height());
			}
			*/
			// Subdivide.
			for(int i = 0; i < m_subdivisionLevels; i++) {
                printf("---   Level: %d\n", i+1);
				mesh = Subdivide::doCatmullClarkSplit(mesh.ptr());
			}
			
			// Compute limit surface normals.
			MeshNormals::computeCatmullClarkNormals(mesh.ptr());
			
			/*
			// Displace.
			if (dispMap != NULL)
			{
				dispMap->scaleBias(0, 1, 8*m_displacementScale, -m_displacementBias);
				Subdivide::displace(mesh, dispMap, 0, Vector2(-0.5, -0.5));
			}
			
			// convert half edge mesh to tri mesh. 
			*/
		
			m_mesh = mesh->toQuadTriMesh();
		}

		// Weld vertices.
		WeldVertices(m_mesh.ptr());

		// Compute mesh normals.
		MeshNormals::computeNormals(m_mesh.ptr(), 
			WeightFaceArea /*| WeightFaceAngle*/);

		// Save mesh to speed up loading next time.
		saveCachedMesh();
	}

	if (m_outputOcclusionMap | m_outputBentNormalMap)
	{
		if (m_rebuildCache || !loadCachedTree())
		{
			printf("Building KD-Tree:\n");
			m_kdTree = buildKDTree(m_mesh.ptr());

			saveCachedTree();
		}
	}

	return true;
}

void DetailedMeshPass::freeMesh()
{
	m_mesh = NULL;
	m_kdTree = NULL;
}


void DetailedMeshPass::initGeometryMap(Vector2 extents)
{
	if (m_superSampling) {
		//nvDebug(" Supersampling detailed mesh.\n");
		extents *= 2;
	}

	uint flags = GeometryImage::PositionFlag | GeometryImage::NormalFlag;

	if (m_outputOcclusionMap)
	{
		flags |= GeometryImage::OcclusionFlag;
	}
	if (m_outputBentNormalMap)
	{
		flags |= GeometryImage::BentNormalFlag;
	}

	m_geometryMap = new GeometryImage(flags, extents);

	m_geometryMask =
		new BitMap(m_geometryMap->width(), m_geometryMap->height());
	m_geometryMask->clearAll();
}

bool DetailedMeshPass::loadGeometryMap()
{
	if ( !m_nmapFileName.empty() )
	{
		AutoPtr<FloatImage> 
			normalMap( ImageIO::loadFloat(m_nmapFileName.c_str()) );
		if (normalMap == NULL) {
			// @@ Do this in ImageIO::loadFloat?
			AutoPtr<Image> img( ImageIO::load(m_nmapFileName.c_str()) );
			normalMap = new FloatImage(img.ptr());
		}
		if (normalMap == NULL) {
			return false;
		}

		uint w = m_geometryMap->width();
		uint h = m_geometryMap->height();

		if (normalMap->width() != w || normalMap->height() != h)
		{
			if (normalMap->width() > w && normalMap->height() > h)
			{
				BoxFilter filter;
				normalMap = normalMap->resize(filter, w, h, FloatImage::WrapMode_Clamp);
			}
			else
			{
				// @@ We need to add support for general resize.
				return false;
			}
		}

		m_geometryMap->setNormalMap(normalMap.ptr());
		m_geometryMask = m_geometryMap->getBitMaskFromCoverage();
	}
	else
	{
		return false;
	}

	return true;
}

void DetailedMeshPass::rasterizeMesh()
{
	nvCheck(m_mesh != NULL);

	Vector2 extents(float(m_geometryMap->width()), 
		float(m_geometryMap->height()));

	nvDebug("--- Rasterizing detailed mesh\n");

	GeometrySampler geometrySampler(*m_geometryMap, *m_geometryMask);

	const uint faceCount = m_mesh->faceCount();

	Vector3 positions[4];
	Vector3 normals[4];
	Vector2 texcoords[4];

	for(uint f = 0, percentage = 0; f < faceCount; f++)
	{
		if (percentage != (100 * f) / faceCount) 
		{
			percentage = (100 * f) / faceCount;
			printf("\r%u%%", percentage);
			fflush(stdout);
		}

		bool isFaceQuad = m_mesh->isQuadFace(f);
		uint vtxCount = 3 + isFaceQuad;

		for (uint k = 0; k < vtxCount; k++) 
		{
			positions[k] = m_mesh->faceVertex(f, k).pos;
			texcoords[k] = scale(m_mesh->faceVertex(f, k).tex, extents);
			normals[k] = m_mesh->faceVertex(f, k).nor;
		}

		geometrySampler.setCurrentFace(vtxCount, positions, normals);

		if (isFaceQuad) 
		{
			Raster::drawQuad(RASTER_ANTIALIAS, extents, texcoords, 
				GeometrySampler::sampleQuadCallback, &geometrySampler);
		}
		else
		{
			Raster::drawTriangle(RASTER_ANTIALIAS, extents, texcoords, 
				GeometrySampler::sampleTriCallback, &geometrySampler);
		}

	}
	printf("\r");

	if (m_outputOcclusionMap | m_outputBentNormalMap) {
		computeOcclusionMap();
	}

	if (m_superSampling)
	{
		nvDebug("--- Downsampling\n");
		m_geometryMap = m_geometryMap->fastDownSample();
		m_geometryMask = m_geometryMap->getBitMaskFromCoverage();
	}
}



/// Occlusion task classes

class DetailedMeshPass::OcclusionTask : public TiledTask {

	// The current pass
	DetailedMeshPass &m_pass;

	// Computed max distance
	const float m_maxDist;

public:
	OcclusionTask(DetailedMeshPass &pass, float maxDist, 
		const Tile & tile, TiledTaskGenerator &generator) :
	  TiledTask(tile, generator),
	  m_pass(pass),
	  m_maxDist(maxDist)
	{}

	virtual void execute() 
	{
		// Calculates the occlusion at each pixel
		//HilbertIterator it(m_tile.rows(), m_tile.columns());
		//ScanlineIterator it(m_tile.rows(), m_tile.columns());
		MortonIterator it(m_tile.rows(), m_tile.columns());
		uint count = 0;
		Vector3 bentNormal;
		while(it.hasNext()) {
			const Index2D index = indexGlobal(it.next());
			const uint &x = index.column;
			const uint &y = index.row;
			const uint counter = index.index(m_pass.m_geometryMap->width());

			if (m_pass.m_geometryMask->bitAt(x, y))
			{
				float coverage = m_pass.m_geometryMap->pixel(x, y, 
					m_pass.m_geometryMap->coverageChannel());

				Vector3 p = m_pass.m_geometryMap->position(x, y) / coverage;
				Vector3 n = 
					normalize(m_pass.m_geometryMap->normal(x, y) / coverage, 
							  0.0f);

				float occlusion = m_pass.sampleOcclusion(p, n, m_maxDist, x, y, 
					counter, &bentNormal);
				nvCheck(occlusion >= 0.0f && occlusion <= 1.0f);
				
				m_pass.m_geometryMap->setPixel(coverage * occlusion, x, y, 
					m_pass.m_geometryMap->occlusionChannel());

				// @@ Read this from here, or add a parameter?
				if (m_pass.m_outputBentNormalMap) {
					m_pass.m_geometryMap->setPixel(bentNormal, x, y,
						m_pass.m_geometryMap->bentNormalChannel());

				}
			}

			++count;
		}
		nvDebugCheck(count == m_tile.rows() * m_tile.columns());
	}
};


class DetailedMeshPass::OcclusionTaskGenerator : public TiledTaskGenerator {

private:
	DetailedMeshPass & m_pass;
	float m_maxDist;

public:
	OcclusionTaskGenerator(DetailedMeshPass & pass, 
	float maxDistRatio = 0.25f,uint8 tileSize = 32) :
	  TiledTaskGenerator(pass.m_geometryMap->width(), 
		  pass.m_geometryMap->height(), tileSize), 
	  m_pass(pass) {

		Box bounds = m_pass.m_kdTree->bounds();
		Vector3 extents = bounds.extents();
		nvCheck( maxDistRatio > 0.0f );
		m_maxDist = maxDistRatio * 
			max(extents.x(), max(extents.y(), extents.z()));
	}

	virtual TiledTask * newTask(const Tile & tile) {
		return new OcclusionTask(m_pass, m_maxDist, tile, *this);
	}
};

/// END Occlusion task classes


void DetailedMeshPass::computeOcclusionMap()
{
	nvDebugCheck(m_geometryMap != NULL);
	nvDebugCheck(m_geometryMask != NULL);

	nvDebug("--- Computing occlusion maps\n");

	// Try to add a task per tile using the global thread pool
	OcclusionTaskGenerator tileGenerator(*this, m_maxDistRatio, 32);
	tileGenerator.runTasks();

	uint percentage = 0;
	for (uint i=0, numTotal=tileGenerator.numTiles(); i<numTotal; ++i) {
		
		Tile tile = tileGenerator.take();	// Blocks until a task has finished
		if (percentage != (100 * i) / numTotal)
		{
			percentage = (100 * i) / numTotal;
			printf("\r%u%%", percentage);
			fflush(stdout);
		}
	}
	printf("\r");
}


// @@ Note that while the occlusion is normalized [0,1] the bent normal
//    isn't: it needs to be normalized and packet. This is done at the
//    end while applying the texture filters.
float DetailedMeshPass::sampleOcclusion(Vector3::Arg pos, Vector3::Arg dir, 
					    float dist, uint x, uint y, uint counter,
						Vector3 * bentNormal) const
{
	Ray ray;
	ray.id = counter;
	ray.dir = dir;
	ray.origin = pos;
	ray.maxt = dist;
	ray.idir = Vector3(1.0f/ray.dir.x(), 1.0f/ray.dir.y(), 1.0f/ray.dir.z());

	// @@ Allocate only once.
	KDTree::Cache * cache = m_kdTree->createCache();
	m_kdTree->initSphericalCache(ray, cache);

	// @@ Precompute distributions, and select one according to x & y.
	// The code below depends on choosing a Cosine-Weighted distribution.
	nvDebugCheck(m_numRays > 0);
	SampleDistribution samples(m_numRays);
	samples.redistribute(SampleDistribution::Method_NRook, 
		SampleDistribution::Distribution_CosineHemisphere);

	float occlusion = 0.0f;
	Vector3 dirSum(zero);

	Basis basis;
	basis.buildFrameForDirection(dir);

	const uint sampleCount = samples.count();
	for (uint i = 0; i < sampleCount; i++)
	{
		Vector3 sampleDir = samples.sampleDir(i);

		// Align sample direction with surface.
		Vector3 alignedSampleDir = basis.transform(sampleDir);
		nvDebugCheck(dot(alignedSampleDir, dir) >= 0.0f);

		// Init ray.
		ray.id   = counter + i;
		ray.dir  = alignedSampleDir;
		ray.idir = 
			Vector3(1.0f/ray.dir.x(), 1.0f/ray.dir.y(), 1.0f/ray.dir.z());

		Hit hit;
		m_kdTree->testRay(*cache, ray, &hit);

		// The integrand is Cos[Theta]/Pi, and the probability of each sample
		// is also Cos[Theta]/Pi because of the cosine weighted distribution,
		// so our occlusion estimate becomes 1/N*[Sum of the visibility].
		// The original ILM definition uses binary visibility, but we default
		// to an expotentially attenuated visibility to get smoother results.
		if (m_useIlmOcclusion) {
			if (hit.t >= ray.maxt) {
				occlusion += 1.0f;
				dirSum    += alignedSampleDir;
			}
		}
		else {
			const float attenuation = powf(hit.t / ray.maxt, 32);
			occlusion        += clamp(attenuation, 0.0f, 1.0f);
			alignedSampleDir *= clamp(attenuation, 0.0f, 1.0f);
			dirSum           += alignedSampleDir;
		}
	}

	m_kdTree->deleteCache(cache);

	if (bentNormal != NULL) {
		*bentNormal = dirSum;
	}

	return occlusion / float(sampleCount);
}


void DetailedMeshPass::applyFilters()
{
	nvDebug("--- Applying filters.\n");

	m_geometryMap->fillSubpixels(m_geometryMask.ptr(), 0.1f);

	//fillExtrapolate(4, m_geometryMap->img(), m_geometryMask.ptr());
	fillQuadraticExtrapolate(4, m_geometryMap->img(), m_geometryMask.ptr(), 
		m_geometryMap->coverageChannel());
	//fillVoronoi(m_geometryMap->img(), m_geometryMask.ptr());
	//fillPullPush(m_geometryMap->img(), m_geometryMask.ptr());
	//fillPullPushLinear(m_geometryMap->img(), m_geometryMask.ptr());

	m_geometryMap->img()->normalize(m_geometryMap->normalChannel());
	m_geometryMap->img()->packNormals(m_geometryMap->normalChannel());

	if (m_outputBentNormalMap) {
		m_geometryMap->img()->normalize(m_geometryMap->bentNormalChannel());
		m_geometryMap->img()->packNormals(m_geometryMap->bentNormalChannel());
	}
}

bool DetailedMeshPass::saveNormalMap(const char * normalMapFile) const
{
	nvDebug("--- Saving normal map.\n");

	if (!ImageIO::saveFloat(normalMapFile, m_geometryMap->img(), 
		                    m_geometryMap->normalChannel(), 3))
	{
		printf("Error, cannot save '%s'\n", normalMapFile);
		return false;
	}

	return true;
}

bool DetailedMeshPass::savePositionMap(const char * positionMapFile) const
{
	nvDebug("--- Saving position map.\n");

	if (!ImageIO::saveFloat(positionMapFile, m_geometryMap->img(), 
		                    m_geometryMap->positionChannel(), 3))
	{
		printf("Error, cannot save '%s'\n", positionMapFile);
		return false;
	}

	return true;
}

bool DetailedMeshPass::saveOcclusionMap(const char * occlusionMapFile) const
{
	nvDebug("--- Saving occlusion map.\n");

	if (!ImageIO::saveFloat(occlusionMapFile, m_geometryMap->img(), 
		                    m_geometryMap->occlusionChannel(), 1))
	{
		printf("Error, cannot save '%s'\n", occlusionMapFile);
		return false;
	}

	return true;
}

bool DetailedMeshPass::saveBentNormalMap(const char * bentNormalMapFile) const
{
	nvDebug("--- Saving bent normal map.\n");

	if (!ImageIO::saveFloat(bentNormalMapFile, m_geometryMap->img(), 
		                    m_geometryMap->bentNormalChannel(), 3))
	{
		printf("Error, cannot save '%s'\n", bentNormalMapFile);
		return false;
	}

	return true;
}


void DetailedMeshPass::saveMaps(const char * name) const
{
	Path fileName;

	// Object space normal map.
	fileName = name;
	fileName.stripExtension();
	fileName.append("_nmap.tif");

	saveNormalMap(fileName.str());
	
	if (m_outputPositionMap) 
	{
		fileName = name;
		fileName.stripExtension();
		fileName.append("_pmap.tif");

		savePositionMap(fileName.str());
	}

	if (m_outputOcclusionMap)
	{
		fileName = name;
		fileName.stripExtension();
		fileName.append("_omap.tga");

		saveOcclusionMap(fileName.str());
	}

	if (m_outputBentNormalMap)
	{
		fileName = name;
		fileName.stripExtension();
		fileName.append("_bnmap.tga");

		saveBentNormalMap(fileName.str());
	}
}
