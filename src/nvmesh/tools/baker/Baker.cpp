// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

// --himesh pixolator/DispTestHead_NG_Low.OBJ --hisub 4 --hidisp pixolator/DispTestHead_Disp16.tif 0.834
// --himesh sword/ZSwordLowRes.OBJ --hisub 4 --hidisp sword/ZSwordDispMap.tif 0.286

// Avoids a warning when including CpuInfo about the attributes of ceil()
#include <math.h>

#include <nvcore/StdStream.h>
#include <nvcore/Ptr.h>
#include <nvcore/Timer.h>
#include <nvcore/CpuInfo.h>

#include <nvimage/FloatImage.h>
#include <nvimage/Image.h>
#include <nvimage/ImageIO.h>
#include <nvimage/HoleFilling.h>

#include <nvmesh/TriMesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdge.h>
#include <nvmesh/raster/Raster.h>

#include <nvmesh/subdiv/Subdivide.h>

#include <nvmesh/raster/ClippedTriangle.h>

#include <IlmThreadPool.h>

#include "Samplers.h"
#include "DetailedMeshPass.h"
#include "BaseMeshPass.h"
#include "CmdOptions.h"

using namespace nv;


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

struct Options : CmdOptions
{
	TCLAP::ValueArg<int> widthArg;
	TCLAP::ValueArg<int> heightArg;
	TCLAP::ValueArg<int> threadsArg;
	TCLAP::ValueArg<std::string> outputArg;

	Options() :
	widthArg("", "width",
		"Output textures width (defaults to 1024).", 
		false, 1024, "size"),
	heightArg("", "height",
		"Output textures height (defaults to 1024).", 
		false, 1024, "size"),
	threadsArg("", "threads",
		"Number of threads to use for processing (0 is auto, default).",
		false, 0, "integer"),
	outputArg("", "output",
		"Filepath prefix for the generated maps (defaults to \"output\").",
		false, "output", "prefix") {}

	void add(TCLAP::CmdLine &cmd) {
		cmd.add(outputArg);
		cmd.add(threadsArg);
		cmd.add(heightArg);
		cmd.add(widthArg);
	}
};


int main(int argc, char *argv[])
{
	printf("Texture Baker - Copyright NVIDIA Corporation 2008\n\n");

	MyMessageHandler messageHandler;
	debug::setMessageHandler(&messageHandler);

	//testIterator2D();

#if _DEBUG
    printf("Press ENTER to continue:");
    fflush(stdout);
    getchar();
#endif

	DetailedMeshPass detailedMeshPass;
	BaseMeshPass baseMeshPass;
	int width, height, numThreads;
	std::string outputPrefix;
	
	try {
		TCLAP::CmdLine cmd("Texture Baker - "
			"Copyright NVIDIA Corporation 2008",
			/*delimiter*/ ' ', /*version*/ "0.0.1");
		AutoPtr<CmdOptions> 
			optDetailed( detailedMeshPass.getCmdOptions() );
		AutoPtr<CmdOptions> 
			optBase( baseMeshPass.getCmdOptions() );
		Options opt;

		optDetailed->add(cmd);
		optBase->add(cmd);
		opt.add(cmd);

		cmd.parse(argc, argv);
		detailedMeshPass.setCmdOptions(optDetailed.ptr());
		baseMeshPass.setCmdOptions(optBase.ptr());
		width  = opt.widthArg.getValue();
		height = opt.heightArg.getValue();
		numThreads = opt.threadsArg.getValue();

		outputPrefix = opt.outputArg.getValue();
		nvDebugCheck( !outputPrefix.empty() );
	}
	catch (TCLAP::ArgException &) { 
		return 1; 
	}

	// Initializes the thread pool with the given number of threads
	if (numThreads < 1) {
		numThreads = CpuInfo::processorCount();
	}
	IlmThread::ThreadPool &pool = IlmThread::ThreadPool::globalThreadPool();
	pool.setNumThreads(numThreads);
	nvDebug("--- Creating %d-by-%d images with %d %s\n",
		width, height, numThreads, 
		(numThreads != 1 ? "threads" : "thread"));

	// detailed mesh pass ............
	Vector2 extents(width, height);
	detailedMeshPass.initGeometryMap(extents);

	if (detailedMeshPass.hasValidMesh())
	{
		if (!detailedMeshPass.loadMesh()) {
			printf("Error loading hi-res mesh.\n");
			return 0;
		};

		Timer timer;
		timer.start();

		detailedMeshPass.rasterizeMesh();

		printf("%.2f sec\n", float(timer.elapsed()) / 1000);

		detailedMeshPass.freeMesh(); // free mesh data ASAP, it's huge!

		detailedMeshPass.applyFilters();
		detailedMeshPass.saveMaps(outputPrefix.c_str());
	}
	else
	{
		if (!detailedMeshPass.loadGeometryMap())
		{
			printf("Error loading geometry map.\n");
		}
	}


	// base mesh pass............
	baseMeshPass.initGeometryMap(extents);
	
	if (!baseMeshPass.loadMesh()) {
		printf("Error loading lo-res mesh.\n");
		return 0;
	};

	baseMeshPass.rasterizeMesh(detailedMeshPass.geometryMap(), 
				   detailedMeshPass.geometryMask());
	baseMeshPass.freeMesh(); 
	baseMeshPass.applyFilters();
	baseMeshPass.saveMaps(outputPrefix.c_str());

	return 0;
}
