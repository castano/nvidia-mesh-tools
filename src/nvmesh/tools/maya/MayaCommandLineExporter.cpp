// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <iostream>

#include <maya/MItDependencyNodes.h>
#include <maya/MLibrary.h>
#include <maya/MFileIO.h>

#include "MayaExporter.h"

using namespace nv;


int main(int argc, char* argv[])
{
	printf("Maya Command Line Exporter - Copyright NVIDIA Corporation 2008\n\n");

	if (argc != 2) {
		printf("Usage options: MayaCommandLineExporter 'fileName'\n");
		return 0;
	}

	char * pluginName = argv[0];
	const char * fileName = argv[1];

	// initialise the maya library
	MLibrary::initialize(pluginName);

	ExportOptions options;

	//options.parse(argc - 1, argv + 1);
	options.setFileName("out.nvmesh");

	if (!options.isValid())
	{
		printf("Invalid export options.\n");
	}
	else
	{
		if (MFileIO::open(fileName) == MS::kSuccess)
		{
			Exporter exporter;

			if (exporter.doExport(options) == MS::kSuccess)
			{
				printf("Export succeeded.\n");
			}
			else
			{
				printf("Export failed.\n");
			}
		}
	}

	printf("\n");

	// cleanup maya
	MLibrary::cleanup();

	return 0;
}