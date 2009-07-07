// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvcore/nvcore.h>

// Required by maya!
#include <iostream>

#include <maya/MPxFileTranslator.h>
#include <maya/MPxCommand.h>

#include <maya/MFnPlugin.h> 

#include "MayaExporter.h"

using namespace nv;


// Export file translator.
class ExportFileTranslator : public MPxFileTranslator
{
public:
	
	bool haveReadMethod() const { return false; }
	MStatus reader(const MFileObject& file, const MString& optionsString, FileAccessMode mode)
	{
		fprintf(stderr, "ExportFileTranslator::reader not implemented\n");
		return MS::kFailure;
	}

	bool haveWriteMethod() const { return true; }

	MStatus writer(const MFileObject & file, const MString & optionsString, FileAccessMode mode)
	{
		ExportOptions options;

		options.parse(optionsString);
		options.setFileObject(file);
		options.setExportSelected(mode == MPxFileTranslator::kExportActiveAccessMode);

		Exporter exporter;

		return exporter.doExport(options);
	}

	MString defaultExtension() const
	{
		return "nvmesh";
	}

    MFileKind identifyFile(const MFileObject & file, const char * buffer, short size) const
	{
		// Check extension.
		MString fileName = file.name();

		int start = fileName.rindex('.');
		int end = fileName.numChars();
		MString extension = fileName.substring(start, end);

		
		if (extension == ".nvmesh") // @@ "nvmesh" ?
		{
			// kIsMyFileType
			return kCouldBeMyFileType;
		}
		else
		{
			return kNotMyFileType;
		}
	}

	static void * creator() { return new ExportFileTranslator; }
};


// Export command.
class ExportCommand : public MPxCommand
{
public:
	MStatus doIt( const MArgList & args );
	static void * creator() { return new ExportCommand; }
};


MStatus ExportCommand::doIt( const MArgList & args )
{
	ExportOptions options;

	options.parse(args);

	if (!options.isValid())
	{
		fprintf(stderr, "Invalid export options\n");
		return MS::kFailure;
	}

	Exporter exporter;

	return exporter.doExport(options);
}


DLL_EXPORT MStatus initializePlugin(MObject obj)
{
	MFnPlugin plugin(obj, "NVIDIA Corp.", "1.0", "Any");

	MStatus status = plugin.registerCommand("nvExport", ExportCommand::creator);
	if (status != MS::kSuccess) {
		status.perror("initializePlugin");
	}

	status = plugin.registerFileTranslator("NVIDIA Exporter", "none", ExportFileTranslator::creator, "NvMayaExportScript");
	if (status != MS::kSuccess) {
		status.perror("initializePlugin");
	}

	return MS::kSuccess; 
}

DLL_EXPORT MStatus uninitializePlugin(MObject obj)
{
	MFnPlugin plugin( obj );

	MStatus status = plugin.deregisterCommand("nvExport");
	if (status != MS::kSuccess) {
		status.perror("uninitializePlugin");
	}

	status = plugin.deregisterFileTranslator("NVIDIA Exporter");
	if (status != MS::kSuccess) {
		status.perror("uninitializePlugin");
	}

	return MS::kSuccess;
}

