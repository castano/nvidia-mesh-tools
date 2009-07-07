// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#include <nvcore/Ptr.h>
#include <nvcore/StdStream.h>

#include <nvmesh/QuadTriMesh.h>
#include <nvmesh/animation/MorphTarget.h>


// Required by maya!
#include <iostream>

#include <maya/MTime.h>
#include <maya/MAnimControl.h>
#include <maya/MGlobal.h>
#include <maya/MString.h>
#include <maya/MFileObject.h>
#include <maya/MArgList.h>

#include "MayaExporter.h"
#include "MayaMeshBuilder.h"
#include "MayaUtils.h"

using namespace nv;


ExportOptions::ExportOptions()
{
	setDefaults();
}

void ExportOptions::setDefaults()
{
	m_exportSelected = false;

	m_exportReferencePose = true;
	m_exportVertexAnimation = true;

	m_filePath = "";

	m_referenceFrame = -1.0f;
	m_beginTime = max(0, (float)MAnimControl::minTime().as(MTime::uiUnit()));
	m_endTime = (float)MAnimControl::maxTime().as(MTime::uiUnit());
	m_sampleFrequency = 0.5f;
}

void ExportOptions::parse(const MArgList & args)
{
	// Parse the arguments. 
	for (uint i = 0; i < args.length(); i++)
	{
		if (args.asString(i) == "-referenceFrame")
		{
			//m_referenceFrame = (float)args.asTime(++i).as(MTime::uiUnit());
			m_referenceFrame = (float)args.asDouble(++i);
		}
		else if (args.asString(i) == "-animationRange")
		{
			//m_beginTime = (float)args.asTime(++i).as(MTime::uiUnit());
			//m_endTime = (float)args.asTime(++i).as(MTime::uiUnit());
			m_beginTime = (float)args.asDouble(++i);
			m_endTime = (float)args.asDouble(++i);
		}
		else if (args.asString(i) == "-sampleFrequency")
		{
			//m_sarmpleFrequency = (float)args.asTime(++i).as(MTime::uiUnit());
			m_sampleFrequency = (float)args.asDouble(++i);
		}
		else if (args.asString(i) == "-fileName")
		{
			m_filePath = args.asString(++i);
		}
	}

	// @@ Make sure totalTime() > 0.
}

void ExportOptions::parse(const MString & optionsString)
{
	MStringArray stringArray;
	optionsString.split(' ', stringArray);

	MArgList args;
	for (uint i = 0; i < stringArray.length(); i++)
	{
		args.addArg(stringArray[i]);
	}

	parse(args);
}

void ExportOptions::parse(int argc, const char ** argv)
{
	MArgList args;
	for (int i = 0; i < argc; i++)
	{
		args.addArg(MString(argv[i]));
	}

	parse(args);
}

void ExportOptions::setFileName(const char * fileName)
{
	m_filePath = fileName;
}

void ExportOptions::setFileObject(const MFileObject & file)
{
	m_filePath = file.fullName();
//	file.resolvedName();
}

void ExportOptions::setExportSelected(bool exportSelected)
{
	m_exportSelected = exportSelected;
}

bool ExportOptions::isValid() const
{
	return m_filePath != "";
}

const char * ExportOptions::filePath() const
{
	return m_filePath.asChar();
}

int ExportOptions::frameCount() const
{
	// @@ How to round?
	return totalTime() / m_sampleFrequency;
}

float ExportOptions::totalTime() const
{
	MTime diff(m_endTime - m_beginTime, MTime::uiUnit());
	return diff.as(MTime::kSeconds);
}


MStatus Exporter::doExport(const ExportOptions & options)
{
	// Open stream.
	StdOutputStream stream(options.filePath());

	// Write header.
	char str[] = "NVMA";
	stream << str[0] << str[1] << str[2] << str[3];

	// Set reference pose.
	MayaTime time(options.m_referenceFrame);

	if (options.m_exportVertexAnimation)
	{
		// Export morph target mesh.
		MorphTargetMesh morphTargetMesh;

		// Set reference target.
		MayaMeshBuilderOptions buildOptions;
		buildOptions.addTexcoords = false;
		buildOptions.addNormals = false;

		MayaMeshBuilder builder(buildOptions);

		if (options.m_exportSelected)
		{
			builder.addSelection();
		}
		else
		{
			builder.addScene();
		}

		builder.done();

		AutoPtr<QuadTriMesh> mesh(builder.buildQuadTriMesh());
		morphTargetMesh.setReferenceTarget(*mesh);

		float frameTime = options.m_beginTime;
		const uint frameCount = options.frameCount();
		for (uint i = 0; i < frameCount; i++)
		{
			time.setCurrentTime(frameTime);

			builder.reset();

			if (options.m_exportSelected)
			{
				builder.addSelection();
			}
			else
			{
				builder.addScene();
			}
			
			builder.done();

			mesh = builder.buildQuadTriMesh();
			morphTargetMesh.addMorphTarget(*mesh);

			frameTime += options.m_sampleFrequency;
		}

		stream << morphTargetMesh;

		// Export morph target sequence.
		MorphTargetSequence morphTargetSequence;
		morphTargetSequence.setLinearSequence(options.frameCount(), options.totalTime());

		stream << morphTargetSequence;
	}

	return MS::kSuccess;
}


