// Copyright NVIDIA Corporation 2008 -- Ignacio Castano <icastano@nvidia.com>

#ifndef NV_MESH_MAYAEXPORTER_H
#define NV_MESH_MAYAEXPORTER_H

#include <maya/MString.h>

class MArgList;
class MFileObject;


namespace nv
{

	// Export options.
	struct ExportOptions
	{
		ExportOptions();

		void setDefaults();

		void parse(const MArgList & args);
		void parse(const MString & optionsString);
		void parse(int argc, const char ** argv);

		void setFileName(const char * fileName);
		void setFileObject(const MFileObject & file);
		
		void setExportSelected(bool exportSelected);

		bool isValid() const;

	private:
		friend class Exporter;

		const char * filePath() const;

		int frameCount() const;
		float totalTime() const;

	private:

		bool m_exportSelected;

		bool m_exportReferencePose;
		bool m_exportVertexAnimation;

		MString m_filePath;

		// Animation export options:
		float m_referenceFrame;
		float m_beginTime;
		float m_endTime;
		float m_sampleFrequency;

	};


	// Exporter.
	class Exporter
	{
	public:

		MStatus doExport(const ExportOptions & options);

	private:

		// @@ Compute scene stats?

	};


} // nv namespace


#endif // NV_MESH_MAYAEXPORTER_H
