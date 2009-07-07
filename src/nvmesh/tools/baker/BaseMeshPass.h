// Copyright NVIDIA Corporation 2006 -- Denis Kovacs <dkovacs@nvidia.com>

#ifndef NV_BAKER_BASEMESHPASS_H
#define NV_BAKER_BASEMESHPASS_H

#include <nvcore/Ptr.h>
#include <nvmath/Vector.h>
#include <nvimage/HoleFilling.h>

#include "GeometryImage.h"
#include "CmdOptions.h"

namespace nv
{
	class BaseSurface;

	class BaseMeshPass : public CmdOptionsProvider
	{
	public:
		BaseMeshPass();

		~BaseMeshPass();

		virtual CmdOptions * getCmdOptions();
		virtual void setCmdOptions(CmdOptions * opt);

		bool loadMesh();
		void freeMesh();

		void initGeometryMap(Vector2::Arg extents);
		const GeometryImage * geometryMap() const { return m_geometryMap.ptr(); }

		void rasterizeMesh(GeometryImage * detailedGeometryMap, BitMap * detailedGeometryMask);
		void applyFilters();

		void saveMaps(const char * name) const;

	private:

		bool saveDisplacementMap(const char * displacementMapFile) const;
		bool saveNormalMap(const char * fileName) const;

	private:
		AutoPtr<BaseSurface> m_mesh;
		AutoPtr<GeometryImage> m_geometryMap;
		AutoPtr<BitMap> m_geometryMask;

		// Options:
		struct Options;

		std::string m_meshFileName;
		int m_subdivisionLevels;
		int m_baseMeshType;
		bool m_displacements3D;
		bool m_tangentSpace;
		

	};
};

#endif // NV_BAKER_BASEMESHPASS_H
