// Copyright NVIDIA Corporation 2006 -- Denis Kovacs <dkovacs@nvidia.com>

#ifndef NV_BAKER_DETAILEDMESHPASS_H
#define NV_BAKER_DETAILEDMESHPASS_H

#include <nvcore/StrLib.h>
#include <nvcore/Array2D.h>

#include <nvmesh/QuadTriMesh.h>
#include <nvmesh/kdtree/MeshKDTree.h>

#include "GeometryImage.h"
#include "CmdOptions.h"

#include <string>
#include <iostream>

namespace nv
{


	class DetailedMeshPass : public CmdOptionsProvider
	{
	public:
		DetailedMeshPass();
		~DetailedMeshPass();

		virtual CmdOptions * getCmdOptions();
		virtual void setCmdOptions(CmdOptions * opt);

		bool hasValidMesh() const;

		bool loadMesh();
		void freeMesh();

		void initGeometryMap(Vector2 extents);

		bool loadGeometryMap();

		QuadTriMesh * getMesh() { return m_mesh.ptr(); }
		GeometryImage * geometryMap() { return m_geometryMap.ptr(); }
		BitMap * geometryMask() { return m_geometryMask.ptr(); } 

		void rasterizeMesh();
		void computeOcclusionMap();
		void applyFilters();
		void saveMaps(const char * name) const;

	private:

		float sampleOcclusion(Vector3::Arg pos, Vector3::Arg dir, 
			float dist, uint x, uint y, uint counter,
			Vector3 * bentNormal = NULL) const;

		bool saveNormalMap(const char * normalMapFile) const;
		bool savePositionMap(const char * positionMapFile) const;
		bool saveOcclusionMap(const char * occlusionMapFile) const;
		bool saveBentNormalMap(const char * bentNormalMapFile) const;

		bool loadCachedMesh();
		void saveCachedMesh() const;

		bool loadCachedTree();
		void saveCachedTree() const;


	private:

		AutoPtr<QuadTriMesh> m_mesh;
		AutoPtr<GeometryImage> m_geometryMap;
		AutoPtr<BitMap> m_geometryMask;
		AutoPtr<KDTree> m_kdTree; 
		
		// Options:
		struct Options;

		std::string m_meshFileName;
		int m_subdivisionLevels;
		std::string m_displacementFileName;
		std::string m_nmapFileName;
		bool m_superSampling;
		bool m_outputPositionMap;
		bool m_outputOcclusionMap;
		bool m_outputBentNormalMap;
		bool m_rebuildCache;
		int m_numRays;
		float m_maxDistRatio;
		Path m_cacheBasePath;
		bool m_useIlmOcclusion;

		// Abient occlusion task classes
		class OcclusionTask;
		class OcclusionTaskGenerator;
	};

}




#endif // NV_BAKER_DETAILEDMESHPASS_H
