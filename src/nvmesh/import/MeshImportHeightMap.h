// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_HEIGHTMAP_H
#define NV_MESH_IMPORT_HEIGHTMAP_H

#include <nvmesh/import/MeshImport.h>

namespace nv
{

	/// Heightmap mesh importer.
	class MeshImportHeightMap : public MeshImport
	{
	public:
		virtual bool import(Stream * stream);
	};

} // nv namespace

#endif // NV_MESH_IMPORT_HEIGHTMAP_H
