// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_PLY_H
#define NV_MESH_IMPORT_PLY_H

#include <nvmesh/import/MeshImport.h>

namespace nv
{

	/// OBJ mesh importer.
	class MeshImportPLY : public MeshImport
	{
	public:
		virtual bool import(Stream * stream);
	};

} // nv namespace

#endif // NV_MESH_IMPORT_PLY_H