// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_OBJ_H
#define NV_MESH_IMPORT_OBJ_H

#include <nvmesh/import/MeshImport.h>

namespace nv
{

	/// OBJ mesh importer.
	class MeshImportOBJ : public MeshImport
	{
	public:
		virtual bool import(Stream * stream);
	};

} // nv namespace

#endif // NV_MESH_IMPORT_OBJ_H
