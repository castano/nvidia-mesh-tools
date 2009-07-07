// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_OFF_H
#define NV_MESH_IMPORT_OFF_H

#include <nvmesh/import/MeshImport.h>

namespace nv
{

	/// OFF mesh importer.
	class MeshImportOFF : public MeshImport
	{
	public:
		virtual bool import(Stream * stream);
	};

} // nv namespace

#endif // NV_MESH_IMPORT_OBJ_H
