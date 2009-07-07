// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_M_H
#define NV_MESH_IMPORT_M_H

#include <nvmesh/import/MeshImport.h>

namespace nv
{

	/// M mesh importer.
	class MeshImportM : public MeshImport
	{
	public:
		virtual bool import(Stream * stream);
	};

} // nv namespace

#endif // NV_MESH_IMPORT_M_H
