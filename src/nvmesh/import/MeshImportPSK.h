// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_PSK_H
#define NV_MESH_IMPORT_PSK_H

#include <nvmesh/import/MeshImport.h>
#include <nvmesh/animation/MeshSkeleton.h>

namespace nv
{

	/// PSK mesh importer.
	class MeshImportPSK : public MeshImport
	{
	public:
		virtual bool import(Stream * stream);
		
		const Skeleton * skeleton() const { return &m_skeleton; };
		
	private:
		Skeleton m_skeleton;
	};

} // nv namespace


#endif // NV_MESH_IMPORT_PSK_H
