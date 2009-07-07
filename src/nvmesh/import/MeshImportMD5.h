// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_MD5_H
#define NV_MESH_IMPORT_MD5_H

#include <nvmesh/import/MeshImport.h>
#include <nvmesh/animation/MeshSkeleton.h>

namespace nv
{
	class Tokenizer;

	/// MD5 mesh importer.
	class MeshImportMD5 : public MeshImport
	{
	public:
		virtual bool import(Stream * stream);
		
		const Skeleton * skeleton() const { return &m_skeleton; };

	private:
		bool parseJoints(Tokenizer &);
		bool parseMesh(Tokenizer &);
		
	private:
		Skeleton m_skeleton;
	};

} // nv namespace


#endif // NV_MESH_IMPORT_MD5_H
