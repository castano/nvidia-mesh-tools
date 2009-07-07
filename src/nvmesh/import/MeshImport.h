// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_IMPORT_H
#define NV_MESH_IMPORT_H

#include <nvcore/Containers.h>
#include <nvmesh/MeshBuilder.h>

// This doesn't work on MSVC!
#define NV_REGISTER_MESH_IMPORT_FACTORY(Factory) \
	namespace { \
		static Factory factory##Factory; \
		struct Factory##Registrar { \
			Factory##Registrar() { MeshImportFactory::addFactory(&factory##Factory); } \
			~Factory##Registrar() { MeshImportFactory::removeFactory(&factory##Factory); } \
		}; \
		static Factory##Registrar registrar##Factory; \
	}

namespace nv
{
	class Stream;
	class Skeleton;

	/// Mesh importer.
	class MeshImport
	{
	public:
		virtual bool import(Stream * stream) = 0;
				
		const MeshBuilder & builder() const { return m_builder; };
		
		// @@ This is not very clean, the importer has ownership of the skeleton.
		virtual const Skeleton * skeleton() const { return NULL; };
		
		static MeshImport * importer(const char * fileName);

	protected:
		MeshBuilder m_builder;
	};

	/// Factory of mesh importers.
	class MeshImportFactory
	{
	public:
		virtual const char * extension() const = 0;
		virtual MeshImport * createImporter() const = 0;
		
	//	static MeshImport * findMeshImportByName(const char * name);
	//	static MeshImport * findMeshImportFor(const char * fileName);
		
		static void addFactory(const MeshImportFactory * factory);
		static void removeFactory(const MeshImportFactory * factory);
	};

} // nv namespace

#endif // NV_MESH_IMPORT_H
