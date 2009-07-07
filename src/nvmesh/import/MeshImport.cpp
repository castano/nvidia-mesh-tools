
#include "MeshImport.h"

#include "MeshImportM.h"
#include "MeshImportOBJ.h"
#include "MeshImportMD5.h"
#include "MeshImportPSK.h"
#include "MeshImportOFF.h"

using namespace nv;

namespace {

	static Array<const MeshImportFactory *> * s_factoryArray = NULL;

} // nanmespace


// static
//MeshImport * findMeshImportByName(const char * name);

// static 
/*MeshImport * findMeshImportFor(const char * fileName)
{
	foreach(i, m_factoryArray) {
		MeshImportFactory * factory = m_factoryArray[i];
		nvCheck(factory != NULL);
		
		if (factory->canImport(filName)) {
			return factory->createMeshImport();
		}
	}
	
	return false;
}
*/

// static 
MeshImport * MeshImport::importer(const char * fileName)
{
	const char * extension = Path::extension(fileName);

	if (strCaseCmp(extension, ".obj") == 0) {
		return new MeshImportOBJ();
	}
	if (strCaseCmp(extension, ".m") == 0) {
		return new MeshImportM();
	}
	if (strCaseCmp(extension, ".md5mesh") == 0) {
		return new MeshImportMD5();
	}
	if (strCaseCmp(extension, ".psk") == 0) {
		return new MeshImportPSK();
	}
	if (strCaseCmp(extension, ".off") == 0) {
		return new MeshImportOFF();
	}

	if (s_factoryArray == NULL) {
		return NULL;
	}

	foreach(i, *s_factoryArray) {
		const MeshImportFactory * factory = (*s_factoryArray)[i];
		nvCheck(factory != NULL);
		
		if (strCaseCmp(extension, factory->extension()) == 0)
		{
			return factory->createImporter();
		}
	}

	return NULL;
}



// static 
void MeshImportFactory::addFactory(const MeshImportFactory * factory)
{
	nvCheck(factory != NULL);
	
	if( s_factoryArray == NULL ) {
		s_factoryArray = new Array<const MeshImportFactory *>(8);
	}
	s_factoryArray->pushBack(factory);
}

// static
void MeshImportFactory::removeFactory(const MeshImportFactory * factory)
{
	nvCheck(factory != NULL);
	
	s_factoryArray->remove(factory);
	
	if( s_factoryArray->isEmpty() ) {
		delete s_factoryArray;
		s_factoryArray = NULL;
	}
}
