// This code is in the public domain -- castanyo@yahoo.es

#include <nvcore/StrLib.h>
#include <nvcore/Tokenizer.h>

#include <nvmath/Vector.h>

#include "MeshImportM.h"

using namespace nv;


/// Import an M mesh.
bool MeshImportM::import(Stream * stream)
{
	Tokenizer parser(stream);
	
	m_builder.reset();

	try {
		HashMap<uint, uint> indexMap(100);
		
		while (parser.nextLine()) {
			if (parser.token() == "{") {
				// Skip comments
				while (parser.nextToken(true)) {
					if (parser.token() == "}") break;
				}
			}
			else if (parser.token() == "Vertex") {
				if (parser.nextToken()) {
					uint index = parser.token().toUnsignedInt();
					indexMap.add(index, indexMap.count());

					float x, y, z;
					if (parser.nextToken()) x = parser.token().toFloat();
					if (parser.nextToken()) y = parser.token().toFloat();
					if (parser.nextToken()) z = parser.token().toFloat();
					m_builder.addPosition(Vector3(x, y, z));
				}
			}
			else if (parser.token() == "Face") {
				
				m_builder.beginPolygon();
				
				parser.nextToken();	// skip face number

				while(parser.nextToken()) {
					uint idx = NIL;
					if (indexMap.get(parser.token().toUnsignedInt(), &idx))
					{
						m_builder.addVertex(idx);
					}
				}
				
				m_builder.endPolygon();
			}
		}
	}
	catch(TokenizerException) {
		return false;
	}
	
	m_builder.done();
	
	return true;
}


class MeshImportFactoryM : public MeshImportFactory
{
public:
	virtual const char * extension() const {
		return ".m";
	}
	virtual MeshImport * createImporter() const {
		return new MeshImportM();
	}
};
NV_REGISTER_MESH_IMPORT_FACTORY(MeshImportFactoryM);
