// This code is in the public domain -- castanyo@yahoo.es

#include "MeshImportOFF.h"

#include <nvcore/StrLib.h>
#include <nvcore/Tokenizer.h>

#include <nvmath/Vector.h>


using namespace nv;


/// Import an OFF mesh.
bool MeshImportOFF::import(Stream * stream)
{
	Tokenizer parser(stream);
	
	m_builder.reset();
	
	int vertexCount = 0;
	int faceCount = 0;
	int edgeCount = 0;

	try 
	{
		// Skip comments and empty lines.
		while (parser.nextLine() && parser.token() == "#") {}

		if (parser.token() == "OFF")
		{
			// Skip token.
			parser.nextToken(true);
		}
		
		vertexCount = parser.token().toInt();
		if (parser.nextToken()) faceCount = parser.token().toInt();
		if (parser.nextToken()) edgeCount = parser.token().toInt();

		m_builder.hintPositionCount(vertexCount);
		m_builder.hintVertexCount(vertexCount);
		m_builder.hintTriangleCount(faceCount);
		
		float x, y, z;
		while (vertexCount--)
		{
			// Skip comments and empty lines.
			while (parser.nextLine() && parser.token() == "#") {}

			x = parser.token().toFloat();
			parser.nextToken();
			y = parser.token().toFloat();
			parser.nextToken();
			z = parser.token().toFloat();
			m_builder.addPosition(Vector3(x, y, z));
		}

		while (faceCount--)
		{
			// Skip comments and empty lines.
			while (parser.nextLine() && parser.token() == "#") {}

			// Ignore token.
			//vertexCount = parser.token().toInt();

			m_builder.beginPolygon();
				
			while (parser.nextToken())
			{
				int v = parser.token().toInt();
				
				m_builder.addVertex(v);
			}
			
			m_builder.endPolygon();
		}

		while (edgeCount--)
		{
			// Skip comments and empty lines.
			while (parser.nextLine() && parser.token() == "#") {}

			// Skip edges.
		}
	}
	catch(TokenizerException) {
		m_builder.reset();
		return false;
	}
	
	m_builder.optimize();
	m_builder.done();
	
	return true;
}

// This doesn't work on MSVC!
class MeshImportFactoryOFF : public MeshImportFactory
{
public:
	virtual const char * extension() const {
		return ".off";
	}
	virtual MeshImport * createImporter() const {
		return new MeshImportOFF();
	}
};
NV_REGISTER_MESH_IMPORT_FACTORY(MeshImportFactoryOFF);

