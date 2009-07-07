// This code is in the public domain -- castanyo@yahoo.es

#include <nvcore/StrLib.h>
#include <nvcore/Tokenizer.h>

#include <nvmath/Vector.h>

#include "MeshImportOBJ.h"

using namespace nv;


/// Import an OBJ mesh.
bool MeshImportOBJ::import(Stream * stream)
{
	Tokenizer parser(stream);
	parser.setDelimiters("/");
	
	m_builder.reset();
	
	int positionCount = 0;
	int normalCount = 0;
	int texCoordCount = 0;
	
	try {
		float x, y, z;
		while (parser.nextLine()) {
			if (parser.token() == "v") {
				x = y = z = 0.0f;
				if (parser.nextToken()) x = parser.token().toFloat();
				if (parser.nextToken()) y = parser.token().toFloat();
				if (parser.nextToken()) z = parser.token().toFloat();
				m_builder.addPosition(Vector3(x, y, z));
				positionCount++;
			}
			else if (parser.token() == "vn") {
				x = y = z = 0.0f;
				if (parser.nextToken()) x = parser.token().toFloat();
				if (parser.nextToken()) y = parser.token().toFloat();
				if (parser.nextToken()) z = parser.token().toFloat();
				m_builder.addNormal(Vector3(x, y, z));
				normalCount++;
			}
			else if (parser.token() == "vt") {
				x = y = z = 0.0f;
				if (parser.nextToken()) x = parser.token().toFloat();
				if (parser.nextToken()) y = parser.token().toFloat();
				m_builder.addTexCoord(Vector2(x, y));
				texCoordCount++;
			}
			else if (parser.token() == "f") {
				m_builder.beginPolygon();
				
				bool tokenAvailable = false;
				while(tokenAvailable || parser.nextToken()) {
					
					int v = 0;
					int t = 0;
					int n = 0;
					
					v = parser.token().toInt();
					tokenAvailable = parser.nextToken();
					
					if (tokenAvailable && parser.token() == "/") {
						tokenAvailable = parser.nextToken();
						
						if (tokenAvailable && parser.token() != "/") {
							t = parser.token().toInt();
							tokenAvailable = parser.nextToken();
						}
						
						if (tokenAvailable && parser.token() == "/") {
							tokenAvailable = parser.nextToken();
							n = parser.token().toInt();
							tokenAvailable = false;
						}
					}
					
					if(v < 0) v = positionCount - v;
					else if (v == 0) v = NIL;
					else v = v - 1;
					
					if(t < 0) t = positionCount - t;
					else if (t == 0) t = NIL;
					else t = t - 1;
					
					if(n < 0) n = positionCount - n;
					else if (n == 0) n = NIL;
					else n = n - 1;
					
					m_builder.addVertex(v, n, t);
				}
				
				m_builder.endPolygon();
			}
			else if (parser.token() == "g") {
				// @@ Ignore groups.
			}
			else if (parser.token() == "s") {
				if (parser.nextToken()) {
					// @@ I'm not sure that a number is provided.
					uint group = parser.token().toUnsignedInt();
					m_builder.beginGroup(group);
				}
			}
			else if (parser.token() == "usemtl") {
			//	if (parser.nextToken()) {
			//		m_builder.beginMaterial(parser.token().toString());
			//	}
			}
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
class MeshImportFactoryOBJ : public MeshImportFactory
{
public:
	virtual const char * extension() const {
		return ".obj";
	}
	virtual MeshImport * createImporter() const {
		return new MeshImportOBJ();
	}
};
NV_REGISTER_MESH_IMPORT_FACTORY(MeshImportFactoryOBJ);

