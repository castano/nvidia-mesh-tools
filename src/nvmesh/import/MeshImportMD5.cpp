// This code is in the public domain -- castanyo@yahoo.es

#include <nvcore/StrLib.h>
#include <nvcore/Tokenizer.h>

#include <nvmath/Vector.h>

#include <nvmesh/animation/Bone.h>

#include "MeshImportMD5.h"

using namespace nv;

namespace 
{
	
	bool parseVector(Tokenizer & parser, float * out, int num)
	{
		if (!parser.nextToken() || parser.token() != "(") return false;
		
		for (int i = 0; i < num; i++)
		{
			if (!parser.nextToken()) return false;
			out[i] = parser.token().toFloat();
		}
		
		if (!parser.nextToken() || parser.token() != ")") return false;
		
		return true;
	}

	bool parseVector3(Tokenizer & parser, Vector3 * vec)
	{
		float v[3];
		if (!parseVector(parser, v, 3)) return false;
		*vec = Vector3(v[0], v[1], v[2]);
		return true;
	}

	bool parseVector2(Tokenizer & parser, Vector2 * vec)
	{
		float v[2];
		if (!parseVector(parser, v, 2)) return false;
		*vec = Vector2(v[0], v[1]);
		return true;
	}

} // namespace



/// Import a MD5 mesh file.
bool MeshImportMD5::import(Stream * stream)
{
	Tokenizer parser(stream);
	parser.setDelimiters("(){}");
	
	// Reset mesh.
	m_builder.reset();

	const char * name = "???";

	// Get version.
	if (!parser.nextLine() || parser.token() != "MD5Version")
	{
		nvDebug("*** Error, invalid MD5 file: '%s'\n", name);
		return false;
	}
	
	if (!parser.nextToken() || parser.token() != "10")
	{
		nvDebug("*** Error, '%s' has an invalid version: %s != 10\n", name, parser.token().toString().str());
		return false;
	}

	while (parser.nextLine())
	{
		if (parser.token() == "commandline")
		{
			parser.nextToken();
			nvDebug("---   Command line: \"%s\"\n", parser.token().toString().str());
		}
		else if (parser.token() == "numJoints")
		{
			parser.nextToken();
			int n = parser.token().toInt();
			nvDebug("---   %d joints.\n", n);
			
			m_skeleton.boneArray().reserve(n);
		}
		else if (parser.token() == "numMeshes")
		{
			parser.nextToken();
			int n = parser.token().toInt();
			nvDebug( "---   %d meshes.\n", n );
		}
		else if (parser.token() == "joints")
		{
			if (!parseJoints(parser))
			{
				nvDebug( "*** Error, parsing joints in file '%s'.\n", name );
				return false;
			}
		}
		else if (parser.token() == "mesh")
		{
			if (!parseMesh(parser))
			{
				nvDebug( "*** Error, parsing mesh in file '%s'.\n", name );
				return false;
			}
		}
		else if (parser.token() == "//")
		{
			// skip line.
		}
		else
		{
			nvDebug( "*** Error, unexpected token '%s' in file '%s'.\n", parser.token().toString().str(), name );
			return false;
		}
	}

	// Init the mesh in the default pose.
	m_skeleton.setDefaultPose();

	const uint vertexCount = m_skeleton.vertexArray().count();
	m_builder.hintPositionCount(vertexCount);
	
	for (uint i = 0; i < vertexCount; i++)
	{
		m_builder.addPosition(m_skeleton.vertexAt(i).pos);
	}

	m_builder.done();

	return true;
}


// Parse joints.
bool MeshImportMD5::parseJoints(Tokenizer & parser)
{
	parser.nextToken(true);
	if (parser.token() != "{")
	{
		return false;
	}

	while (parser.nextLine())
	{
		if (parser.token() == "}")
		{
			// done;
			return true;
		}
		else if (parser.token() == "//")
		{
			// skip line.
		}
		else
		{
			String name = parser.token().toString();

			if (!parser.nextToken()) break;
			int parent = parser.token().toInt();

			Vector3 offset;
			if (!parseVector3(parser, &offset)) break;

			Vector3 tmp;
			if (!parseVector3(parser, &tmp)) break;
			
			Quaternion rotation;
			float r = 1.0f - length_squared(tmp);
			if( r < 0.0f ) {
				rotation = Quaternion(tmp.x(), tmp.y(), tmp.z(), 0.0f);
			}
			else {
				rotation = Quaternion(tmp.x(), tmp.y(), tmp.z(), sqrtf(r));
			}
			// No need to normalize!
			//rotation = normalize(rotation);

			m_skeleton.addAbsoluteBone(name, parent, rotation, offset);
		}
	}

	// @@ Premature end of file!
	return false;
}



/// Parse the mesh.
bool MeshImportMD5::parseMesh(Tokenizer & parser)
{
	// Count vertices and links before adding this mesh.
	uint base_vertex = m_skeleton.vertexArray().count();
	uint base_link = m_skeleton.linkArray().count();

	parser.nextToken(true);
	if (parser.token() != "{")
	{
		return false;
	}

	while (parser.nextLine())
	{
		if (parser.token() == "shader")
		{
			parser.nextToken();
			m_builder.beginMaterial(parser.token().toString());
		}
		else if (parser.token() == "numverts")
		{
			parser.nextToken();
			int n = parser.token().toInt();
			m_builder.hintVertexCount(n);
			m_builder.hintTexCoordCount(n);
			m_skeleton.vertexArray().reserve(m_skeleton.vertexArray().count() + n);
		}
		else if (parser.token() == "vert")
		{
			// skip index.
			parser.nextToken();
			int idx = parser.token().toInt();
			
			Vector2 texcoord;
			parseVector2(parser, &texcoord);
			m_builder.addTexCoord(texcoord);

			Skeleton::Vertex vertex;

			parser.nextToken();
			vertex.link_first = base_link + parser.token().toInt();
			
			parser.nextToken();
			vertex.link_num = parser.token().toInt();

			m_skeleton.vertexArray().append(vertex);
		}
		else if (parser.token() == "numtris")
		{
			parser.nextToken();
			int n = parser.token().toInt();
			m_builder.hintTriangleCount(n);
		}
		else if (parser.token() == "tri")
		{
			// skip index.
			parser.nextToken();
			int idx = parser.token().toInt();

			m_builder.beginPolygon();

			parser.nextToken();
			idx = base_vertex + parser.token().toInt();
			m_builder.addVertex(idx, NIL, idx);

			parser.nextToken();
			idx = base_vertex + parser.token().toInt();
			m_builder.addVertex(idx, NIL, idx);

			parser.nextToken();
			idx = base_vertex + parser.token().toInt();
			m_builder.addVertex(idx, NIL, idx);

			m_builder.endPolygon();
		}
		else if (parser.token() == "numweights")
		{
			parser.nextToken();
			int n = parser.token().toInt();
			m_skeleton.linkArray().reserve(m_skeleton.linkArray().count() + n);
		}
		else if (parser.token() == "weight")
		{
			// skip index.
			parser.nextToken();
			int idx = parser.token().toInt();
			
			Skeleton::Link link;

			parser.nextToken();
			link.bone = parser.token().toInt();

			parser.nextToken();
			link.weight = parser.token().toFloat();

			parseVector3(parser, &link.offset);

			m_skeleton.linkArray().append(link);
		}
		else if (parser.token() == "}")
		{
			return true;
		}
		else if (parser.token() == "//")
		{
			// skip line.
		}
		else
		{
			nvDebug("*** Error, unexpected token '%s'.\n", parser.token().toString().str());
			break;
		}
	}

	// @@ Premature end of file!
	return false;
}


