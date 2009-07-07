// This code is in the public domain -- castanyo@yahoo.es

#include <nvmath/TypeSerialization.h>

#include "BaseMesh.h"


namespace nv
{
	static Stream & operator<< (Stream & s, BaseMesh::Vertex & vertex)
	{
		return s << vertex.id << vertex.pos << vertex.nor << vertex.tex;
	}

	Stream & operator<< (Stream & s, BaseMesh & mesh)
	{
		return s << mesh.m_vertexArray;
	}
}
