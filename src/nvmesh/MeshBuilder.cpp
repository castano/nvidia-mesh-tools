// This code is in the public domain -- castanyo@yahoo.es

#include "MeshBuilder.h"
#include "TriMesh.h"
#include "QuadTriMesh.h"
#include "halfedge/HalfEdgeMesh.h"
#include "halfedge/HalfEdgeVertex.h"

#include "weld/Weld.h"

//#include "MeshGroup.h"
//#include "MeshMaterial.h"

using namespace nv;

/*
By default the mesh builder creates 3 streams (position, normal, texcoord), I'm planning to add support for extra streams as follows:

enum StreamType { StreamType_Float, StreamType_Vector2, StreamType_Vector3, StreamType_Vector4 };

uint addStream(const char *, uint idx, StreamType);

uint addAttribute(float)
uint addAttribute(Vector2)
uint addAttribute(Vector3)
uint addAttribute(Vector4)

struct Vertex
{
	uint pos;
	uint nor;
	uint tex;
	uint * attribs;	// NULL or NIL terminated array?
};

All streams must be added before hand, so that you know the size of the attribs array.

The vertex hash function could be kept as is, but the == operator should be extended to test 
the extra atributes when available.

That might require a custom hash implementation, or an extension of the current one. How to
handle the variable number of attributes in the attribs array?

	bool operator()(const Vertex & a, const Vertex & b) const
	{ 
		if (a.pos != b.pos || a.nor != b.nor || a.tex != b.tex) return false;
		if (a.attribs == NULL && b.attribs == NULL) return true;
		return 0 == memcmp(a.attribs, b.attribs, ???);
	}

We could use a NIL terminated array, or provide custom user data to the equals functor.

vertexMap.setUserData((void *)vertexAttribCount);

bool operator()(const Vertex & a, const Vertex & b, void * userData) const { ... }

*/



namespace 
{

struct Group
{
	Group() : id(NIL) {}
	Group(uint id) : id(id) {}
	
	uint id;
	uint first;
	uint last;
};

struct Material
{
	Material() {}
	Material(const String & str) : name(str) {}
	
	String name;
	uint first;
	uint last;
};

struct Vertex
{
	Vertex() {}
	Vertex(uint p, uint n, uint t) : pos(p), nor(n), tex(t) {}
	
	friend bool operator==(const Vertex & a, const Vertex & b)
	{
		return a.pos == b.pos && a.nor == b.nor && a.tex == b.tex;
	}
	
	uint pos;
	uint nor;
	uint tex;
};

} // namespace


namespace nv
{
	// This is a much better hash than the default and greatly improves performance!
	template <> struct hash<Vertex>
	{
		uint operator()(const Vertex & v) const { return v.pos + v.nor + v.tex; }
	};
}


struct MeshBuilder::PrivateData
{
	PrivateData() : currentGroup(NIL), currentMaterial(NIL), lastIndex(NIL), maxCount(0), faceCount(0) {}
	
	uint pushVertex(uint p, uint n, uint t);
	uint pushVertex(const Vertex & v);
		
	Array<Vector3> posArray;
	Array<Vector3> norArray;
	Array<Vector2> texArray;
	
	Array<Vertex> vertexArray;
	HashMap<Vertex, uint> vertexMap;
	
	Array<Group> groupArray;
	Array<Material> materialArray;
	
	uint currentGroup;	
	uint currentMaterial;
	
	Array<uint> indexArray;
	
	uint lastIndex;
	uint maxCount;
	uint faceCount;
};


uint MeshBuilder::PrivateData::pushVertex(uint p, uint n, uint t)
{
	Vertex v(p, n, t);
	return pushVertex(v);
}

uint MeshBuilder::PrivateData::pushVertex(const Vertex & v)
{
	// Lookup vertex v in map.
	uint idx;
	if (vertexMap.get(v, &idx))
	{
		return idx;
	}
	
	idx = vertexArray.count();
	vertexArray.pushBack(v);
	vertexMap.add(v, idx);

	return idx;
}


MeshBuilder::MeshBuilder() : d(new PrivateData())
{
}

MeshBuilder::~MeshBuilder()
{
	nvDebugCheck(d != NULL);
	delete d;
}

// Builder methods.
uint MeshBuilder::addPosition(const Vector3 & v)
{
	d->posArray.pushBack(v);
	return d->posArray.count() - 1;
}

uint MeshBuilder::addNormal(const Vector3 & v)
{
	d->norArray.pushBack(v);
	return d->norArray.count() - 1;
}

uint MeshBuilder::addTexCoord(const Vector2 & v)
{
	d->texArray.pushBack(v);
	return d->texArray.count() - 1;
}

void MeshBuilder::beginGroup(uint id)
{
	if (d->currentGroup != NIL) {
		endGroup();
	}

	d->currentGroup = d->groupArray.count(); 
	d->groupArray.pushBack(Group(id));
	d->groupArray[d->currentGroup].first = d->indexArray.count();
}

void MeshBuilder::endGroup()
{
	d->groupArray[d->currentGroup].last = d->indexArray.count();
	d->currentGroup = NIL;
}

uint MeshBuilder::beginMaterial(const String & name)
{
	if (d->currentMaterial != NIL) {
		endMaterial();
	}

	d->currentMaterial = d->groupArray.count();
	d->materialArray.pushBack(Material(name));
	d->materialArray[d->currentMaterial].first = d->indexArray.count();

	return d->currentMaterial;
}

void MeshBuilder::beginMaterial(uint id)
{
	if (d->currentMaterial != NIL) {
		endMaterial();
	}

	nvCheck(id < d->materialArray.count());
	
	d->currentMaterial = id;
	d->materialArray.pushBack(d->materialArray[id]);
	d->materialArray[d->currentMaterial].first = d->indexArray.count();
}

void MeshBuilder::endMaterial()
{
	d->materialArray[d->currentMaterial].last = d->indexArray.count();
	d->currentMaterial = NIL;
}

void MeshBuilder::beginPolygon()
{
	d->lastIndex = d->indexArray.count();
	d->indexArray.pushBack(0);
}

uint MeshBuilder::addVertex(uint p, uint n /*= NIL*/, uint t /*= NIL*/)
{
	uint idx = d->pushVertex(p, n, t);
	d->indexArray.pushBack(idx);
	d->indexArray[d->lastIndex]++;
	return idx;
}

uint MeshBuilder::addVertex(const Vector3 & pos)
{
	uint p = addPosition(pos);
	return addVertex(p, NIL, NIL);
}

uint MeshBuilder::addVertex(const Vector3 & pos, const Vector3 & nor, const Vector2 & tex)
{
	uint p = addPosition(pos);
	uint n = addNormal(nor);
	uint t = addTexCoord(tex);
	return addVertex(p, n, t);
}

void MeshBuilder::endPolygon()
{
	uint count = d->indexArray[d->lastIndex];
	nvCheck(count > 2);
	
	d->maxCount = max(d->maxCount, count);
	d->faceCount++;
}

void MeshBuilder::optimize()
{
	if (d->vertexArray.count() == 0)
	{
		return;
	}

	Array<uint> posxrefs;
	Array<uint> norxrefs;
	Array<uint> texxrefs;
	Weld<Vector3> weldVector3;
	Weld<Vector2> weldVector2;
	
	// Weld vertex attributes.
	if (d->posArray.count()) weldVector3(d->posArray, posxrefs);
	if (d->norArray.count()) weldVector3(d->norArray, norxrefs);
	if (d->texArray.count()) weldVector2(d->texArray, texxrefs);
	
	// Remap vertex indices.
	const uint vertexCount = d->vertexArray.count();
	for (uint v = 0; v < vertexCount; v++)
	{
		Vertex & vertex = d->vertexArray[v];
		if (d->posArray.count()) vertex.pos = posxrefs[vertex.pos];
		if (d->norArray.count()) vertex.nor = norxrefs[vertex.nor];
		if (d->texArray.count()) vertex.tex = texxrefs[vertex.tex];
	}


	Array<uint> xrefs;
	Weld<Vertex> weldVertex;
	
	// Weld vertices.
	weldVertex(d->vertexArray, xrefs);
	
	// Remap face indices.
	const uint indexCount = d->indexArray.count();
	for (uint i = 0; i < indexCount; i++)
	{
		uint count = d->indexArray[i];
		while(count--)
		{
			i++;
			d->indexArray[i] = xrefs[d->indexArray[i]];
		}
	}
	
	// Remap vertex map.
	foreach(i, d->vertexMap)
	{
		d->vertexMap[i].value = xrefs[d->vertexMap[i].value];
	}
}

void MeshBuilder::reset()
{
	nvDebugCheck(d != NULL);
	delete d;
	d = new PrivateData();
}

void MeshBuilder::done()
{
	if (d->currentGroup != NIL) {
		endGroup();
	}
	
	if (d->currentMaterial != NIL) {
		endMaterial();
	}
}

// Hints.
void MeshBuilder::hintTriangleCount(uint count)
{
	d->indexArray.reserve(d->indexArray.count() + count * 4);
}

void MeshBuilder::hintVertexCount(uint count)
{
	d->vertexArray.reserve(d->vertexArray.count() + count);
	d->vertexMap.resize(d->vertexMap.count() + count);
}

void MeshBuilder::hintPositionCount(uint count)
{
	d->posArray.reserve(d->posArray.count() + count);
}

void MeshBuilder::hintNormalCount(uint count)
{
	d->norArray.reserve(d->norArray.count() + count);
}

void MeshBuilder::hintTexCoordCount(uint count)
{
	d->texArray.reserve(d->texArray.count() + count);
}



// Helpers.
void MeshBuilder::addTriangle(uint v0, uint v1, uint v2)
{
	beginPolygon();
	addVertex(v0);
	addVertex(v1);
	addVertex(v2);
	endPolygon();
}

void MeshBuilder::addQuad(uint v0, uint v1, uint v2, uint v3)
{
	beginPolygon();
	addVertex(v0);
	addVertex(v1);
	addVertex(v2);
	addVertex(v3);
	endPolygon();
}


// Get tri mesh.
TriMesh * MeshBuilder::buildTriMesh(uint group /*= NIL*/, uint material /*= NIL*/) const
{
	// @@ group and material are not used yet.
	NV_UNUSED(group);
	NV_UNUSED(material);

	// Count triangles.
	uint triangleCount = 0;
	const uint indexCount = d->indexArray.count();
	for(uint i = 0; i < indexCount; i++)
	{
		int count = d->indexArray[i];
		triangleCount += count-2;
		i += count;
	}

	const uint vertexCount = d->vertexArray.count();
	TriMesh * mesh = new TriMesh(triangleCount, vertexCount);

	// Build faces.
	Array<TriMesh::Face> & faces = mesh->faces();

	for(uint i = 0; i < indexCount; i++)
	{
		int count = d->indexArray[i];

		int v0 = d->indexArray[i+1];
		int v1 = d->indexArray[i+2];

		// @@ Use strip order triangulation.
		for(int t = 0; t < count-2; t++) {
			int v2 = d->indexArray[i+t+3];

			TriMesh::Face face;
			face.id = faces.count();
			face.v[0] = v0;
			face.v[1] = v1;
			face.v[2] = v2;
			faces.append(face);

			v1 = v2;
		}

		i += count;
	}

	// Build vertices.
	Array<BaseMesh::Vertex> & vertices = mesh->vertices();

	for(uint i = 0; i < vertexCount; i++)
	{
		BaseMesh::Vertex vertex;
		vertex.id = i;
		if (d->vertexArray[i].pos != NIL) vertex.pos = d->posArray[d->vertexArray[i].pos];
		if (d->vertexArray[i].nor != NIL) vertex.nor = d->norArray[d->vertexArray[i].nor];
		if (d->vertexArray[i].tex != NIL) vertex.tex = d->texArray[d->vertexArray[i].tex];

		vertices.append(vertex);
	}

	return mesh;
}

// Get quad/tri mesh.
QuadTriMesh * MeshBuilder::buildQuadTriMesh(uint group /*= NIL*/, uint material /*= NIL*/) const
{
	// @@ group and material are not used yet.
	NV_UNUSED(group);
	NV_UNUSED(material);

	// Only triangle/quad meshes supported.
	nvCheck(d->maxCount == 4 || d->maxCount == 3);

	const uint faceCount = d->faceCount;
	const uint vertexCount = d->vertexArray.count();
	QuadTriMesh * mesh = new QuadTriMesh(faceCount, vertexCount);
	
	// Build faces.
	Array<QuadTriMesh::Face> & faces = mesh->faces();

	const uint indexCount = d->indexArray.count();
	for(uint i = 0; i < indexCount; i++) 
	{
		int count = d->indexArray[i];

		QuadTriMesh::Face face;
		face.id = i;

		face.v[0] = d->indexArray[i+1];
		face.v[1] = d->indexArray[i+2];
		face.v[2] = d->indexArray[i+3];

		if (count == 3) {
			face.v[3] = NIL;
		}
		else if (count == 4) {
			face.v[3] = d->indexArray[i+4];
		}
		faces.append(face);

		i += count;
	}

	// Build vertices.
	Array<BaseMesh::Vertex> & vertices = mesh->vertices();

	for(uint i = 0; i < vertexCount; i++)
	{
		BaseMesh::Vertex vertex;
		vertex.id = i;
		if (d->vertexArray[i].pos != NIL) vertex.pos = d->posArray[d->vertexArray[i].pos];
		if (d->vertexArray[i].nor != NIL) vertex.nor = d->norArray[d->vertexArray[i].nor];
		if (d->vertexArray[i].tex != NIL) vertex.tex = d->texArray[d->vertexArray[i].tex];

		vertices.append(vertex);
	}

	return mesh;
}

// Get half edge mesh.
HalfEdge::Mesh * MeshBuilder::buildHalfEdgeMesh(uint group /*= NIL*/, uint material /*= NIL*/) const
{
	NV_UNUSED(group);
	NV_UNUSED(material);

	const uint vertexCount = d->vertexArray.count();
	HalfEdge::Mesh * mesh = new HalfEdge::Mesh();

	for(uint v = 0; v < vertexCount; v++)
	{
		HalfEdge::Vertex * vertex = mesh->addVertex(d->posArray[d->vertexArray[v].pos]);
		if (d->vertexArray[v].nor != NIL) vertex->setNor(d->norArray[d->vertexArray[v].nor]);
		if (d->vertexArray[v].tex != NIL) vertex->setTex(d->texArray[d->vertexArray[v].tex]);
	}

	mesh->linkColocals();
	
	const uint indexCount = d->indexArray.count();
	for(uint i = 0; i < indexCount; i++) 
	{
		int count = d->indexArray[i];

		if (mesh->addFace(d->indexArray, i+1, count) != NULL) {
			// @@ Display warning?
		}

		i += count;
	}

	if (!mesh->linkBoundary())
	{
		return NULL;
	}

	return mesh;
}


// @@ TBD
MeshGroup * MeshBuilder::buildMeshGroup() const
{
	return NULL;
}

// @@ TBD
MeshMaterial * MeshBuilder::buildMeshMaterial() const
{
	return NULL;
}

uint MeshBuilder::vertexCount() const
{
	return d->vertexArray.count();
}

uint MeshBuilder::positionIndex(uint vertex) const
{
	return d->vertexArray[vertex].pos;
}
uint MeshBuilder::normalIndex(uint vertex) const
{
	return d->vertexArray[vertex].nor;
}
uint MeshBuilder::texcoordIndex(uint vertex) const
{
	return d->vertexArray[vertex].tex;
}
