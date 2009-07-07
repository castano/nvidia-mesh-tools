// Copyright NVIDIA Corporation 2007 -- Ignacio Castano <icastano@nvidia.com>

#include "MayaMeshBuilder.h"

// Required by maya!
#include <iostream>

#include <maya/MSelectionList.h>
#include <maya/MGlobal.h>
#include <maya/MItSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MFnMesh.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MPointArray.h>

using namespace nv;


MayaMeshBuilderOptions::MayaMeshBuilderOptions() :
	addTexcoords(false),
	addNormals(false)
{
}

MayaMeshBuilderOptions::MayaMeshBuilderOptions(const MayaMeshBuilderOptions & options)
{
	this->addTexcoords = options.addTexcoords;
	this->addNormals = options.addNormals;
}

const MayaMeshBuilderOptions & MayaMeshBuilderOptions::operator=(const MayaMeshBuilderOptions & options)
{
	this->addTexcoords = options.addTexcoords;
	this->addNormals = options.addNormals;
	return *this;
}



MayaMeshBuilder::MayaMeshBuilder(const MayaMeshBuilderOptions & options) :
	m_options(options)
{
}


bool MayaMeshBuilder::addScene()
{
	// @@ TODO!
	return false;
}

bool MayaMeshBuilder::addSelection()
{
	MSelectionList list; 
	MGlobal::getActiveSelectionList(list); 

	MItSelectionList iter(list, MFn::kMesh);
    if (iter.isDone())
	{
        return false;
    }

	bool succeed = false;

	for (; !iter.isDone(); iter.next()) 
	{
		MDagPath node;
		MObject component;
		iter.getDagPath(node, component);

		if (addNode(node))
		{
			succeed = true;
		}
	}

	return succeed;
}

bool MayaMeshBuilder::addNode(const MDagPath & dagPath)
{
	MStatus status;
	MFnMesh meshFn(dagPath, &status);

	MItMeshPolygon polyIt(dagPath, MObject::kNullObj, &status);
	if (MS::kSuccess != status) return false; // @@ Not a polygonal mesh?

	// Add positions.
	MPointArray positionArray;
	status = meshFn.getPoints(positionArray, MSpace::kObject);

	const uint positionCount = positionArray.length();
	m_builder.hintPositionCount(positionCount);
	
	for (uint i = 0; i < positionCount; i++)
	{
		MPoint point = positionArray[i];
		m_builder.addPosition(Vector3(point.x, point.y, point.z));
	}

	if (m_options.addNormals)
	{
		// Add normals.
		MFloatVectorArray normalArray;
		status = meshFn.getNormals(normalArray, MSpace::kObject);

		const uint normalCount = normalArray.length();
		m_builder.hintNormalCount(normalCount);

		for (uint i = 0; i < normalCount; i++)
		{
			MFloatVector normal = normalArray[i];
			m_builder.addNormal(Vector3(normal.x, normal.y, normal.z));
		}
	}

	if (m_options.addTexcoords)
	{
		// Add texcoords. @@ Only the default UV set.
		MFloatArray uArray;
		MFloatArray vArray;

		status = meshFn.getUVs(uArray, vArray);

		const uint uCount = uArray.length();
		const uint vCount = vArray.length();
		nvCheck(uCount == vCount);
		m_builder.hintTexCoordCount(uCount);

		for (uint i = 0; i < uCount; ++i)
		{
			m_builder.addTexCoord(Vector2(uArray[i], vArray[i]));
		}
	}

	
	// Add polygons.
	for (; !polyIt.isDone(); polyIt.next())
	{
		const uint vertexCount = polyIt.polygonVertexCount();

		if (vertexCount >= 3)
		{
			m_builder.beginPolygon();

			for (uint i = 0; i < vertexCount; i++)
			{
				uint normalIndex = NIL;
				int texCoordIndex = NIL;

				if (m_options.addNormals) normalIndex = polyIt.normalIndex((int)i);
				if (m_options.addTexcoords) polyIt.getUVIndex((int)i, /*ref*/texCoordIndex);

				m_builder.addVertex(polyIt.vertexIndex((int)i), normalIndex, (uint)texCoordIndex);
			}

			m_builder.endPolygon();
		}
	}

	return true;
}

void MayaMeshBuilder::done()
{
	m_builder.done();
}

void MayaMeshBuilder::reset()
{
	m_builder.reset();
}

TriMesh * MayaMeshBuilder::buildTriMesh() const
{
	return m_builder.buildTriMesh();
}

QuadTriMesh * MayaMeshBuilder::buildQuadTriMesh() const
{
	return m_builder.buildQuadTriMesh();
}

/*PolyMesh * MayaMeshBuilder::buildPolyMesh() const
{
	return m_builder.buildPolyMesh();
}*/

HalfEdge::Mesh * MayaMeshBuilder::buildHalfEdgeMesh() const
{
	return m_builder.buildHalfEdgeMesh();
}


uint MayaMeshBuilder::vertexCount() const
{
	return m_builder.vertexCount();
}

uint MayaMeshBuilder::positionIndex(uint vertex) const
{
	return m_builder.positionIndex(vertex);
}
uint MayaMeshBuilder::normalIndex(uint vertex) const
{
	return m_builder.normalIndex(vertex);
}
uint MayaMeshBuilder::texcoordIndex(uint vertex) const
{
	return m_builder.texcoordIndex(vertex);
}
