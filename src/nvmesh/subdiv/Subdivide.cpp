// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include "Subdivide.h"

#include <nvcore/Ptr.h>

#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdge.h>

#include <nvimage/FloatImage.h>

using namespace nv;

HalfEdge::Mesh * nv::Subdivide::doCatmullClarkSplit(const HalfEdge::Mesh * mesh)
{
	nvDebug("--- Applying Catmull Clark split.\n");
	
	HalfEdge::Mesh * result = new HalfEdge::Mesh();
	
	// Add face vertices.
	const uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++) {
		const HalfEdge::Face * face = mesh->faceAt(f);
		
		// Compute face centroid.
		Vector3 pos(zero);
		Vector2 tex(zero);
		
		int num = 0;
		for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {
			const HalfEdge::Vertex * vertex = it.current()->from();
			pos += vertex->pos();
			tex += vertex->tex();
			num++;
		}
		pos *= 1.0f / num;
		tex *= 1.0f / num;
		
		HalfEdge::Vertex * vertex = result->addVertex(pos);
		nvDebugCheck(vertex != NULL);
		vertex->setTex(tex);
	}
	nvCheck(result->vertexCount() == faceCount);

	// Add edge vertices. Takes texture seams into account.
	uint seamEdgeCount = 0;
	const uint edgeCount = mesh->edgeCount();
	HashMap<const HalfEdge::Edge *, HalfEdge::Vertex *> edgeVertexMap(edgeCount);

	for(uint e = 0; e < edgeCount; e++) {
		const HalfEdge::Edge * edge = mesh->edgeAt(e);
		
		// Compute edge midpoint.
		Vector3 pos(zero);

		uint f0 = edge->face()->id();
		Vector3 c0 = result->vertexAt(f0)->pos();

		if (edge->isBoundary())
		{
			pos = edge->midPoint();
		}
		else 
		{
			uint f1 = edge->pair()->face()->id();
			Vector3 c1 = result->vertexAt(f1)->pos();

			pos = (edge->midPoint() + (c0 + c1) * 0.5) * 0.5;
		}

		HalfEdge::Vertex * vertex = result->addVertex(pos);
		nvDebugCheck(vertex != NULL);
		vertex->setTex(0.5 * (edge->from()->tex() + edge->to()->tex()));

		edgeVertexMap.add(edge, vertex);

		if (!edge->isBoundary()) 
		{
			if (edge->isSeam())
			{
				// Duplicate seam vertices.
				seamEdgeCount++;

				HalfEdge::Vertex * pairVertex = result->addVertex(pos);
				nvDebugCheck(pairVertex != NULL);
				pairVertex->setTex(0.5 * (edge->pair()->from()->tex() + edge->pair()->to()->tex()));

				edgeVertexMap.add(edge->pair(), pairVertex);
			}
			else 
			{
				// Do not add vertices for edges without face.
				edgeVertexMap.add(edge->pair(), vertex);
			}
		}
	}
	nvCheck(result->vertexCount() == faceCount + edgeCount + seamEdgeCount);


	const uint vertexCount = mesh->vertexCount();
	for(uint v = 0; v < vertexCount; v++) {
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		HalfEdge::Vertex * newVertex = NULL;

		if (vertex->isBoundary()) {
			Vector3 pos = 3.0f/4.0f * vertex->pos();

			const HalfEdge::Edge * edge0 = vertex->edge();
			nvCheck(edge0->face() == NULL);
			nvCheck(edge0->to() != vertex);

			pos += 1.0f/8.0f * edge0->to()->pos();

			const HalfEdge::Edge * edge1 = vertex->edge()->prev();
			nvCheck(edge1->face() == NULL);
			nvCheck(edge1->from() != vertex);

			pos += 1.0f/8.0f * edge1->from()->pos();

			newVertex = result->addVertex(pos);
		}
		else {
			const uint valence = vertex->valence();

			if (valence == 0)
			{
				// Copy the isolated vertex as is.
				newVertex = result->addVertex(vertex->pos());
			}
			else if (valence == 2 && false)
			{
				// "Smooth Geometry Images" explains how to deal with valence two vertices.
				// constrain vertex position to the centroid of its four neighbors.
				// The perturbed surface is C1.
				
				// @@ So, what if the neighbor faces are not quads? 
				// -- I average all the neighbor vertices instead. That's essentially the midpoint of the face centroids.
				// That's different to what Losasso et al. do, because they don't consider the current vertex in the 
				// average, only the neighbors. I tried that, but caused undesired rippling, while this method produces
				// smooth results.
				
				Vector3 pos(zero);
				
				const HalfEdge::Edge * edge0 = vertex->edge();
				const HalfEdge::Face * face0 = edge0->face();
				
				uint f0 = face0->id();
				Vector3 c0 = result->vertexAt(f0)->pos();
				
				const HalfEdge::Edge * edge1 = vertex->edge()->pair()->next();
				const HalfEdge::Face * face1 = edge1->face();
				
				uint f1 = face1->id();
				Vector3 c1 = result->vertexAt(f1)->pos();

				newVertex = result->addVertex((c0 + c1) * 0.5);
			}
			else
			{
				// You can only use the stencil weights when all the faces are quads.
				//float gamma = 1.0f / (4.0f * valence);
				//float beta = 3.0f / (2.0f * valence);
				// w0 = (1 - gamma - beta);
				// w1 = gamma / valence;
				// w2 = beta / valence;
		
				// So we use the recursive formula instead.
		
				Vector3 Q(zero);	// Average of face vertices around vertex.
				Vector3 R(zero);	// Average of edge midpoints around vertex.
				Vector3 S = vertex->pos();
		
				// Compute Q and R.
				for(HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance()) {
					const HalfEdge::Edge * edge = it.current();
					const HalfEdge::Face * face = edge->face();
		
					if (face != NULL)
					{
						Q += result->vertexAt(face->id())->pos();
					}
		
					R += edge->midPoint();
				}
		
				Q *= 1.0f / valence;
				R *= 1.0f / valence;
		
				Vector3 pos = (Q + (2.0f * R) + S * (float(valence) - 3)) / float(valence);
				
				newVertex = result->addVertex(pos);
			}
		}
		
		nvDebugCheck(newVertex != NULL);
		newVertex->setTex(vertex->tex());
	}
	nvCheck(result->vertexCount() == faceCount + edgeCount + seamEdgeCount + vertexCount);
	
	result->linkColocals();

	// Add faces.
	for(uint f = 0; f < faceCount; f++) {
		const HalfEdge::Face * face = mesh->faceAt(f);
		nvDebugCheck(face->id() == f);
		
		// Traverse face vertices.
		for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance()) {

			const HalfEdge::Edge * edge0 = it.current();
			const HalfEdge::Edge * edge1 = it.current()->prev();

			nvCheck(edge1->next() == edge0);
			
			nvCheck(edge0->from() == edge1->to());
			
			/*if (edge0->isBoundary() || edge1->isBoundary()) {
				// @@ Skip boundary faces.
				continue;
			}*/
			
			// There's a direct mapping for face and corner vertices.
			// But edge vertices need to be looked up.

			const HalfEdge::Vertex * vertex = edge0->from();
			nvDebugCheck(vertex == edge1->to());

			HalfEdge::Vertex * vertex0 = NULL;
			edgeVertexMap.get(edge0, &vertex0);
			nvDebugCheck(vertex0 != NULL);

			HalfEdge::Vertex * vertex1 = NULL;
			edgeVertexMap.get(edge1, &vertex1);
			nvDebugCheck(vertex1 != NULL);

			uint v0 = f;
			uint v1 = vertex1->id(); //faceCount + edge1->id();
			uint v2 = faceCount + edgeCount + seamEdgeCount + vertex->id();
			uint v3 = vertex0->id(); // faceCount + edge0->id();
			
			result->addFace(v0, v1, v2, v3);
		}
	}
	
	result->linkBoundary();

	return result;
}


void nv::Subdivide::displace(HalfEdge::Mesh * mesh, const FloatImage * img, int channel/*=0*/, Vector2::Arg offset/*=Vector2(zero)*/)
{
	nvCheck(mesh != NULL);
	nvCheck(img != NULL);

	nvDebug("--- Applying displacement map.\n");
	
	const Vector2 texOffset(offset.x()/ img->width(), offset.y()/img->height());
	
	bool isSeamless = true;
	int errorCount = 0;

	const uint vertexCount = mesh->vertexCount();
	for(uint v = 0; v < vertexCount; v++) {
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		// Use volatile to force stack storage and make sure that comparison is exact.

		Vector2 tex = vertex->tex() + texOffset;
		volatile float displacement = img->sampleLinear(tex.x(), tex.y(), channel, FloatImage::WrapMode_Clamp);

		const uint colocalCount = vertex->colocalCount();
		if (colocalCount != 1) {
			float averageDisplacement = 0.0f;

			for(HalfEdge::Vertex::VertexIterator it(vertex->colocals()); !it.isDone(); it.advance()) {
				HalfEdge::Vertex * colocal = it.current();
	
				Vector2 colocalTex = colocal->tex() + texOffset;
				volatile float colocalDisplacement = img->sampleLinear(colocalTex.x(), colocalTex.y(), channel, FloatImage::WrapMode_Clamp);
				averageDisplacement += colocalDisplacement;
				
				if (colocalDisplacement != displacement) {
					isSeamless = false;
					errorCount++;
					//displacement = -1.0f;
				}
			}

			vertex->setPos(vertex->pos() + vertex->nor() * averageDisplacement/float(colocalCount));
		}
		else
		{
			vertex->setPos(vertex->pos() + vertex->nor() * displacement);
		}
	}

	if (!isSeamless) {
		nvDebug("---   Mesh parameterization is *not* seamless. (%d errors)\n", errorCount);
	}
	else {
		nvDebug("---   Mesh parameterization is seamless.\n");
	}
}


HalfEdge::Mesh * nv::Subdivide::doQuadTriangleSplit(const HalfEdge::Mesh * mesh)
{
	// @@ Implement quad/triangle subdivision here.
	// http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/qtEG.pdf
	return NULL;
}


// Current implementation simply removed non-quad faces from the mesh.
HalfEdge::Mesh * nv::Subdivide::reduceToQuadOnly(const HalfEdge::Mesh * mesh)
{
	nvDebug("--- Reduce to Quad Only.\n");
	
	AutoPtr<HalfEdge::Mesh> result(new HalfEdge::Mesh());

	// Add vertices.
	const uint vertexCount = mesh->vertexCount();
	for (uint v = 0; v < vertexCount; v++) {
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);

		HalfEdge::Vertex * newVertex = result->addVertex(vertex->pos());
		nvDebugCheck(newVertex != NULL);

		newVertex->setNor(vertex->nor());
		newVertex->setTex(vertex->tex());
	}
	nvCheck(result->vertexCount() == vertexCount);

	result->linkColocals();

	// Add faces.
	const uint faceCount = mesh->faceCount();
	for (uint f = 0; f < faceCount; f++) {
		const HalfEdge::Face * face = mesh->faceAt(f);
		nvDebugCheck(face->id() == f);

		// Add only quad faces.
		if (face->edgeCount() == 4)
		{
			uint v[4];
			int i = 0;
			for(HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++) {
				v[i] = it.current()->vertex()->id();
			}
			
			result->addFace(v[0], v[1], v[2], v[3]);
		}
	}

	result->linkBoundary();

	// Unlink isolated vertices.
	/*for (uint v = 0; v < vertexCount; v++)
	{
		HalfEdge::Vertex * vertex = result->vertexAt(v);

		if (vertex->edge() == NULL)
		{
			// @@ Ideally we should remove the vertex, but that is a mess right now.
			//vertex->unlink();
			//mesh->removeVertex(vertex);
		}
	}*/

	return result.release();
}


// This method assumes the given mesh is a quad mesh!
void nv::Subdivide::projectToCatmullClarkLimitSurface(HalfEdge::Mesh * mesh)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	Array<Vector3> vertexArray;
	vertexArray.resize(vertexCount);

	for(uint v = 0; v < vertexCount; v++)
	{
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		Vector3 pos = Vector3(zero);

		// Evaluate limit surface pos.
		if (vertex->isBoundary())
		{
			/*if (m_boundaryMode != BoundaryMode_Linear)*/
			{
				// Compute limit vertex position.
				pos = 2.0f / 3.0f * vertex->pos();

				const HalfEdge::Edge * edge0 = vertex->edge();
				nvCheck(edge0->face() == NULL);
				nvCheck(edge0->to() != vertex);

				pos += 1.0f / 6.0f * edge0->to()->pos();

				const HalfEdge::Edge * edge1 = vertex->edge()->prev();
				nvCheck(edge1->face() == NULL);
				nvCheck(edge1->from() != vertex);

				pos += 1.0f / 6.0f * edge1->from()->pos();
			}
		}
		else
		{
			const float valence = float(vertex->valence());

			pos = valence * valence * vertex->pos();

			int i = 0;
			for(HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance(), i++)
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(vertex->pos() == edge->from()->pos());

				pos += 4 * edge->to()->pos();
				pos += edge->next()->to()->pos();
			}

			pos /= valence * valence + valence * 5;
		}

		vertexArray[v] = pos;
	}

	for(uint v = 0; v < vertexCount; v++)
	{
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		vertex->setPos(vertexArray[v]);
	}
}



void nv::Subdivide::projectToLoopLimitSurface(HalfEdge::Mesh * mesh)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	Array<Vector3> vertexArray;
	vertexArray.resize(vertexCount);

	for(uint v = 0; v < vertexCount; v++)
	{
		const HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		Vector3 pos = Vector3(zero);

		// Evaluate limit surface pos.
		if (vertex->isBoundary())
		{
			/*if (m_boundaryMode != BoundaryMode_Linear)*/
			{
				// Compute limit vertex position.
				pos = 2.0f / 3.0f * vertex->pos();

				const HalfEdge::Edge * edge0 = vertex->edge();
				nvCheck(edge0->face() == NULL);
				nvCheck(edge0->to() != vertex);

				pos += 1.0f / 6.0f * edge0->to()->pos();

				const HalfEdge::Edge * edge1 = vertex->edge()->prev();
				nvCheck(edge1->face() == NULL);
				nvCheck(edge1->from() != vertex);

				pos += 1.0f / 6.0f * edge1->from()->pos();
			}
		}
		else
		{
			const uint valence = vertex->valence();

			float beta = 3.0f / (8.0f * valence);
			if (valence == 3) beta = 3.0f / 16.0f;

			pos = vertex->pos();

			int i = 0;
			for(HalfEdge::Vertex::ConstEdgeIterator eit(vertex->edges()); !eit.isDone(); eit.advance(), i++)
			{
				const HalfEdge::Edge * edge = eit.current();
				nvDebugCheck(vertex->pos() == edge->from()->pos());

				pos += (beta * 8.0f / 3.0f) * edge->to()->pos();
				//pos += (1.0f / valence) * edge->to()->pos();
			}

			pos *= 1.0f / (1 + beta * 8 * valence / 3);
			//pos *= 0.5f;
		}

		vertexArray[v] = pos;
	}

	for(uint v = 0; v < vertexCount; v++)
	{
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		vertex->setPos(vertexArray[v]);
	}
}

