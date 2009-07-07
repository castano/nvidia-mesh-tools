// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include <nvmesh/subdiv/RemapFaces.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>

using namespace nv;

namespace
{

	struct Orientation
	{
		Orientation() : u(0), v(0), uSet(false), vSet(false) {}

		uint8 u : 1;
		uint8 v : 1;
		uint8 uSet : 1;
		uint8 vSet : 1;

		void setU(bool b)
		{
			nvDebugCheck(uSet == false);
			u = b;
			uSet = true;
		}
		
		void setV(bool b)
		{
			nvDebugCheck(vSet == false);
			v = b;
			vSet = true;
		}
	};


	static HalfEdge::Edge * faceEdge(HalfEdge::Face * face, int idx)
	{
		int i = 0;
		HalfEdge::Face::EdgeIterator it(face->edges());
		
		while (i != idx)
		{
			nvCheck(!it.isDone());
			i++;
			it.advance();
		}
		
		return it.current();
	}

	static int faceEdgeIndex(const HalfEdge::Face * face, const HalfEdge::Edge * edge)
	{
		int i = 0;
		for (HalfEdge::Face::ConstEdgeIterator it(face->edges()); !it.isDone(); it.advance(), i++)
		{
			if (edge == it.current()) break;
		}
		return i;
	}

	static void propagate(Array<Orientation> & orientationArray, const HalfEdge::Edge * edge, bool dir)
	{
		nvDebugCheck(edge != NULL);
		
		while(true)
		{
			const HalfEdge::Edge * faceEdge = edge->pair();
			const HalfEdge::Face * face = faceEdge->face();
			
			if (face == NULL || face->edgeCount() != 4)
			{
				// Stop at mesh boundaries or non-quad faces.
				break;
			}
			
			int edgeIndex = faceEdgeIndex(face, faceEdge);
			
			Orientation & faceOrientation = orientationArray[face->id()];
			
			if (edgeIndex == 1 || edgeIndex == 3)
			{
				if (faceOrientation.uSet) {
					nvDebugCheck(faceOrientation.u == ((edgeIndex == 1) ^ dir));
					break;
				}
				faceOrientation.setU((edgeIndex == 1) ^ dir);
			}
			else // if (edgeIndex == 0 || edgeIndex == 2)
			{
				if (faceOrientation.vSet) {
					nvDebugCheck(faceOrientation.v == ((edgeIndex == 0) ^ dir));
					break;
				}
				faceOrientation.setV((edgeIndex == 0) ^ dir);
			}
			
			edge = faceEdge->next()->next();
		}
	}

} // namespace


// Remap faces in order to avoid parametric discontinuities.
// Only works with quad meshes.
void nv::RemapFaces::consistentPatchOrientation(HalfEdge::Mesh * mesh)
{
	nvDebugCheck(mesh != NULL);
	
	const uint faceCount = mesh->faceCount();
	
	// Assign patch edges. Make sure that orientation of patches is consistent to avoid parametric discontinuities.
	Array<Orientation> orientationArray;
	orientationArray.resize(faceCount);

	for(uint f = 0; f < faceCount; f++)
	{
		const HalfEdge::Face * face = mesh->faceAt(f);
		
		if (face->edgeCount() != 4)
		{
			// Skip non-quad faces.
			continue;
		}
		
		const HalfEdge::Edge * edges[4];
		edges[0] = face->edge();
		edges[1] = edges[0]->next();
		edges[2] = edges[1]->next();
		edges[3] = edges[2]->next();
		nvDebugCheck(edges[3]->next() == edges[0]);
		
		if (!orientationArray[f].uSet)
		{
			orientationArray[f].setU(false);
			
			// Propagate to neighbors. That is through edges 1 and 3.
			propagate(orientationArray, edges[1], false);
			propagate(orientationArray, edges[3], true);
		}
		if (!orientationArray[f].vSet)
		{
			orientationArray[f].setV(false);
			
			propagate(orientationArray, edges[0], false);
			propagate(orientationArray, edges[2], true);
		}
	}

	Array<const HalfEdge::Edge *> edgeArray;
	edgeArray.resize(faceCount, NULL);

    for(uint f = 0; f < faceCount; f++)
	{
		HalfEdge::Face * face = mesh->faceAt(f);
		
		// Determine edge from orientation flags.
		const Orientation & o = orientationArray[f];
		static const int index[4] = {0, 3, 1, 2};
		const int idx = index[(o.v << 1) + o.u];
		
		face->setEdge(faceEdge(face, idx));
	}
}


static uint topologyId(const HalfEdge::Face * face, const HalfEdge::Edge * edge)
{
	nvDebugCheck(face != NULL);
	nvDebugCheck(edge != NULL);
	
	const uint edgeCount = face->edgeCount();
	
	// Only valid for quads.
	nvCheck(edgeCount == 4 || edgeCount == 3);
	
	uint id = 0;
	
	for (HalfEdge::Face::ConstEdgeIterator it(face->edges(edge)); !it.isDone(); it.advance())
	{
		const HalfEdge::Edge * edge = it.current();
		const HalfEdge::Vertex * vertex = edge->from();
		
		uint valence = vertex->valence();
		nvCheck(valence < 128);
		
		if (vertex->isBoundary())
		{
			//valence += 128;
			valence = ~valence + 1;
		}
		
		id = (id << 8) + (valence & 0xFF);
	}
	
	return id;
}

uint nv::RemapFaces::topologyId(const HalfEdge::Face * face)
{
	return ::topologyId(face, face->edge());
}



void nv::RemapFaces::minimizeTopologyCount(HalfEdge::Mesh * mesh)
{
	const uint faceCount = mesh->faceCount();
	
    for(uint f = 0; f < faceCount; f++)
	{
		HalfEdge::Face * face = mesh->faceAt(f);
		
		HalfEdge::Edge * bestEdge = NULL;
		uint bestId = uint(~0);
		
		const uint edgeCount = face->edgeCount();
		if (edgeCount != 3 && edgeCount != 4) 
		{
			// Only handle triangles and polygons.
			continue;
		}

		for (HalfEdge::Face::EdgeIterator it(face->edges()); !it.isDone(); it.advance())
		{
			HalfEdge::Edge * edge = it.current();
			
			uint id = ::topologyId(face, edge);
			
			if (id < bestId)
			{
				bestId = id;
				bestEdge = edge;
			}
		}
		
		// Set edge with unique topology id.
		face->setEdge(bestEdge);
	}
	
	// @@ Sort faces according to topology.
	
}


