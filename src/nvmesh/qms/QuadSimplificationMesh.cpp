
#include <qms/QuadSimplificationMesh.h>

using namespace nv;
using namespace QuadMeshSimplification;

namespace nv
{
	namespace QuadMeshSimplification
	{
		enum
		{
			State_ChordCollapse,
			State_QuadCollapse,
			State_DoubletCollapse,
		};

		struct SimplificationState
		{
			HalfEdge::Mesh * mesh;
			float alpha_v;
			float alpha_q;
			float alpha_d;
			uint maxQuad;
			float maxQEM;
			int state;
		};
	}
}


void QuadMeshSimplification::simplify(HalfEdge::Mesh * mesh, uint maxQuad, float maxQEM, float v/*= 0.9f*/, float q/*= 0.05f*/, float d/*= 0.05f*/)
{
	SimplificationState * state = simplifyStart(mesh, maxQuad, maxQEM, v, q, d);
	while (simplifyStep(state)) {
		// @@ Print debug messages.
	}
	simplifyEnd(state);
}

void QuadMeshSimplification::smooth(HalfEdge::Mesh * mesh)
{
	// @@ Not implemented.
}

void QuadMeshSimplification::polyChordCollapse(HalfEdge::Edge * edge)
{
}

void QuadMeshSimplification::quadCollapse(HalfEdge::Face * face)
{
}

void QuadMeshSimplification::doubletCollapse(HalfEdge::Vertex * vertex)
{

}


/*
bool QuadMeshSimplification::isDoublet(const HalfEdge::Face * face0, const HalfEdge::Face * face1)
{
	for (HalfEdge::Face::ConstEdgeIterator it(face0->edges()); !it.isDone(); it.advance())
	{
		it
	}
}
*/


QuadMeshSimplification::SimplificationState * QuadMeshSimplification::simplifyStart(HalfEdge::Mesh * mesh, uint maxQuad, float maxQEM, float v/*= 0.9f*/, float q/*= 0.05f*/, float d/*= 0.05f*/)
{
	SimplificationState * state = new SimplificationState;

	state->mesh = mesh;
	state->alpha_v = v;
	state->alpha_q = q;
	state->alpha_d = d;
	state->maxQuad = maxQuad;
	state->maxQEM = maxQEM;

	return state;
}

bool QuadMeshSimplification::simplifyStep(SimplificationState * state)
{
	bool done = false;

	while (!done)
	{
		if (state->state == State_ChordCollapse)
		{
		}
		else if (state->state == State_ChordCollapse)
		{
		}
		else if (state->state == State_DoubletCollapse)
		{
		}
	}

	return done;
}

void QuadMeshSimplification::simplifyEnd(SimplificationState * state)
{
	delete state;
}

