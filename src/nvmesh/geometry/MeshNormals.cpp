// This code is in the public domain -- castanyo@yahoo.es

#include <nvmesh/TriMesh.h>
#include <nvmesh/QuadTriMesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>
#include <nvmesh/halfedge/HalfEdgeVertex.h>
#include <nvmesh/halfedge/HalfEdgeFace.h>
#include "MeshNormals.h"

using namespace nv;

namespace
{
/*
	/// A vertex normal to compute object normals with smooth groups.
	struct VertexNormal
	{
		/// Ctor.
		VertexNormal() : m_normal(zero), m_averageNormal(zero), m_group(NIL)
		{
		}
	
		/// Add normal.
		void addNormal(const Vector3 & n, uint g)
		{
			// This *must* be true.
			// You have to weld the vertices after computing the smooth groups and before computing the normals!
			nvCheck(m_group == NIL || m_group == g);
			
			m_normal += n;
			m_averageNormal += n;
			m_group = g;
		}
	
		/// Add two vertex normals.
		void addNormal(const VertexNormal & n)
		{
			//if( ptr->group & group ) {	// Max style smooth groups.
			if( n.m_group == m_group ) {
				m_averageNormal += n.m_normal;
			}
		}
	
		/// Normalize the normals.
		void normalize()
		{
			m_averageNormal = ::normalize(m_averageNormal);
		}
	
		/// Get the normal.
		const Vector3 & normal() const
		{
			return m_averageNormal;
		}
	
	public:
	
		Vector3 m_normal;
		Vector3 m_averageNormal;
		uint m_group;
	
	};
*/
} // namespace


void /*nv::MeshNormals::*/computeNormals(HalfEdge::Mesh * mesh, uint flags, MeshSmoothGroup * meshSmoothGroup/*= NULL*/)
{
	nvCheck(mesh != NULL);
	nvCheck(meshSmoothGroup == NULL);
	
	const bool weightFaceArea = (flags & WeightFaceArea) != 0;
	const bool weightFaceAngle = (flags & WeightFaceAngle) != 0;

	const uint vertexCount = mesh->vertexCount();
	for(uint v = 0; v < vertexCount; v++ ) {
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);
		
		Vector3 normal(zero);
		const Vector3 & v0 = vertex->pos();
		
		for(HalfEdge::Vertex::EdgeIterator it(vertex->edges()); !it.isDone(); it.advance()) {
			HalfEdge::Face * face = it.current()->face();
			
			if (face != NULL) {
				
				const Vector3 & v1 = it.current()->to()->pos();
				const Vector3 & v2 = it.current()->prev()->from()->pos(); 
				
				Vector3 e0 = v1 - v0;
				Vector3 e1 = v2 - v0;
				
				// @@ Make sure the edges are not degenerate.
				
				// Compute face normal. @@ Assumes the face is planar.
				normal = normalize(cross(e0, e1));
				
				// Compute face area.
				if (weightFaceArea) {
					normal *= face->area();
				}
				
				// @@ Compute face angle.
				if (weightFaceAngle) {
					//normal = ;
				}
			}
		}
		
		vertex->setNor(normal);
	}
}


bool nv::MeshNormals::hasNormals(const BaseMesh * mesh)
{
	nvCheck(mesh != NULL);

	const uint vertexCount = mesh->vertexCount();
	for(uint v = 0; v < vertexCount; v++ )
	{
		const BaseMesh::Vertex & vertex = mesh->vertexAt(v);

		if (equal(vertex.nor, Vector3(zero)))
		{
			return false;
		}
	}

	return true;
}

// Fast normal computation without smooth groups nor creases. Does not require colocal information.
void nv::MeshNormals::computeNormals(TriMesh * mesh, uint flags)
{
	nvCheck(mesh != NULL);

	nvDebug("--- Computing normals.\n");

	const bool weightFaceArea = (flags & WeightFaceArea) != 0;
	const bool weightFaceAngle = (flags & WeightFaceAngle) != 0;

	const uint vertexCount = mesh->vertexCount();
	HashMap<Vector3, Vector3> vertexNormalMap(vertexCount);
	
	const uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++ ) {

		const TriMesh::Face & face = mesh->faceAt(f);
		
		Vector3 p[3];
		for(int i = 0; i < 3; i++) {
			p[i] = mesh->vertexAt(face.v[i]).pos;
		}
		
		const Vector3 A = p[1] - p[0];
		const Vector3 B = p[2] - p[0];
		const Vector3 C = p[2] - p[1];

		Vector3 faceNormal = cross(A, B);

		if (!weightFaceArea) {
			faceNormal = normalizeSafe(faceNormal, Vector3(zero), 0.0f);
		}
		
		float a[3] = {1.0f, 1.0f, 1.0f};
		
		if( weightFaceAngle ) {
			float l2 = length(A);
			float l1 = length(B);
			float l0 = length(C);

			// Heron formulas. They do not work very well for sliver triangles.
			a[0] = acosf( (l2*l2 + l1*l1 - l0*l0) / (2*l2*l1) );
			a[1] = acosf( (l0*l0 + l2*l2 - l1*l1) / (2*l0*l2) );
			//a[2] = acosf( (l1*l1 + l0*l0 - l2*l2) / (2*l1*l0) );
			a[2] = PI - a[0] - a[1];
		}
		
		for(int i = 0; i < 3; i++) {
			Vector3 n(zero);
			vertexNormalMap.get(p[i], &n);
			vertexNormalMap.set(p[i], n + faceNormal * a[i]);
		}
	}

	// Add normals.
	for(uint v = 0; v < vertexCount; v++ ) {

		TriMesh::Vertex & vertex = mesh->vertexAt(v);
		
		Vector3 normal;
		if (vertexNormalMap.get(vertex.pos, &normal))
		{
			vertex.nor = normalizeSafe(normal, Vector3(zero), 0.0f);
		}
	}
}


// Fast normal computation without smooth groups nor creases. 
// fast hack, needs review
void nv::MeshNormals::computeNormals(QuadTriMesh * mesh, uint flags)
{
	nvCheck(mesh != NULL);

	nvDebug("--- Computing normals.\n");

	const bool weightFaceArea = (flags & WeightFaceArea) != 0;
	const bool weightFaceAngle = (flags & WeightFaceAngle) != 0;

	const uint vertexCount = mesh->vertexCount();
	HashMap<Vector3, Vector3> vertexNormalMap(vertexCount);
	
	const uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++ )
	{
		const QuadTriMesh::Face & face = mesh->faceAt(f);
		
		uint vertexCount = 3;
		if (mesh->isQuadFace(f)) vertexCount++ ;

		Vector3 p[4];
		
		for(uint i = 0; i < vertexCount; i++) {
			p[i] = mesh->vertexAt(face.v[i]).pos;
		}

		for (uint i = 0; i < vertexCount; i++)
		{
			const uint prev = (i + (vertexCount - 1)) % vertexCount;
			const uint next = (i + 1) % vertexCount;

		//	nvDebug("%u %u %u\n", i, prev, next);

			Vector3 A = p[next] - p[i];
			Vector3 B = p[prev] - p[i];

			Vector3 faceNormal = cross(A, B);

			if (!weightFaceArea) {
				faceNormal = normalizeSafe(faceNormal, Vector3(zero), 0.0f);
			}

			if( weightFaceAngle ) {
				Vector3 nA = normalizeSafe(A, Vector3(zero), 0.0f);
				Vector3 nB = normalizeSafe(B, Vector3(zero), 0.0f);
				faceNormal *= acosf(dot(nA, nB));
			}

			Vector3 n(zero);
			vertexNormalMap.get(p[i], &n);
			vertexNormalMap.set(p[i], n + faceNormal);
		}
	}

	// Add normals.
	for(uint v = 0; v < vertexCount; v++ ) {

		QuadTriMesh::Vertex & vertex = mesh->vertexAt(v);
		
		Vector3 normal;
		if (vertexNormalMap.get(vertex.pos, &normal))
		{
			vertex.nor = normalizeSafe(normal, Vector3(zero), 0.0f);
		}
	}
}



#if 0
// @@ I think this code is very crappy, it needs to be reviewed.
// TriMesh needs to have VertexData with colocal vertex information.
void nv::MeshNormals::computeNormals(TriMesh * mesh, bool weightFaceArea, bool weightFaceAngle, MeshSmoothGroup * meshSmoothGroup/*= NULL*/)
{
	nvCheck(mesh != NULL);
	nvCheck(meshSmoothGroup == NULL);	// @@ smooth groups not supported yet.

	const uint vertexCount = mesh->vertexCount();
	
	HashMap<Vector3, Vector3> vertexNormal;
	
	Array<VertexNormal> vertexNormalArray;
	vertexNormalArray.resize(vertexCount);

	const uint faceCount = mesh->faceCount();
	for(uint f = 0; f < faceCount; f++ ) {

		// Get face group.
		uint faceGroup = 0;
		//if( meshSmoothGroup != NULL && meshSmoothGroup->HasNormalFaceGroups() ) {
		//	faceGroup = mesh_smooth_group->GetNormalFaceGroup(f);
		//}

		const TriMesh::Face & face = mesh->faceAt(f);
		const uint v0 = face.v[0];
		const uint v1 = face.v[1];
		const uint v2 = face.v[2];
		
		Vector3 A, B, C;
		A = mesh->faceVertex(f, 1).pos - mesh->faceVertex(f, 0).pos;
		B = mesh->faceVertex(f, 2).pos - mesh->faceVertex(f, 0).pos;
		C = mesh->faceVertex(f, 2).pos - mesh->faceVertex(f, 1).pos;

		Vector3 faceNormal = cross(A, B);

		if (!weightFaceArea) {
			faceNormal = normalizeSafe(faceNormal, Vector3(zero), 0.0f);
		}
		
		if( weightFaceAngle ) {
			float l2 = length(A);
			float l1 = length(B);
			float l0 = length(C);

			// Heron formulas. They do not work very well for sliver triangles.
			float a0 = acos( (l2*l2 + l1*l1 - l0*l0) / (2*l2*l1) );
			float a1 = acos( (l0*l0 + l2*l2 - l1*l1) / (2*l0*l2) );
			//a2 = acos( (l1*l1 + l0*l0 - l2*l2) / (2*l1*l0) );
			float a2 = PI - a0 - a1;

			vertexNormalArray[v0].addNormal(faceNormal * a0, faceGroup);
			vertexNormalArray[v1].addNormal(faceNormal * a1, faceGroup);
			vertexNormalArray[v2].addNormal(faceNormal * a2, faceGroup);
		}
		else {
			vertexNormalArray[v0].addNormal(faceNormal, faceGroup);
			vertexNormalArray[v1].addNormal(faceNormal, faceGroup);
			vertexNormalArray[v2].addNormal(faceNormal, faceGroup);
		}
	}

	// Add normals.
	for(uint v = 0; v < vertexCount; v++ ) {

		//const PiMeshAdjacency::Vertex & vertex = mesh_adjacency->GetVertex( v );

		// Add the normals of the proximals.
		//foreach( p, vertex.proximal ) {
		//	vnormal_array[v].AddNormal( vnormal_array[vertex.proximal[p]] );
		//}

		mesh->vertexAt(v).nor = normalizeSafe(vertexNormalArray[v].normal(), Vector3(zero), 0.0f);
	}
}
#endif

#if 0

void nv::MeshNormals::computeCreasedNormals(HalfEdge::Mesh * mesh, bool weightFaceArea, bool weightFaceAngle, MeshSmoothGroup * meshSmoothGroup/*= NULL*/)
{
	nvCheck(mesh != NULL);
	nvCheck(meshSmoothGroup == NULL);
	
	/* This should be quite easy:
	- Traverse every vertex.
	- Traverse each face around the vertex.
	- Measure angle between adjacent faces around vertex.
	- If below crease angle, then add to vertex.
	- If over crease angle, start a new vertex, assign new smooth group.
	- Vertex smooth groups?
	*/
}

#endif


/**
 * Dave Moore wrote:
 *   Just to throw another technique out there, if you know that your mesh
 * is a manifold, or a manifold with a boundary, you can get the vertex
 * normals by computing the tangent plane to the Loop subdivision limit
 * surface, which is well defined at all verts.  This is actually much
 * faster than it sounds, and has some nice properties:
 * 
 * - You only need to be able to find the verts in your 1-ring which
 *   means you only need a very simple connectivity structure
 * - No intermediate (i.e, "face") normal calculations, the whole thing
 *   is dependant on vertex positions.
 * - Hordes of brainy Ph.D's have verified that the solution exists  :)
 * 
 *   Being practical though, I couldn't tell the difference (visually)
 * between normals calculated with Loop, Face area weighted, or Face angle
 * weighted.  Especially if you have to deal with non-manifold meshes, stay
 * away from this one...:)  You can find the formula for this in "Real Time
 * Rendering, 2nd" which I unfortunately don't have handy.
 * 
 *   Quick follow up: If p_0..p_(n-1) are the (ordered) neighbor vertices
 * for vertex p, the two tangents and the limit normal of the Loop surface
 * are:
 * 
 * t_u = sum[j=0..n-1]( cos(2 PI j / n) * p_j )
 * t_v = sum[j=0..n-1]( sin(2 PI j / n) * p_j )
 * n   = cross(t_u, t_v);
 * 
 *   Sorry to clutter the list with stuff that can be easily looked up
 * (RTR2, p536), but that's a fun set of formulas.  Looking at them, it's
 * kind of hard to believe that they're actually doing anything interesting
 * at all.  The cos/sin values are something I'd expect to find in a DFT,
 * rather than a surface normal calculation, in particular.
 *
 * ---
 * Ignacio Castano:
 * I added support for proper handling of boundaries.
 *
 * The stencils for the limit tangents are described in the Appendix A of:
 * http://www.mrl.nyu.edu/publications/piecewise-smooth/
 *
 */
void nv::MeshNormals::computeLoopNormals(HalfEdge::Mesh * mesh)
{
	nvCheck(mesh != NULL);

	// Traverse vertices.	
	const uint vertexCount = mesh->vertexCount();
	for (uint v = 0; v < vertexCount; v++ )
	{
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);
		
		const uint valence = vertex->valence();

		Vector3 tu(zero);
		Vector3 tv(zero);

		if (vertex->isBoundary())
		{
			// Compute limit vertex tangents.
			if (valence == 2)
			{
				tu = -1.0f * vertex->pos();

				const HalfEdge::Edge * edge0 = vertex->edge();
				nvCheck(edge0->face() == NULL);
				nvCheck(edge0->to() != vertex);

				tu += 1.0f / 2.0f * edge0->to()->pos();
				tv += 1.0f / 2.0f * edge0->to()->pos();

				const HalfEdge::Edge * edge1 = vertex->edge()->prev();
				nvCheck(edge1->face() == NULL);
				nvCheck(edge1->from() != vertex);

				tu += 1.0f / 2.0f * edge1->from()->pos();
				tv -= 1.0f / 2.0f * edge1->from()->pos();
			}
			else
			{
				const int k = valence - 1;
				const float omega = PI / k;
				const float s = sinf(omega);
				const float c = cosf(omega);
				float foo = acosf(c - 1);
				float a = ((1 + c) / 4.0f) / (3.0f / 2.0f - 3.0f / 4.0f * c);
				float b = (2.0f / 3.0f - a) / cosf(k * foo / 2.0f);
				float sigma1 = s / (1 - c);
				float sigma3 = cosf(k * foo / 2.0f) * s;

				tv = -2.0f / k * ((2.0f / 3.0f - a) * sigma1  - b * sigma3) * vertex->pos();
				
				const HalfEdge::Edge * edge0 = vertex->edge();
				nvCheck(edge0->face() == NULL);
				nvCheck(edge0->to() != vertex);

				tu += 1.0f / 2.0f * edge0->to()->pos();
				tv -= 2.0f / k * ((a / 2.0f + 1.0f / 6.0f) * sigma1 + b * sigma3 / 2.0f) * edge0->to()->pos();

				const HalfEdge::Edge * edge1 = vertex->edge()->prev();
				nvCheck(edge1->face() == NULL);
				nvCheck(edge1->from() != vertex);

				tu -= 1.0f / 2.0f * edge1->from()->pos();
				tv -= 2.0f / k * ((a / 2.0f + 1.0f / 6.0f) * sigma1 + b * sigma3 / 2.0f) * edge1->from()->pos();

				int j = 0;
				for(HalfEdge::Vertex::EdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), j++) {
					const Vector3 & p = it.current()->to()->pos();
					
					// This is zero for j = 0 and j = k.
					tv += 2.0f / k * sinf(omega * j) * p;
				}
			}
		}
		else
		{
			int j = 0;
			for(HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), j++)
			{
				const Vector3 & p = it.current()->to()->pos();
				tu += p * cosf(2 * PI * j / valence);
				tv += p * sinf(2 * PI * j / valence);
			}
		}

		Vector3 normal = normalizeSafe(cross(tv, tu), Vector3(zero), 0.0f);

		for(HalfEdge::Vertex::VertexIterator it(vertex->colocals()); !it.isDone(); it.advance())
		{
			HalfEdge::Vertex * colocal = it.current();
			colocal->setNor(normal);
		}
	}
}

namespace 
{
	inline static float alphaWeightU(int i, int n)
	{
		const float cp = cosf(PI/n);
		const float denom = (n * sqrtf(4 + cp*cp));
		return (1.0f / n + cp / denom) * cosf((2.0f * PI * i) / n);  
	}

	inline static float betaWeightU(int i, int n)
	{
		const float cp = cosf(PI/n);
		const float denom = (n * sqrtf(4 + cp*cp));
		return (1.0f / denom) * cosf((2.0f * PI * i + PI) / n);  
	}

	inline static float alphaWeightV(int i, int n)
	{
		const float cp = cosf(PI/n);
		const float denom = (n * sqrtf(4 + cp*cp));
		return (1.0f / n + cp / denom) * sinf((2.0f * PI * i) / n);  
	}

	inline static float betaWeightV(int i, int n)
	{
		const float cp = cosf(PI/n);
		const float denom = (n * sqrtf(4 + cp*cp));
		return (1.0f / denom) * sinf((2.0f * PI * i + PI) / n);  
	}

} // namespace


Vector3 nv::MeshNormals::computeCatmullClarkNormal(const HalfEdge::Vertex * vertex)
{
	nvDebugCheck(vertex != NULL);

	const uint valence = vertex->valence();

	if (vertex->isBoundary())
	{
		if (valence == 2)
		{
			Vector3 tu = 6.0f * vertex->pos();
			Vector3 tv(zero);

			const HalfEdge::Edge * edge0 = vertex->edge();
			nvCheck(edge0->face() == NULL);
			nvCheck(edge0->to() != vertex);

			tu -= 3.0f * edge0->to()->pos();
			tv -= 1.0f * edge0->to()->pos();

			const HalfEdge::Edge * edge1 = vertex->edge()->prev();
			nvCheck(edge1->face() == NULL);
			nvCheck(edge1->from() != vertex);

			tu -= 3.0f * edge1->from()->pos();
			tv += 1.0f * edge1->from()->pos();

			return normalizeSafe(cross(tu, tv), Vector3(zero), 0.0f);
		}
		else
		{
			const int k = valence - 1;
			const float omega = PI / k;
			const float s = sinf(omega);
			const float c = cosf(omega);
			const float R = (c + 1) / (k * s * (3 + c));

			Vector3 tu = (-4 * s) / (3 * k + c) * vertex->pos();
		//	Vector3 tu = 4.0f * R * (c - 1.0f) * vertex->pos();
			Vector3 tv(zero);
			
			const HalfEdge::Edge * edge0 = vertex->edge();
			nvCheck(edge0->face() == NULL);
			nvCheck(edge0->to() != vertex);

			tu -= ((1 + 2 * c) * sqrtf(1 + c)) / ((3 * k + c) * sqrtf(1 - c)) * edge0->to()->pos();
		//	tu -= R * (1.0f + 2.0f * c) * edge0->to()->pos();
			tv += 1.0f / 2.0f * edge0->to()->pos();

			const HalfEdge::Edge * edge1 = vertex->edge()->prev();
			nvCheck(edge1->face() == NULL);
			nvCheck(edge1->from() != vertex);

			tu -= ((1 + 2 * c) * sqrtf(1 + c)) / ((3 * k + c) * sqrtf(1 - c)) * edge1->from()->pos();
		//	tu -= R * (1.0f + 2.0f * c) * edge1->from()->pos();
			tv -= 1.0f / 2.0f * edge1->from()->pos();

			int j = 0;
			for (HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), j++)
			{
				if (j == 0) continue;
				if (j == k) break;

				const HalfEdge::Edge * edge = it.current();

				Vector3 p = edge->to()->pos();
				Vector3 q = edge->pair()->prev()->from()->pos();

				tu += 4 * sinf(j * omega) / (3 * k + c) * p;
				tu += 4 * (sinf(j * omega) * sinf((j + 1) * omega)) / (3 * k + c) * q;
			//	tu += 4.0f * sinf(j * omega) / ((3.0f + c) * k) * p;
			//	tu += 1.0f * (sinf(j * omega) + sinf((j + 1.0f) * omega)) / ((3.0f + c) * k) * q;
			}
			
			return normalizeSafe(cross(tu, tv), Vector3(zero), 0.0f);
		}
	}
	else
	{
		Vector3 tu(zero);
		Vector3 tv(zero);
		
		int i = 0;
		for(HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), i++)
		{
			const HalfEdge::Edge * edge = it.current();
			
			Vector3 p0 = edge->to()->pos();
			Vector3 p1 = edge->pair()->prev()->from()->pos();

			tu += p0 * alphaWeightV(i, valence);
			tu += p1 * betaWeightV(i, valence);
			
			tv += p0 * alphaWeightU(i, valence);
			tv += p1 * betaWeightU(i, valence);
		}
		
		return normalizeSafe(cross(tu, tv), Vector3(zero), 0.0f);
	}
}

// Compute tangent in the direction of the given edge.
Vector3 nv::MeshNormals::computeCatmullClarkTangent(const HalfEdge::Edge * firstEdge)
{
	nvDebugCheck(firstEdge != NULL);

	const HalfEdge::Vertex * vertex = firstEdge->vertex();
	const uint valence = vertex->valence();
	
	// Edge direction for debug purposes; tangent is roughly in the same direction.
	const Vector3 dir = normalize(firstEdge->to()->pos() - firstEdge->from()->pos());

	Vector3 tangent(zero);

	if (vertex->isBoundary())
	{
		if (valence == 2)
		{
			tangent = firstEdge->to()->pos() - firstEdge->from()->pos();
		}
		else
		{
			const int k = valence - 1;
			const float omega = PI / k;
			const float s = sinf(omega);
			const float c = cosf(omega);

			Vector3 tu = (-4 * s) / (3 * k + c) * vertex->pos();
			Vector3 tv(zero);

			const HalfEdge::Edge * edge0 = vertex->edge();
			nvDebugCheck(edge0->face() == NULL);
			nvDebugCheck(edge0->to() != vertex);

			tu -= ((1 + 2 * c) * sqrtf(1 + c)) / ((3 * k + c) * sqrtf(1 - c)) * edge0->to()->pos();
			tv += 1.0f / 2.0f * edge0->to()->pos();

			const HalfEdge::Edge * edge1 = vertex->edge()->prev();
			nvDebugCheck(edge1->face() == NULL);
			nvDebugCheck(edge1->from() != vertex);

			tu -= ((1 + 2 * c) * sqrtf(1 + c)) / ((3 * k + c) * sqrtf(1 - c)) * edge1->from()->pos();
			tv -= 1.0f / 2.0f * edge1->from()->pos();

			int i = -1;
			int j = 0;
			for (HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges()); !it.isDone(); it.advance(), j++)
			{
				const HalfEdge::Edge * edge = it.current();

				if (firstEdge == edge) i = j;
				if (j == k) break;

				Vector3 p = edge->to()->pos();
				Vector3 q = edge->pair()->prev()->from()->pos();

				if (j != 0) tu += 4 * sinf(j * omega) / (3 * k + c) * p;
				tu += (sinf(j * omega) + sinf((j + 1) * omega)) / (3 * k + c) * q;
			}

			tangent = sinf(i * omega) * tu + cosf(i * omega) * tv;
		}
	}
	else
	{
		int i = 0;
		for (HalfEdge::Vertex::ConstEdgeIterator it(vertex->edges(firstEdge)); !it.isDone(); it.advance(), i++)
		{
			const HalfEdge::Edge * edge = it.current();

			Vector3 p0 = edge->to()->pos();
			Vector3 p1 = edge->pair()->prev()->from()->pos();

			float alpha = alphaWeightU(i, valence);
			float beta = betaWeightU(i, valence);

			tangent += p0 * alpha;
			tangent += p1 * beta;
		}
	}

//	tangent = normalizeSafe(tangent, Vector3(zero), 0.0f);

	// Make sure tangent is aligned properly.
	//nvDebugCheck(dot(tangent, dir) > 0);

	return tangent;
}



// Same as above, but using the Catmull Clark limit normal. Assumes a quad mesh.
void nv::MeshNormals::computeCatmullClarkNormals(HalfEdge::Mesh * mesh)
{
	nvCheck(mesh != NULL);

	// Traverse vertices.	
	const uint vertexCount = mesh->vertexCount();
	for (uint v = 0; v < vertexCount; v++ )
	{
		HalfEdge::Vertex * vertex = mesh->vertexAt(v);
		nvDebugCheck(vertex != NULL);

		Vector3 normal = computeCatmullClarkNormal(vertex);

		// Make sure all colocals get exactly the same normal.
		for (HalfEdge::Vertex::VertexIterator it(vertex->colocals()); !it.isDone(); it.advance())
		{
			HalfEdge::Vertex * colocal = it.current();
			colocal->setNor(normal);
		}
	}
}






