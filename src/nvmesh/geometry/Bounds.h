// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_MESHBOUNDS_H
#define NV_MESH_MESHBOUNDS_H

//#include <nvmath/Sphere.h>
#include <nvmath/Box.h>

#include <nvmesh/nvmesh.h>

namespace nv
{
	class BaseMesh;
	namespace HalfEdge { class Mesh; }

	/// Bounding volumes computation.
	namespace MeshBounds 
	{
	//	Sphere sphere(const Mesh & mesh);
		Box box(const BaseMesh * mesh);
		Box box(const HalfEdge::Mesh * mesh);
	}

} // nv namespace

#endif // NV_MESH_MESHBOUNDS_H
