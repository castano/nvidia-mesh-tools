// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_MESHTRANSFORM_H
#define NV_MESH_MESHTRANSFORM_H

#include <nvmesh/nvmesh.h>

namespace nv
{
	class BaseMesh;
	namespace HalfEdge { class Mesh; }
	class Vector3;
	class Matrix;
	class Box;

	namespace MeshTransform
	{
		NVMESH_API void transform(HalfEdge::Mesh * mesh, const Matrix & matrix);
		NVMESH_API float fitBox(HalfEdge::Mesh * mesh, const Box & box);
		NVMESH_API void translate(HalfEdge::Mesh * mesh, const Vector3 & v);
		//NVMESH_API void flipAxis(HalfEdge::Mesh * mesh, int a, int b, int c);
		
		NVMESH_API void transform(BaseMesh * mesh, const Matrix & matrix);
		NVMESH_API float fitBox(BaseMesh * mesh, const Box & box);
		NVMESH_API void translate(BaseMesh * mesh, const Vector3 & v);
	}

} // nv namespace


#endif // NV_MESH_MESHTRANSFORM_H
