// This code is in the public domain -- Ignacio Castaño <castanyo@yahoo.es>

#ifndef NV_MESH_RAYTRACING_H
#define NV_MESH_RAYTRACING_H

#include <nvcore/Ptr.h>
#include <nvmath/Vector.h>
#include <nvmesh/nvmesh.h>


namespace nv
{
	struct Ray
	{
		uint id;
		Vector3 origin;
		Vector3 idir;
		Vector3 dir;
		float maxt;
	};
	
	struct Hit
	{
		float t;
		uint face;
		float u, v;
	};
	
} // nv namespace


#endif // NV_MESH_RAYTRACING_H

