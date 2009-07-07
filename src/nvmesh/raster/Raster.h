// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_RASTER_H
#define NV_MESH_RASTER_H

/** @file Raster.h
 * @brief Rasterization library.
 *
 * This is just a standard scanline rasterizer that I took from one of my old
 * projects. The perspective correction wasn't necessary so I just removed it.
**/

#include <nvmesh/nvmesh.h>
#include <nvmath/Vector.h>

#define RASTER_ANTIALIAS true
#define RASTER_NOAA false

namespace nv
{
	namespace Raster 
	{
		/// A callback to sample the environment.
		typedef void (NV_CDECL * SamplingCallback)(void * param, int x, int y, Vector3::Arg bar, Vector3::Arg dx, Vector3::Arg dy, float coverage);

		// Process the given triangle.
		NVMESH_API bool drawTriangle(bool antialias, Vector2::Arg extents, const Vector2 v[3], SamplingCallback cb, void * param);

		// Process the given quad.
		NVMESH_API bool drawQuad(bool antialias, Vector2::Arg extents, const Vector2 vertex[4], SamplingCallback, void * param);

	}
}


#endif // NV_MESH_RASTER_H
