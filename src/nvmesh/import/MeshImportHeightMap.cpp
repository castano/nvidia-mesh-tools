// This code is in the public domain -- castanyo@yahoo.es

#include <nvcore/StrLib.h>
#include <nvcore/Tokenizer.h>

#include <nvmath/Vector.h>

#include "MeshImportHeightMap.h"


#define SCALE_X		0.1f
#define SCALE_Y		0.1f
#define SCALE_Z		0.03f


using namespace nv;


/// Import a heightmap mesh.
bool MeshImportHeightMap::import(Stream * stream)
{
	m_builder.reset();
	
	return false;
	
	// @@ Need image file name!
	Image * height_map = ImageIO::load("???", *stream);

	if (height_map == NULL)
	{
		// Image type not recognized.
		return false;
	}
	
	const uint w = height_map->width;
	const uint h = height_map->height;

	for (uint y = 0; y < h; y++)
	{
		for (uint x = 0; x < w; x++)
		{
			float z = height_map->pixel(x, y).a;
			m_builder.addPosition(Vector3(SCALE_X * x, SCALE_Y * y, SCALE_Z * z));
		}
	}

	for (uint y = 0; y < h-1; y++)
	{
		for (uint x = 0; x < w-1; x++)
		{
			m_builder.addTriangle(y*w+x, y*w+x+1, (y+1)*w+x);
			m_builder.addTriangle((y+1)*w+x, y*w+x+1, (y+1)*w+x+1);
		}
	}

	m_builder.done();
	
	return true;
}


