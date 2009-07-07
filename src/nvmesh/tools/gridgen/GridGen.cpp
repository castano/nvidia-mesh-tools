// Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>

#include <stdio.h>
#include <float.h> // FLT_MAX

#include <nvcore/Timer.h>
#include <nvcore/Ptr.h>
#include <nvcore/Array2D.h>

#include <nvmesh/render/VertexCache.h>
#include <nvmesh/render/NvTriStrip.h>
#include <nvmesh/render/MeshOptimizer.h>

using namespace	nv;


static void gridGen(int x0, int x1, int y0, int y1, int width, int cacheSize, Array<uint> & indices)
{
	int w = x1 - x0;
	int h = y1 - y0;
	
	if (w + 1 < cacheSize)
	{
	//	printf("gen: %dx%d - %dx%d\n", x0, y0, x1, y1);

		// Prefetch only when the strip is small.
		if (w * 2 + 1 > cacheSize)
		{
			for (int x = x0; x < x1; x++) 
			{
				indices.append(x + 0);
				indices.append(x + 0);
				indices.append(x + 1);
			}
		}
		
		for (int y = y0; y < y1; y++)
		{
			for (int x = x0; x < x1; x++) 
			{
				indices.append((width + 1) * (y + 0) + (x + 0));
				indices.append((width + 1) * (y + 1) + (x + 0));
				indices.append((width + 1) * (y + 0) + (x + 1));

				indices.append((width + 1) * (y + 0) + (x + 1));
				indices.append((width + 1) * (y + 1) + (x + 0));
				indices.append((width + 1) * (y + 1) + (x + 1));
			}
		}
	}
	else
	{
		int xm = x0 + cacheSize - 2;
		
		gridGen(x0, xm, y0, y1, width, cacheSize, indices);
		gridGen(xm, x1, y0, y1, width, cacheSize, indices);
	}
}


struct MyMessageHandler : public MessageHandler
{
	void log(const char * str, va_list arg)
	{
		va_list tmp;
		va_copy(tmp, arg);
		vprintf(str, arg);

#if _DEBUG && NV_OS_WIN32
//		static StringBuilder buffer(1024);
//		buffer.format(str, arg);
//		OutputDebugStringA(buffer.str());
#endif

		va_end(tmp);
	}
};


int main(int argc, const char *	argv[])
{
	MyMessageHandler messageHandler;
	debug::setMessageHandler(&messageHandler);

	const char * fileName =	NULL;
    int width = 16;
    int pattern = 0;
	VertexCache::Mode cacheMode = VertexCache::Mode_Batch;
	int cacheSize = 32;
	int targetCacheSize = 0;
	int batchSize = 32;
    bool showHelp = false;
	
	for	(int i = 1;	i <	argc; i++)
	{
		if (strcmp(argv[i], "-w") == 0 || strcmp(argv[i], "--width") == 0)
		{
			width = atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--pattern") == 0)
		{
			pattern = atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--cacheMode") == 0)
		{
			cacheMode = (VertexCache::Mode)atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-n") == 0)
		{
			cacheSize = atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-t") == 0)
		{
			targetCacheSize = atoi(argv[i+1]);
			i++;
		}
		else if (strcmp(argv[i], "-b") == 0)
		{
			batchSize = atoi(argv[i+1]);
			i++;
		}
        else
        {
            showHelp = true;
        }
	}

	if (showHelp)
	{
		printf("Usage: nvgridgen [options]\n");
		printf("options are:\n");
		printf("\t-w, --width: Grid width [default=16]\n");
        printf("\t-p, --pattern: Grid pattern\n");
		printf("\t\t0 = strip [default]\n");
		printf("\t\t1 = prefetch strip\n");
		printf("\t\t2 = morton\n");
		printf("\t\t3 = hilbert\n");
		printf("\t\t4 = nvtristrip\n");
		printf("\t\t5 = Forsyth\n");
		printf("\t\t6 = Tipsy\n");
		printf("\t\t7 = KCacheReorder\n");
		printf("\t\t8 = FermiBatcher\n");
		printf("\t\t9 = Castano\n");
		printf("\t\t10 = Tootle\n");
		printf("\t-c, --cacheMode: Cache mode\n");
		printf("\t\t0 =	batch [default]\n");
		printf("\t\t1 =	FIFO\n");
		printf("\t\t2 =	LRU\n");
		printf("\t-n: Cache size [default=32]\n");
		printf("\t-t: Target cache size\n");
		return 0;
	}

	Array<uint> indices(6 * width * width);
	const uint vertexCount = (width + 1) * (width + 1);


    if (pattern == 1)
    {
		gridGen(0, width, 0, width, width, cacheSize, indices);
    }
	else
	{
		AutoPtr<Iterator2D> it;

		if (pattern == 2)
		{
			it = new MortonIterator (width, width);
		}
		else if (pattern == 3)
		{
			it = new HilbertIterator (width, width);
		}
		else // 1, 4, 5 & 6
		{
			it = new ScanlineIterator (width, width);
		}

		// Generate grid indices.
		while(it->moveNext())
		{
			Index2D idx = it->current();

			indices.append((width + 1) * (idx.y + 0) + (idx.x + 0));
			indices.append((width + 1) * (idx.y + 1) + (idx.x + 0));
			indices.append((width + 1) * (idx.y + 0) + (idx.x + 1));

			indices.append((width + 1) * (idx.y + 0) + (idx.x + 1));	
			indices.append((width + 1) * (idx.y + 1) + (idx.x + 0));
			indices.append((width + 1) * (idx.y + 1) + (idx.x + 1));
		}
	}

    if (pattern == 4) // Optimize with nvtristrip
    {
        TriStrip::SetCacheSize(cacheSize);
        TriStrip::SetListsOnly(true);

        uint count = indices.count();
        Array<uint16> shortIndices;
        shortIndices.resize(count);

        for (uint i = 0; i < count; i++)
        {
            shortIndices[i] = short(indices[i]);
        }

        TriStrip::PrimitiveGroup * groupArray;
        uint16 groupCount;

        TriStrip::GenerateStrips(shortIndices.buffer(), shortIndices.count(), &groupArray, &groupCount);

        nvCheck (groupCount == 1);

        count = groupArray[0].numIndices;
        indices.resize(count);

        for (uint i = 0; i < count; i++)
        {
            indices[i] = groupArray[0].indices[i];
        }

        delete [] groupArray;
    }
	else // Optimize with MeshOptimizer
    {
		uint indexCount = indices.count();
		MeshOptimizer optimizer(vertexCount, indices, 3, cacheSize, cacheMode);

		if (pattern == 5) optimizer.optimize(MeshOptimizer::Method_Forsyth, &indices, NULL);
		else if (pattern == 6) optimizer.optimize(MeshOptimizer::Method_Tipsy, &indices, NULL);
		else if (pattern == 7) optimizer.optimize(MeshOptimizer::Method_KCacheReorder, &indices, NULL);
		else if (pattern == 8) optimizer.optimize(MeshOptimizer::Method_FermiBatcher, &indices, NULL);
		else if (pattern == 9) optimizer.optimize(MeshOptimizer::Method_Castano, &indices, NULL);
		else if (pattern == 10) optimizer.optimize(MeshOptimizer::Method_AMDTootle, &indices, NULL);
	}


	if (targetCacheSize != 0) cacheSize = targetCacheSize;

	const VertexCache::Stats stats  = MeshOptimizer::processIndices(indices, 3, cacheSize, cacheMode, batchSize);
	const uint primitiveCount = 2 * width * width;

	const char * methodName[] = {
		"Scanline",
		"Optimal",
		"Morton",
		"Hilbert",
		"NV TriStrip",
		"Forsyth",
		"Tipsy",
		"K-Cache-Reorder",
		"FermiBatcher",
		"Castano",
		"AMD Tootle",
	};
	nvCheck(pattern < sizeof(methodName)/sizeof(methodName[0]));

	printf("Batcher	result using method: %s\n", methodName[pattern]);
	printf("    Transform count = %u\n", stats.transformCount);
	printf("    Batch count = %u\n", stats.batchCount);
	printf("    ACMR = %.3f\n",	float(stats.transformCount) / primitiveCount);
	printf("    Transform per vertex = %.3f\n", float(stats.transformCount) / vertexCount);
	printf("    Threads per vertex = %.3f\n", float(stats.batchCount * cacheSize) / vertexCount);

	return 0;
}
