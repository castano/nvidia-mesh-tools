// Copyright NVIDIA Corporation 2008 -- Edgar Velazquez-Armendariz <edgarv@nvidia.com>

#include "TiledTask.h"

#include <nvcore/Algorithms.h>

using namespace nv;
using namespace IlmThread;


TiledTask::TiledTask(const Tile & tile, TiledTaskGenerator & generator) :
Task( & generator.m_taskGroup ),
m_tile(tile), m_generator(generator) {}

TiledTask::~TiledTask() {
	m_generator.m_queueLock.acquire();
	nvDebugCheck( (uint)m_generator.m_finishQueue.size() < m_generator.numTiles() );
	m_generator.m_finishQueue.push(m_tile);
	m_generator.m_queueCountSem.post();
	m_generator.m_queueLock.release();
}


TiledTaskGenerator::TiledTaskGenerator(uint width, uint height, uint8 tileSize) : 
  m_width(width), m_height(height),
  m_tileSize(tileSize),
  m_tiles(NULL),
  m_queueLock(m_queueMutex, false),
  m_queueCountSem(0)
{
	nvCheck(tileSize > 0);

	const uint & w = width;
	const uint & h = height;

	m_numTilesX = w/tileSize + (w % tileSize > 0 ? 1 : 0);
	m_numTilesY = h/tileSize + (h % tileSize > 0 ? 1 : 0);
	nvCheck(m_numTilesX > 0 && m_numTilesY > 0);

	m_tiles = new Array2D<Tile>(m_numTilesY /* rows */, m_numTilesX /* columns */);
	

	// @@ Generates the tiles
	for (uint i = 0; i < m_numTilesY; ++i) {
		const uint yLo = i*m_tileSize;
		const uint yHi = min((yLo+tileSize), h);
		for (uint j = 0; j < m_numTilesX; ++j) {
			const uint xLo = j*m_tileSize;
			const uint xHi = min((xLo+tileSize), w);
			nvCheck(xHi > xLo && yHi > yLo && yLo >= 0 && xLo >= 0);

			// Create tile with its (i,j) coordinates and the range
			// of pixels it takes
			Tile & tile = m_tiles->get(i,j);
			tile.m = i;
			tile.n = j;
			tile.xLo = xLo;
			tile.xHi = xHi;
			tile.yLo = yLo;
			tile.yHi = yHi;
		}
	}
}

TiledTaskGenerator::~TiledTaskGenerator() {
	if (m_tiles != NULL) {
		delete m_tiles;
	}
}


void TiledTaskGenerator::runTasks(IlmThread::ThreadPool &pool) {

	//ScanlineIterator it(m_tiles->rows(), m_tiles->columns());
	HilbertIterator it(m_tiles->rows(), m_tiles->columns());
	while(it.hasNext()) {
		pool.addTask( newTask(m_tiles->get(it.next())) );
	}
}

Tile TiledTaskGenerator::take() {
	
	// Waits until there is something
	m_queueCountSem.wait();

	// And retrieves the corresponding tile
	m_queueLock.acquire();
	nvDebugCheck( m_finishQueue.size() > 0 );
	Tile tile = m_finishQueue.front();
	m_finishQueue.pop();
	m_queueLock.release();

	return tile;
}
