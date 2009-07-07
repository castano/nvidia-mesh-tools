// Copyright NVIDIA Corporation 2007 -- Denis Kovacs <den.kovacs@gmail.com>

#ifndef NV_MESH_CUBICGREGORYMESH_H
#define NV_MESH_CUBICGREGORYMESH_H

#include <nvcore/Containers.h>
#include <nvmath/Vector.h>
#include <nvmesh/nvmesh.h>
#include <nvmesh/halfedge/HalfEdgeMesh.h>

namespace nv
{
  /// Piecewise cubic Gregory surface.
  class CubicGregoryMesh
  {
  public:
    struct Patch;
    CubicGregoryMesh(uint patchCount)
    {
      m_patchArray.resize(patchCount);
      vVertsTexObj = iVertsTexObj = fWeightsTexObj = 0xffffffff;
    }
    uint patchCount() const { return m_patchArray.count(); }
    const Patch &patchAt(uint i) const { return m_patchArray[i]; }
    Patch &patchAt(uint i) { return m_patchArray[i]; }
    Patch &topologyPatchAt(uint i) { return m_patchArray[m_PatchesSortedByTopology[i]]; }
    const Array<Patch> &patchArray() const { return m_patchArray; }

    unsigned vVertsTexObj, iVertsTexObj, fWeightsTexObj;

  public:
    void FreeStencils();
  public:
    Array<Patch> m_patchArray;
  public:
    Array<int> m_PatchesSortedByTopology;
    int nTopologies;
  public:
    unsigned nVertsPerTopology;
  };

  struct Stencil
  {
    Stencil()
    {
      memset(weights, 0, sizeof(weights));
    }
    int operator <(const Stencil &s) const
    {
      if (fSum < s.fSum)
        return -1;
      if (fSum > s.fSum)
        return 1;
      for (int iw = 0; iw < 20; ++iw)
      {
        if (weights[iw] < s.weights[iw])
          return -1;
        if (weights[iw] > s.weights[iw])
          return 1;
      }
      return 0; // they are equal
    }
    bool AddWeight(int i, float weight)
    {
      nvCheck(i >= 0 && i <= 19);
      weights[i] += weight;
	  return weights[i] == weight;
    }
    void UpdateSum()
    {
      fSum = weights[0];
      for (int iw = 1; iw < 20; ++iw)
      {
        fSum += weights[iw];
      }
    }
    float weights[20];
    double fSum; // sum of all weights (used to sort stencils)
    Vector3 pos; // vertex the stencil is for
    unsigned indByPos; // index of stencils if stencils are sorted by pos
  };

  class StencilsArray
  {
  public:
    void Clear()
	{
		for (unsigned is = 0; is < stencils.count(); ++is)
		{
			delete stencils[is];
			stencils[is] = NULL;
		}
		stencils.clear();
	}
    Stencil &operator[](const Vector3 &pos)
    {
      for (unsigned i = 0; i < stencils.count(); ++i)
      {
        if (stencils[i]->pos == pos)
          return *stencils[i];
      }
      stencils.push_back(new Stencil);
      stencils[stencils.count() - 1]->pos = pos;
      return *stencils[stencils.count() - 1];
    }
    int operator <(const StencilsArray &s) const
    {
      if (stencils.count() < s.stencils.count())
        return -1;
      if (stencils.count() > s.stencils.count())
        return 1;
      for (unsigned is = 0; is < stencils.count(); ++is)
      {
        int r = *stencils[is] < *s.stencils[is];
        if (r != 0)
        {
          return r;
        }
      }
      return 0;
    }
    void Divide(int i, float val)
    {
      nvCheck(i >= 0 && i <= 19);
      for (unsigned j = 0; j < stencils.count(); ++j)
      {
        stencils[j]->weights[i] /= val;
      }
    }
    void CopyWeights(int src, int dst)
    {
      nvCheck(src >= 0 && src <= 19);
      nvCheck(dst >= 0 && dst <= 19);
      for (unsigned i = 0; i < stencils.count(); ++i)
      {
        stencils[i]->weights[dst] = stencils[i]->weights[src];
      }
    }
    void AddWeights(int src, int dst)
    {
      nvCheck(src >= 0 && src <= 19);
      nvCheck(dst >= 0 && dst <= 19);
      for (unsigned i = 0; i < stencils.count(); ++i)
      {
        stencils[i]->weights[dst] += stencils[i]->weights[src];
      }
    }
    void CopyWeights(int src, StencilsArray &d, int dst);
    Vector3 EvalPoint(int ip, bool bQuad);
    void Sort();

  public:
    int iTopology; // topology index
	bool bRegular; // true if patch is regular
    Array<Stencil *> stencils;
    // when doing EvalPoint(), verts must be sumed up in order prescribed by this array (otherwise
    // we will get non water-tight results due to fp-precision issues)
    struct lkjadfg
    {
      int iv, is[2];
    };
    Array<lkjadfg> iVerts;
  };

  /// Bicubic patch or quartic triangle patch produced by Loop's approximation.
  struct CubicGregoryMesh::Patch
  {
    // Vertex indices:

    // quad:
    //   8    13     14      9
    //  12   0\1     2/3    15
    //  19   4/5     6\7    16
    //  11    18     17     10
    // triangle:
    //               15   16
    //                  6
    //              
    //             9   0/1   10
    //                          
    //        14   5/4     3\2   11
    //               
    //  20  8      13       12      7  17
    //     19                      18
    Vector3 positionArray[36];
	Vector3 vW[4];
	Vector3 vSeed[4];
	int iValence[4];
	int iWT0[4];
	int iWT1[4];
	unsigned iOwnershipMask;

	int iEdge0(int index);
	int iEdge1(int index);
	int iCorner(int index);

    // tangents for triangle:
    //             0
    //
    //        8         1
    //
    //    7        9         2
    //               
    // 6      5         4        3
    Vector3 along[9];
    Vector3 across[12];
    Vector3 tangents[16];
    Vector3 bitangents[16];

    // Texture coordinates.
	Vector2 texcoords[16]; // 4 * 4 for quads,  4 * 3 for triangles

	// Tangent basis.
	Vector3 normals[4];
	//Vector3 texTangents[4];
	//Vector3 texBitangents[4];

    float valenceArray[4];
    StencilsArray stencils;

    bool isQuad; // true for quads, false for triangles

    Vector3 evalP(float u, float v) const;
    Vector3 evalTV(float u, float v) const;
    //Vector3 evalTV0(double u, double v) const;
    Vector3 evalTU(float u, float v) const;
    //Vector3 evalTU0(double u, double v) const;
    void RewritePositionArrayFromStencils();
    void CheckPositionsWithStencil();
	void SetupIndices(int idst_corner, Patch &src_patch, int isrc_corner);

    int face_id; // index of face for the patch in the base mesh
  };

} // nv namespace

#endif // NV_MESH_CUBICGREGORYMESH_H
