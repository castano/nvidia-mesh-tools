// Copyright NVIDIA Corporation 2007 -- Denis Kovacs <den.kovacs@gmail.com>

#include <nvmesh/subdiv/CubicGregoryMesh.h>

using namespace nv;

namespace nv {

	Vector2 operator+(float a, Vector2 b)
	{
		return Vector2(a+b.x(), a+b.y());
	}

	Vector2 operator-(float a, Vector2 b)
	{
		return Vector2(a-b.x(), a-b.y());
	}

}

static const int iQuadIndex[36] =
{ 1,  1,  1,  1,  1,  1, 1,  1,  1,  1,  1, 1, 1,  1,  1,  1,  1,  1,  1, 1, // 20 threads all go from first 5 bits of ivert.y
  1,  1,  1,  1,  1,  2, 2,  2,  2,  2,  2, 3, 3,  3,  3,  3 }; // then we have 6 indices per i32
static const int iQuadShift[36] = 
{ 0,  0,  0,  0,  0,  0, 0,  0,  0,  0,  0, 0, 0,  0,  0,  0,  0,  0,  0, 0,
  5, 10, 15, 20, 25,  0, 5, 10, 15, 20, 25, 0, 5, 10, 15, 20 };

static const int iTrgIndex[27] =
{ 1,  1,  1,  1,  1, 1, 1,  1,  1,  1,  1, 1, 1,  1,  1, // 15 threads all go from first 5 bits of ivert.y
  1,  1,  1,  1,  1, 2, 2,  2,  2,  2,  2, 3 }; // then we have 6 indices per i32
static const int iTrgShift[27] = 
{ 0,  0,  0,  0,  0, 0, 0,  0,  0,  0,  0, 0, 0,  0,  0,
  5, 10, 15, 20, 25, 0, 5, 10, 15, 20, 25, 0 };

static const int iQuadCornerDevotedEdges[] = { 9, 13, 15, 17 };

static const int iTrgCornerDevotedEdges[] = { 10, 12, 14 };

void CubicGregoryMesh::Patch::RewritePositionArrayFromStencils()
{
	for (int ip = 0; ip < 15; ++ip)
	{
		Vector3 p = stencils.EvalPoint(ip, isQuad);
		positionArray[ip] = p;
	}
	if (isQuad)
	{
		for (int ip = 15; ip < 20; ++ip)
		{
			Vector3 p = stencils.EvalPoint(ip, isQuad);
			positionArray[ip] = p;
		}
	}
}
void CubicGregoryMesh::Patch::CheckPositionsWithStencil()
{
	for (int ip = 0; ip < (isQuad ? 20 : 15); ++ip)
	{
		Vector3 p = stencils.EvalPoint(ip, isQuad);
		nvCheck(positionArray[ip] == p);
	}
}
void CubicGregoryMesh::Patch::SetupIndices(int idstCorner, Patch &srcPatch, int isrcCorner)
{
	int isrcSeed = srcPatch.iEdge0(isrcCorner);
	int isrcW = srcPatch.iEdge1(isrcCorner);
	int idstSeed = iEdge0(idstCorner);
	int idstW = iEdge1(idstCorner);
	if (&srcPatch == this) // source and destination pathches are the same
	{
		for (int is = 0; is < (int)stencils.stencils.size(); ++is)
		{
			Stencil *dstStencil = stencils.stencils[is];
			if (iValence[idstCorner] != 4)
			  dstStencil->weights[idstW] = (dstStencil->weights[idstW] - dstStencil->weights[idstSeed] * cos(2 * PI / iValence[idstCorner])) / sin(2 * PI / iValence[idstCorner]);
		}
	}
	else
	{
		int iindex;
		int ishift;
		if (isQuad)
		{
			iindex = (idstCorner < 2) ? 0 : 1;
			ishift = (iindex == 0 ? 5 : 0) + (idstCorner & 1) * 10;
		}
		else
		{
			iindex = 1;
			ishift = idstCorner * 10;
		}
		// fill everything with zero weights
		int izeroW = -1, izeroSeed = -1;
		for (unsigned is = 0; is < stencils.stencils.size(); ++is)
		{
			if (stencils.stencils[is]->weights[idstSeed] == 0)
				izeroSeed = is;
			if (stencils.stencils[is]->weights[idstW] == 0)
				izeroW = is;
			if (izeroSeed != -1 && izeroW != -1)
				break;
		}
		for (unsigned is = 0; is < stencils.stencils.size(); ++is)
		{
			stencils.iVerts[is].is[iindex] &= ~(1023 << ishift);
			stencils.iVerts[is].is[iindex] |= izeroSeed << ishift;
			stencils.iVerts[is].is[iindex] |= izeroW << (ishift + 5);
		}
		for (unsigned is = 0; is < srcPatch.stencils.stencils.size(); ++is)
		{
			const Stencil *pStencilSrc = srcPatch.stencils.stencils[srcPatch.stencils.iVerts[is].is[0] & 31];
			float fSeed = pStencilSrc->weights[isrcSeed];
			int isSeedNew = (fSeed == 0) ? izeroSeed : -1;
			float fW = pStencilSrc->weights[isrcW];
			int isWNew = (fW == 0) ? izeroW : -1;
			Vector3 vPos = pStencilSrc->pos;
			int isPosNew = (fW == 0 && fSeed == 0) ? -2 : -1;
			for (int is = 0; ; ++is)
			{
				if (isSeedNew != -1 && isWNew != -1 && isPosNew != -1)
					break;
			    const Stencil *pStencilDst = stencils.stencils[stencils.iVerts[is].is[0] & 31];
				if (stencils.stencils[is]->weights[idstSeed] == fSeed)
					isSeedNew = is;
				if (stencils.stencils[is]->weights[idstW] == fW)
					isWNew = is;
				if (pStencilDst->pos == vPos)
				    isPosNew = is;
			}
			if (isPosNew != -2)
			{
				stencils.iVerts[isPosNew].is[iindex] &= ~(1023 << ishift);
				stencils.iVerts[isPosNew].is[iindex] |= isSeedNew << ishift;
				stencils.iVerts[isPosNew].is[iindex] |= isWNew << (ishift + 5);
			}
		}
		{
          Vector3 p0 = stencils.EvalPoint(idstSeed, isQuad);
		  Vector3 p1 = srcPatch.stencils.EvalPoint(isrcSeed, srcPatch.isQuad);
		  nvCheck(p0 == p1);
		}
		{
          Vector3 p0 = stencils.EvalPoint(idstW, isQuad);
		  Vector3 p1 = srcPatch.stencils.EvalPoint(isrcW, srcPatch.isQuad);
		  nvCheck(p0 == p1);
		}
	}
}
Vector3 StencilsArray::EvalPoint(int ip, bool bQuad)
{
	Vector3 p(0, 0, 0);
	int iindex = 0;
	int ishift = 0;
	if (bQuad)
	{
		if (ip >= 12)
		{
			if (ip < 16)
			{
				iindex = 0;
				ishift = (ip - 11) * 5;
			}
			else
			{
				iindex = 1;
				ishift = (ip - 16) * 5;
			}
		}
	}
	else
	{
		if (ip >= 9)
		{
			iindex = 1;
			ishift = (ip - 9) * 5;
		}
	}
	for (unsigned is = 0; is < stencils.count(); ++is)
	{
		const Stencil *pStencilPos = stencils[iVerts[is].is[0] & 31];
		const Stencil *pStencilWgt = stencils[(iVerts[is].is[iindex] >> ishift) & 31];
		p += pStencilPos->pos * pStencilWgt->weights[ip];
	}
	return p;
}
static int compare_stencil_ind(const void *_s1, const void *_s2)
{
	const Stencil *s1 = *(Stencil **)_s1;
	const Stencil *s2 = *(Stencil **)_s2;
	nvCheck(s1->indByPos != -1 && s2->indByPos != -1);
	return s1->indByPos - s2->indByPos;
}
static int lexical_compare(const Vector3 &v1, const Vector3 &v2)
{
	if (v1.x() < v2.x())
		return -1;
	if (v1.x() > v2.x())
		return 1;
	if (v1.y() < v2.y())
		return -1;
	if (v1.y() > v2.y())
		return 1;
	if (v1.z() < v2.z())
		return -1;
	if (v1.z() > v2.z())
		return 1;
	return 0;
}
static int compare_stencil_pos(const void *_s1, const void *_s2)
{
	const Stencil *s1 = *(Stencil **)_s1;
	const Stencil *s2 = *(Stencil **)_s2;
	return lexical_compare(s1->pos, s2->pos);
}
void StencilsArray::Sort()
{
	// backup current stencils sorting
	iVerts.resize(stencils.count());
	memcpy(&iVerts[0], &stencils[0], sizeof(stencils[0]) * stencils.count());
	// sort stencils according to 3d pos of their verts
	qsort(&stencils[0], stencils.count(), sizeof(stencils[0]), compare_stencil_pos);
	// memorize positions at which stencils were when sorted by 3d pos
	for (unsigned is = 0; is < stencils.count(); ++is)
	{
		stencils[is]->UpdateSum();
		stencils[is]->indByPos = is;
	}
	// restore backed up stencils sorting
	memcpy(&stencils[0], &iVerts[0], sizeof(stencils[0]) * stencils.count());
	// create order in which we take stencils when performing position evaluation
	for (unsigned is = 0; is < stencils.count(); ++is)
	{
		int iTmp = is;
		iTmp |= iTmp << 5;
		iTmp |= iTmp << 10;
		iTmp |= iTmp << 10;
		iVerts[stencils[is]->indByPos].is[0] = iTmp;
		iVerts[stencils[is]->indByPos].is[1] = iTmp;
	}
}

int CubicGregoryMesh::Patch::iEdge0(int index)
{
	if (isQuad)
	{
		static const int iEdge0[] = { 12, 14, 16, 18 };
		return iEdge0[index];
	}
	else
	{
		static const int iEdge0[] = { 9, 11, 13 };
		return iEdge0[index];
	}
}
int CubicGregoryMesh::Patch::iEdge1(int index)
{
	if (isQuad)
	{
		static const int iEdge1[] = { 13, 15, 17, 19 };
		return iEdge1[index];
	}
	else
	{
		static const int iEdge1[] = { 10, 12, 14 };
		return iEdge1[index];
	}
}
int CubicGregoryMesh::Patch::iCorner(int index)
{
	if (isQuad)
	{
		static const int iCorner[] = { 8, 11, 19, 16 };
		return iCorner[index];
	}
	else
	{
		static const int iCorner[] = { 6, 7, 8 };
		return iCorner[index];
	}
}

// Quick and dirty position evaluation.
Vector3 CubicGregoryMesh::Patch::evalP(float u, float v) const
{
	if (isQuad) 
	{
		const Vector3 * p = positionArray;
		float U, V;

		Vector2 uv = Vector2(2*u - 1, 2*v - 1);
		Vector2 L0 = 0.5 * (1 - uv);
		Vector2 L1 = 0.5 * (1 + uv);

		Vector2 T0 = 0.25 * (1 - uv) * (1 - uv);
		Vector2 T1 = 0.50 * (1 - uv) * (1 + uv);
		Vector2 T2 = 0.25 * (1 + uv) * (1 + uv);

		Vector2 B0 = 0.125 * (1 - uv) * (1 - uv) * (1 - uv);
		Vector2 B1 = 0.375 * (1 + uv) * (1 - uv) * (1 - uv);
		Vector2 B2 = 0.375 * (1 - uv) * (1 + uv) * (1 + uv);
		Vector2 B3 = 0.125 * (1 + uv) * (1 + uv) * (1 + uv);

		u = L0.x(); v=L0.y();
		U = L1.x(); V=L1.y();

		float uu  = u * u,  vv  = v * v;
		float uuu = uu * u, vvv = vv * v;
		float UU  = U * U,  VV  = V * V;
		float UUU = UU * U, VVV = VV * V;

		float denom;

		denom = u + v; if (denom==0) denom = 1; Vector3 b11 = (v*p[0] + u*p[1]) / denom; 
		denom = U + v; if (denom==0) denom = 1; Vector3 b12 = (v*p[3] + U*p[2]) / denom; 
		denom = u + V; if (denom==0) denom = 1; Vector3 b21 = (V*p[4] + u*p[5]) / denom; 
		denom = U + V; if (denom==0) denom = 1; Vector3 b22 = (V*p[7] + U*p[6]) / denom; 

		return 
			(((B3.x() * p[ 8] + B0.x() * p[11]) + (B1.x() * p[10] + B2.x() * p[ 9])) * B3.y()  +
			((B3.x() * p[16] + B0.x() * p[19]) + (B1.x() * p[18] + B2.x() * p[17])) * B0.y()) +
			(((B3.x() * p[12] + B0.x() * p[13]) + (B1.x() *        b12 + B2.x() * b11       )) * B2.y()  +
			((B3.x() * p[14] + B0.x() * p[15]) + (B1.x() *        b22 + B2.x() * b21       )) * B1.y());

	} 
	else
	{
		float w = 1 - u - v;
		float uu = u * u;
		float vv = v * v;
		float ww = w * w;
		float uuu = u * u * u;
		float vvv = v * v * v;
		float www = w * w * w;

		const Vector3 * p = positionArray;

		float  d0 = ( v*(1-u) + u*(1-v) ); if ( (d0<0.001) ) d0 = 1; // avoid NaNs at the corners (value of d* doesn't matter there)
		float  d1 = ( w*(1-v) + v*(1-w) ); if ( (d1<0.001) ) d1 = 1;
		float  d2 = ( u*(1-w) + w*(1-u) ); if ( (d2<0.001) ) d2 = 1;

		Vector3 C0 = ( v*(1-u)*p[5] + u*(1-v)*p[4] ) / d0;
		Vector3 C1 = ( w*(1-v)*p[3] + v*(1-w)*p[2] ) / d1;
		Vector3 C2 = ( u*(1-w)*p[1] + w*(1-u)*p[0] ) / d2;

		return  (p[12]*www + 3*p[11]*ww*u + 3*p[10]*w*uu + p[ 9]*uuu)*(w+u)
			+ (p[ 9]*uuu + 3*p[ 8]*uu*v + 3*p[ 7]*u*vv + p[ 6]*vvv)*(u+v)
			+ (p[ 6]*vvv + 3*p[14]*vv*w + 3*p[13]*v*ww + p[12]*www)*(v+w)
			- (p[12]*www*w + p[9]*uuu*u + p[6]*vvv*v) 
			+  12*(C0*u*v*ww + C1*uu*v*w   + C2*u*vv*w);
	}
}




// Quick and dirty tangent evaluation.
Vector3 CubicGregoryMesh::Patch::evalTV(float u, float v) const
{
	const Vector3 * p = positionArray;

	if (isQuad) {
		float U, V;
		Vector2 uv = Vector2(2*u-1, 2*v-1);
		Vector2 L0 = 0.5 * (1 - uv);
		Vector2 L1 = 0.5 * (1 + uv);

		Vector2 T0 = 0.25 * (1 - uv) * (1 - uv);
		Vector2 T1 = 0.50 * (1 - uv) * (1 + uv);
		Vector2 T2 = 0.25 * (1 + uv) * (1 + uv);

		Vector2 B0 = 0.125 * (1 - uv) * (1 - uv) * (1 - uv);
		Vector2 B1 = 0.375 * (1 + uv) * (1 - uv) * (1 - uv);
		Vector2 B2 = 0.375 * (1 - uv) * (1 + uv) * (1 + uv);
		Vector2 B3 = 0.125 * (1 + uv) * (1 + uv) * (1 + uv);

		u = L0.x(), v=L0.y();
		U = L1.x(), V=L1.y();

		float uu  = u * u,  vv  = v * v;
		float uuu = uu * u, vvv = vv * v;
		float UU  = U * U,  VV  = V * V;
		float UUU = UU * U, VVV = VV * V;

		float denom;

		denom = u + v; if (denom==0) denom = 1; Vector3 b11 = (v*p[0] + u*p[1]) / denom; 
		denom = U + v; if (denom==0) denom = 1; Vector3 b12 = (v*p[3] + U*p[2]) / denom; 
		denom = u + V; if (denom==0) denom = 1; Vector3 b21 = (V*p[4] + u*p[5]) / denom; 
		denom = U + V; if (denom==0) denom = 1; Vector3 b22 = (V*p[7] + U*p[6]) / denom; 

		Vector3 tv00 = (p[12] - p[8 ]);
		Vector3 tv01 = (p[14] - p[12]);
		Vector3 tv02 = (p[16] - p[14]);

		Vector3 tv10 = (       b11 - p[9 ]);
		Vector3 tv11 = (       b21 -        b11);
		Vector3 tv12 = (p[17] -        b21);

		Vector3 tv20 = (       b12 - p[10]);
		Vector3 tv21 = (       b22 -        b12);
		Vector3 tv22 = (p[18] -        b22);

		Vector3 tv30 = (p[13] - p[11]);
		Vector3 tv31 = (p[15] - p[13]);
		Vector3 tv32 = (p[19] - p[15]);

		return 3 * (
			(((tv00 * T2.y() + tv02 * T0.y()) + tv01 * T1.y()) * B3.x() +
			( (tv30 * T2.y() + tv32 * T0.y()) + tv31 * T1.y()) * B0.x() ) +
			(((tv10 * T2.y() + tv12 * T0.y()) + tv11 * T1.y()) * B2.x() +
			( (tv20 * T2.y() + tv22 * T0.y()) + tv21 * T1.y()) * B1.x() ) );

	} else {

		float w = 1 - u - v;
		float uu = u * u;
		float vv = v * v;
		float ww = w * w;
		float uuu = uu * u;
		float vvv = vv * v;
		float www = ww * w;

		float d;

		d = ( v*(1-u) + u*(1-v) ); if (d==0) d = 1; Vector3 C0 = ( v*(1-u)*p[5] + u*(1-v)*p[4] ) / d;
		d = ( w*(1-v) + v*(1-w) ); if (d==0) d = 1; Vector3 C1 = ( w*(1-v)*p[3] + v*(1-w)*p[2] ) / d;
		d = ( u*(1-w) + w*(1-u) ); if (d==0) d = 1; Vector3 C2 = ( u*(1-w)*p[1] + w*(1-u)*p[0] ) / d;

		Vector3 E1 = (p[12]*www + 3*p[11]*ww*u + 3*p[10]*w*uu + p[ 9]*uuu);
		Vector3 E2 = (p[ 9]*uuu + 3*p[ 8]*uu*v + 3*p[ 7]*u*vv + p[ 6]*vvv);
		Vector3 E3 = (p[ 6]*vvv + 3*p[14]*vv*w + 3*p[13]*v*ww + p[12]*www);

		Vector3 E1v = 3*( -p[12]*ww - 2*p[11]*w*u        -   p[10]*uu         );
		Vector3 E2v = 3*(               p[ 8]*uu         + 2*p[ 7]*u*v        + p[ 6]*vv);
		Vector3 E3v = 3*( p[ 6]*vv +    p[14]*(2*w*v-vv) +   p[13]*(ww-2*w*v) - p[12]*ww);
		Vector3 Sv  = 4*( -p[12]*www + p[6]*vvv);
		Vector3 Cv  = 12*(C0*(u*ww-2*u*v*w) + C1*uu*(w-v) + C2*(2*u*v*w-u*vv));

		return - ((E1v*(w+u)-E1) + (E2+E2v*(u+v)) + E3v*(v+w) - Sv + Cv );
	}

}



// Quick and dirty tangent evaluation.
Vector3 CubicGregoryMesh::Patch::evalTU(float u, float v) const
{
	const Vector3 * p = positionArray;

	if (isQuad)
	{
		float U, V;
		Vector2 uv = Vector2(2*u-1, 2*v-1);
		Vector2 L0 = 0.5 * (1 - uv);
		Vector2 L1 = 0.5 * (1 + uv);

		Vector2 T0 = 0.25 * (1 - uv) * (1 - uv);
		Vector2 T1 = 0.50 * (1 - uv) * (1 + uv);
		Vector2 T2 = 0.25 * (1 + uv) * (1 + uv);

		Vector2 B0 = 0.125 * (1 - uv) * (1 - uv) * (1 - uv);
		Vector2 B1 = 0.375 * (1 + uv) * (1 - uv) * (1 - uv);
		Vector2 B2 = 0.375 * (1 - uv) * (1 + uv) * (1 + uv);
		Vector2 B3 = 0.125 * (1 + uv) * (1 + uv) * (1 + uv);

		u = L0.x(); v=L0.y();
		U = L1.x(); V=L1.y();

		float uu  = u * u,  vv  = v * v;
		float uuu = uu * u, vvv = vv * v;
		float UU  = U * U,  VV  = V * V;
		float UUU = UU * U, VVV = VV * V;

		float denom;

		denom = u + v; if (denom==0) denom = 1; Vector3 b11 = (v*p[0] + u*p[1]) / denom; 
		denom = U + v; if (denom==0) denom = 1; Vector3 b12 = (v*p[3] + U*p[2]) / denom; 
		denom = u + V; if (denom==0) denom = 1; Vector3 b21 = (V*p[4] + u*p[5]) / denom; 
		denom = U + V; if (denom==0) denom = 1; Vector3 b22 = (V*p[7] + U*p[6]) / denom; 


		Vector3 tu00 = (p[9 ] - p[ 8]);  
		Vector3 tu01 = (p[10] - p[ 9]);
		Vector3 tu02 = (p[11] - p[10]);

		Vector3 tu10 = (       b11 - p[12]);
		Vector3 tu11 = (       b12 -        b11);
		Vector3 tu12 = (p[13] -        b12);

		Vector3 tu20 = (       b21 - p[14]);
		Vector3 tu21 = (       b22 -        b21);
		Vector3 tu22 = (p[15] -        b22);

		Vector3 tu30 = (p[17] - p[16]);
		Vector3 tu31 = (p[18] - p[17]);
		Vector3 tu32 = (p[19] - p[18]);

		return 3 * (
			(((tu00 * T2.x() + tu02 * T0.x()) + tu01 * T1.x()) * B3.y() +
			( (tu30 * T2.x() + tu32 * T0.x()) + tu31 * T1.x()) * B0.y() ) +
			(((tu10 * T2.x() + tu12 * T0.x()) + tu11 * T1.x()) * B2.y() +
			( (tu20 * T2.x() + tu22 * T0.x()) + tu21 * T1.x()) * B1.y() ) );




		/*		
		float uu = u * u;
		float vv = v * v;
		float uuu = u * u * u;
		float vvv = v * v * v;
		float U = 1-u;
		float V = 1-v;
		float UU = U * U;
		float VV = V * V;
		float VVV = VV * V;

		Vector3 b11 = (v*p[0] + u*p[1]) / (u + v); if (u+v==0) b11=Vector3(0,0,0);
		Vector3 b12 = (v*p[3] + U*p[2]) / (U + v); if (U+v==0) b12=Vector3(0,0,0);
		Vector3 b21 = (V*p[4] + u*p[5]) / (u + V); if (u+V==0) b21=Vector3(0,0,0);
		Vector3 b22 = (V*p[7] + U*p[6]) / (U + V); if (U+V==0) b22=Vector3(0,0,0); 

		Vector3 b11 = (V*p[0] + U*p[1]) / (U + V); if (U+V==0) b11=Vector3(0,0,0);
		Vector3 b12 = (V*p[3] + u*p[2]) / (u + V); if (u+V==0) b12=Vector3(0,0,0);
		Vector3 b21 = (v*p[4] + U*p[5]) / (U + v); if (v+U==0) b21=Vector3(0,0,0);
		Vector3 b22 = (v*p[7] + u*p[6]) / (u + v); if (u+v==0) b22=Vector3(0,0,0); 

		Vector3 tu00 = (p[9 ] - p[ 8]);  
		Vector3 tu01 = (p[10] - p[ 9]);
		Vector3 tu02 = (p[11] - p[10]);

		Vector3 tu10 = (  b11 - p[12]);
		Vector3 tu11 = (  b12 -   b11);
		Vector3 tu12 = (p[13] -   b12);

		Vector3 tu20 = (  b21 - p[14]);
		Vector3 tu21 = (  b22 -   b21);
		Vector3 tu22 = (p[15] -   b22);

		Vector3 tu30 = (p[17] - p[16]);
		Vector3 tu31 = (p[18] - p[17]);
		Vector3 tu32 = (p[19] - p[18]);

		return 3*(
		(tu00 * UU + tu01 * u*U + tu02 * uu) * VVV +
		(tu10 * UU + tu11 * u*U + tu12 * uu) * VV * v * 3 +
		(tu20 * UU + tu21 * u*U + tu22 * uu) * V * vv * 3 +
		(tu31 * UU + tu32 * u*U + tu32 * uu) * vvv ); */

	}
	else
	{
		float w = 1 - u - v;

		float uu = u * u;
		float vv = v * v;
		float ww = w * w;
		float uuu = u * u * u;
		float vvv = v * v * v;
		float www = w * w * w;

		float d;

		d = ( v*(1-u) + u*(1-v) ); if (d==0) d = 1; Vector3 C0 = ( v*(1-u)*p[5] + u*(1-v)*p[4] ) / d;
		d = ( w*(1-v) + v*(1-w) ); if (d==0) d = 1; Vector3 C1 = ( w*(1-v)*p[3] + v*(1-w)*p[2] ) / d;
		d = ( u*(1-w) + w*(1-u) ); if (d==0) d = 1; Vector3 C2 = ( u*(1-w)*p[1] + w*(1-u)*p[0] ) / d;

		Vector3 E1 = (p[12]*www + 3*p[11]*ww*u + 3*p[10]*w*uu + p[ 9]*uuu);
		Vector3 E2 = (p[ 9]*uuu + 3*p[ 8]*uu*v + 3*p[ 7]*u*vv + p[ 6]*vvv);
		Vector3 E3 = (p[ 6]*vvv + 3*p[14]*vv*w + 3*p[13]*v*ww + p[12]*www);

		Vector3 E1u = 3*( -p[12]*ww + p[11]*(ww-2*u*w) +   p[10]*(2*u*w-uu) + p[ 9]*uu);
		Vector3 E2u = 3*(  p[ 9]*uu + 2*p[ 8]*u*v      +   p[ 7]*vv         );
		Vector3 E3u = 3*(           -   p[14]*vv       - 2*p[13]*v*w        - p[12]*ww);
		Vector3 Su  = 4*( -p[12]*www + p[9]*uuu);
		Vector3 Cu  = 12*( C0*(ww*v-2*u*v*w) + C1*(2*u*v*w-uu*v) + C2*vv*(w-u));

		return  E1u*(w+u) + (E2+E2u*(u+v)) + (E3u*(v+w)-E3) - Su + Cu;
	} 
}

/*
// de casteljau's algorithm for triangles
static Vector3 DeCasteljau(float u, float v, const Vector3 values[10])
{
	float w = 1.0f - float(u) - float(v);
	Vector3 t0[6], t1[3];
	t0[0] = values[0] * v + values[1] * u + values[8] * w;
	t0[1] = values[1] * v + values[2] * u + values[9] * w;
	t0[2] = values[2] * v + values[3] * u + values[4] * w;
	t0[3] = values[9] * v + values[4] * u + values[5] * w;
	t0[4] = values[7] * v + values[5] * u + values[6] * w;
	t0[5] = values[8] * v + values[9] * u + values[7] * w;
	t1[0] = t0[0] * v + t0[1] * u + t0[5] * w;
	t1[1] = t0[1] * v + t0[2] * u + t0[3] * w;
	t1[2] = t0[5] * v + t0[3] * u + t0[4] * w;
	return t1[0] * v + t1[1] * u + t1[2] * w;
}

Vector3 CubicGregoryMesh::Patch::evalTU0(float u, float v) const
{
	return DeCasteljau(u, v, tangents);
}
Vector3 CubicGregoryMesh::Patch::evalTV0(float u, float v) const
{
	return DeCasteljau(u, v, bitangents);
}
*/