// This code is in the public domain -- castanyo@yahoo.es

#include "Bezier.h"


using namespace nv;


static void deCasteljau(float u, Vector3::Arg p0, Vector3::Arg p1, Vector3::Arg p2, Vector3::Arg p3, Vector3 * p, Vector3 * dp)
{
	Vector3 q0 = lerp(p0, p1, u);
	Vector3 q1 = lerp(p1, p2, u);
	Vector3 q2 = lerp(p2, p3, u);
	Vector3 r0 = lerp(q0, q1, u);
	Vector3 r1 = lerp(q1, q2, u);
	
	if (dp != NULL) *dp = r0 - r1;
	if (p != NULL) *p = lerp(r0, r1, u);
}


void nv::evaluateCubicBezierPatch(float u, float v, const Vector3 * cp, Vector3 * pos, Vector3 * du, Vector3 * dv)
{
#if 0
	Vector2 L0(1-u,1-v);
	Vector2 L1(u,v);

	Vector2 Q0 =     L0 * L0;
	Vector2 Q1 = 2 * L0 * L1;
	Vector2 Q2 =     L1 * L1;

	Vector2 B0 =     L0 * L0 * L0;
	Vector2 B1 = 3 * L0 * L0 * L1;
	Vector2 B2 = 3 * L1 * L1 * L0;
	Vector2 B3 =     L1 * L1 * L1;

	*pos = 
		(B0.x() * cp[ 0] + B1.x() * cp[ 1] + B2.x() * cp[ 2] + B3.x() * cp[ 3]) * B0.y() +
		(B0.x() * cp[ 4] + B1.x() * cp[ 5] + B2.x() * cp[ 6] + B3.x() * cp[ 7]) * B1.y() +
		(B0.x() * cp[ 8] + B1.x() * cp[ 9] + B2.x() * cp[10] + B3.x() * cp[11]) * B2.y() +
		(B0.x() * cp[12] + B1.x() * cp[13] + B2.x() * cp[14] + B3.x() * cp[15]) * B3.y();

	*du =
		((cp[0]-cp[1]) * B0.y() + (cp[4]-cp[5]) * B1.y() + (cp[ 8]-cp[ 9]) * B2.y() + (cp[12]-cp[13]) * B3.y()) * Q0.x() +
		((cp[1]-cp[2]) * B0.y() + (cp[5]-cp[6]) * B1.y() + (cp[ 9]-cp[10]) * B2.y() + (cp[13]-cp[14]) * B3.y()) * Q1.x() +
		((cp[2]-cp[3]) * B0.y() + (cp[6]-cp[7]) * B1.y() + (cp[10]-cp[11]) * B2.y() + (cp[14]-cp[15]) * B3.y()) * Q2.x();

	*dv =
		((cp[0]-cp[ 4]) * B0.x() + (cp[1]-cp[ 5]) * B1.x() + (cp[ 2]-cp[ 6]) * B2.x() + (cp[ 3]-cp[ 7]) * B3.x()) * Q0.y() +
		((cp[4]-cp[ 8]) * B0.x() + (cp[5]-cp[ 9]) * B1.x() + (cp[ 6]-cp[10]) * B2.x() + (cp[ 7]-cp[11]) * B3.x()) * Q1.y() +
		((cp[8]-cp[12]) * B0.x() + (cp[9]-cp[13]) * B1.x() + (cp[10]-cp[14]) * B2.x() + (cp[11]-cp[15]) * B3.x()) * Q2.y();
#else
	Vector3 t0, t1, t2, t3;
	Vector3 q0, q1, q2, q3;

	deCasteljau(u, cp[ 0], cp[ 1], cp[ 2], cp[ 3], &q0, &t0);
	deCasteljau(u, cp[ 4], cp[ 5], cp[ 6], cp[ 7], &q1, &t1);
	deCasteljau(u, cp[ 8], cp[ 9], cp[10], cp[11], &q2, &t2);
	deCasteljau(u, cp[12], cp[13], cp[14], cp[15], &q3, &t3);

	deCasteljau(v, q0, q1, q2, q3, pos, dv);
	deCasteljau(v, t0, t1, t2, t3, du, NULL);
#endif
}

static float noZeroDiv(float f)
{
	if (equal(f, 0, 0.0f)) return 1.0f;
	return f;
}

void nv::evaluateQuadGregoryPatch(float u, float v, const Vector3 cp[20], Vector3 * pos, Vector3 * du, Vector3 * dv)
{
    Vector3 q[16];

    float U = 1 - u;
    float V = 1 - v;

    //  8     9     10     11
    // 12   0\1     2/3    13
    // 14   4/5     6\7    15
    // 16    17     18     19

    q[ 5] = (u * cp[1] + v * cp[0]) / noZeroDiv(u + v);
    q[ 6] = (U * cp[2] + v * cp[3]) / noZeroDiv(U + v);
    q[ 9] = (u * cp[5] + V * cp[4]) / noZeroDiv(u + V);
    q[10] = (U * cp[6] + V * cp[7]) / noZeroDiv(U + V);

    // Map gregory control points to bezier control points.
    q[ 0] = cp[8];
    q[ 1] = cp[9];
    q[ 2] = cp[10];
    q[ 3] = cp[11];
    q[ 4] = cp[12];
    q[ 7] = cp[13];
    q[ 8] = cp[14];
    q[11] = cp[15];
    q[12] = cp[16];
    q[13] = cp[17];
    q[14] = cp[18];
    q[15] = cp[19];

    evaluateCubicBezierPatch(u, v, q, pos, du, dv);
}


void nv::evaluateTriangleGregoryPatch(float u, float v, const Vector3 p[15], Vector3 * pos, Vector3 * du, Vector3 * dv)
{
#pragma message(NV_FILE_LINE "This needs to be cleaned up optimized!")

    {
	    float w = 1 - u - v;
	    float uu = u * u;
	    float vv = v * v;
	    float ww = w * w;
	    float uuu = u * u * u;
	    float vvv = v * v * v;
	    float www = w * w * w;

	    float U = 1 - u;
	    float V = 1 - v;
	    float W = 1 - w;

	    //              6
        //              
	    //        14   0/1   7
	    //                          
	    //    13   5/4     3\2   8
	    //               
	    // 12      11       10      9

	    Vector3 C0 = ( v*U * p[5] + u*V * p[4] ) / noZeroDiv(v*U + u*V);
	    Vector3 C1 = ( w*V * p[3] + v*W * p[2] ) / noZeroDiv(w*V + v*W);
	    Vector3 C2 = ( u*W * p[1] + w*U * p[0] ) / noZeroDiv(u*W + w*U);

	    *pos =
		    (p[12] * www + 3*p[11] * ww*u + 3*p[10] * w*uu + p[ 9]*uuu) * (w + u) +
		    (p[ 9] * uuu + 3*p[ 8] * uu*v + 3*p[ 7] * u*vv + p[ 6]*vvv) * (u + v) +
		    (p[ 6] * vvv + 3*p[14] * vv*w + 3*p[13] * v*ww + p[12]*www) * (v + w) -
		    (p[12] * www*w + p[ 9] * uuu*u + p[ 6] * vvv*v) +
	        12*(C0 * u*v*ww + C1 * uu*v*w   + C2 * u*vv*w);
    }
    {
		float w = 1 - u - v;
		float uu = u * u;
		float vv = v * v;
		float ww = w * w;
		float uuu = uu * u;
		float vvv = vv * v;
		float www = ww * w;

		float U = 1 - u;
		float V = 1 - v;
		float W = 1 - w;

		Vector3 C0 = ( v*U * p[5] + u*V * p[4] ) / noZeroDiv(v*U + u*V);
		Vector3 C1 = ( w*V * p[3] + v*W * p[2] ) / noZeroDiv(w*V + v*W);
		Vector3 C2 = ( u*W * p[1] + w*U * p[0] ) / noZeroDiv(u*W + w*U);

		Vector3 E1 = (p[12]*www + 3*p[11]*ww*u + 3*p[10]*w*uu + p[ 9]*uuu);
		Vector3 E2 = (p[ 9]*uuu + 3*p[ 8]*uu*v + 3*p[ 7]*u*vv + p[ 6]*vvv);
		Vector3 E3 = (p[ 6]*vvv + 3*p[14]*vv*w + 3*p[13]*v*ww + p[12]*www);

		Vector3 E1u = 3*( - p[12]*ww + p[11]*(ww-2*u*w) +   p[10]*(2*u*w-uu) + p[ 9]*uu);
		Vector3 E2u = 3*(   p[ 9]*uu + 2*p[ 8]*u*v      +   p[ 7]*vv         );
		Vector3 E3u = 3*(            - p[14]*vv         - 2*p[13]*v*w        - p[12]*ww);
		Vector3 Su  = 4*( -p[12]*www + p[9]*uuu);
		Vector3 Cu  = 12*( C0*(ww*v-2*u*v*w) + C1*(2*u*v*w-uu*v) + C2*vv*(w-u) );

		*du = E1u*(w+u) + (E2+E2u*(u+v)) + (E3u*(v+w)-E3) - Su + Cu;
    }
    {
		float w = 1 - u - v;
		float uu = u * u;
		float vv = v * v;
		float ww = w * w;
		float uuu = uu * u;
		float vvv = vv * v;
		float www = ww * w;

		float U = 1 - u;
		float V = 1 - v;
		float W = 1 - w;

		Vector3 C0 = ( v*U * p[5] + u*V * p[4] ) / noZeroDiv(v*U + u*V);
		Vector3 C1 = ( w*V * p[3] + v*W * p[2] ) / noZeroDiv(w*V + v*W);
		Vector3 C2 = ( u*W * p[1] + w*U * p[0] ) / noZeroDiv(u*W + w*U);

		Vector3 E1 = (p[12]*www + 3*p[11]*ww*u + 3*p[10]*w*uu + p[ 9]*uuu);
		Vector3 E2 = (p[ 9]*uuu + 3*p[ 8]*uu*v + 3*p[ 7]*u*vv + p[ 6]*vvv);
		Vector3 E3 = (p[ 6]*vvv + 3*p[14]*vv*w + 3*p[13]*v*ww + p[12]*www);

		Vector3 E1v = 3*(-p[12]*ww  - 2*p[11]*w*u       -   p[10]*uu         );
		Vector3 E2v = 3*(              p[ 8]*uu         + 2*p[ 7]*u*v        + p[ 6]*vv);
		Vector3 E3v = 3*( p[ 6]*vv  +  p[14]*(2*w*v-vv) +   p[13]*(ww-2*w*v) - p[12]*ww);
		Vector3 Sv  = 4*(-p[12]*www +  p[ 6]*vvv);
		Vector3 Cv  = 12*(C0*(u*ww-2*u*v*w) + C1*uu*(w-v) + C2*(2*u*v*w-u*vv));

		*dv = ((E1v*(w+u)-E1) + (E2+E2v*(u+v)) + E3v*(v+w) - Sv + Cv );
    }
}


static void evaluateQuarticTriangle(float u, float v, float w, const Vector3 q[15], Vector3 * pos, Vector3 * du, Vector3 * dv)
{
    Vector3 p0[10];
    
    uint k = 0;
    for (int j = 0; j < 4; j++) {
        p0[k++] = u*q[j] + v*q[j+1] + w*q[j+5];
    }
    for (int j = 5; j < 8; j++) {
        p0[k++] = u*q[j] + v*q[j+1] + w*q[j+4];
    }
    for (int j = 9; j < 11; j++) {
        p0[k++] = u*q[j] + v*q[j+1] + w*q[j+3];
    }	 
    p0[9] = u*q[12] + v*q[13] + w*q[14];
    
    Vector3 p1[6];
    k = 0;
    for (int j = 0; j < 3; j++) {
        p1[k++] = u*p0[j] + v*p0[j+1] + w*p0[j+4];
    } 
    for (int j = 4; j < 6; j++) {
        p1[k++] = u*p0[j] + v*p0[j+1] + w*p0[j+3];
    } 	
    p1[5] = u*p0[7] + v*p0[8] + w*p0[9]; 

    Vector3 p2[3];
    for (int j = 0; j < 2; j++) {
        p2[j] = u*p1[j] + v*p1[j+1] + w*p1[j+3];
    }  
    p2[2] = u*p1[3] + v*p1[4] + w*p1[5];
    
    if (pos) *pos = u * p2[0] + v * p2[1] + w * p2[2];
    if (du) *du = p2[0] - p2[2];
    if (dv) *dv = p2[1] - p2[0];
}


// Get control points of a quartic triangle.
static void getPmControlPoints(const Vector3 p[24], Vector3 q[15], int sector)
{
    int next = (sector + 1) % 4;
    int prev = (sector + 3) % 4;
    int opposite = (sector + 2) % 4;
    
    q[ 0] = p[6 * sector];
	q[ 1] = lerp(p[6*sector + 0], p[6*sector+1], 0.75f);
	q[ 2] = lerp(p[6*sector + 1], p[6*sector+2], 0.5f);
	q[ 3] = lerp(p[6*next      ], p[6*sector+2], 0.75f);
	q[ 4] = p[6 * next];
	q[ 6] = p[6 * sector + 3];
	q[ 7] = p[6 * sector + 4];
	q[10] = p[6 * sector + 5];

	q[ 5] = (q[ 1] + lerp(q[0], p[6 * prev + 2], 0.75f)) * 0.5f;
	q[ 8] = (q[ 3] + lerp(q[4], p[6 * next + 1], 0.75f)) * 0.5f;
	q[ 9] = (q[ 6] + p[6*prev+4]) * 0.5f;
	q[11] = (q[ 7] + p[6*next+3]) * 0.5f;
	q[12] = (q[10] + p[6*prev+5]) * 0.5f;
	q[13] = (q[10] + p[6*next+5]) * 0.5f;
	q[14] = (q[10] + p[6*next+5] + p[6*prev+5] + p[6*opposite+5]) * 0.25f;
}


void nv::evaluateQuadPmPatch(float u, float v, const Vector3 cp[24], Vector3 * pos, Vector3 * du, Vector3 * dv)
{
	Vector3 q[15];

    float x, y, z;

	if (v <= u && u+v <= 1) {
        x = 1-u-v;
        y = u-v;
        z = 2*v;
        getPmControlPoints(cp, q, 0);
	}
	else if (v <= u && u+v > 1) {
        x = u-v;
        y = u+v-1;
        z = 2-2*u;
        getPmControlPoints(cp, q, 1);
	}
	else if (v > u && u+v <= 1) {
        x = v-u;
        y = 1-v-u;
        z = 2*u;
	    getPmControlPoints(cp, q, 3);
	}
	else {
	    x = u+v-1;
        y = v-u;
        z = 2-2*v;
	    getPmControlPoints(cp, q, 2);
	}

	evaluateQuarticTriangle(x, y, z, q, pos, du, dv);
}

void nv::evaluateTrianglePmPatch(float u, float v, const Vector3 cp[18], Vector3 * pos, Vector3 * du, Vector3 * dv)
{
}
