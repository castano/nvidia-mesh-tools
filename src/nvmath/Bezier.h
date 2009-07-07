// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MATH_BEZIER_H
#define NV_MATH_BEZIER_H

#include <nvmath/nvmath.h>
#include <nvmath/Vector.h>

namespace nv
{

	void evaluateCubicBezierPatch(float u, float v, const Vector3 cp[16], Vector3 * pos, Vector3 * du, Vector3 * dv);

	void evaluateQuadGregoryPatch(float u, float v, const Vector3 cp[20], Vector3 * pos, Vector3 * du, Vector3 * dv);

	void evaluateTriangleGregoryPatch(float u, float v, const Vector3 cp[15], Vector3 * pos, Vector3 * du, Vector3 * dv);

	void evaluateQuadPmPatch(float u, float v, const Vector3 cp[24], Vector3 * pos, Vector3 * du, Vector3 * dv);

	void evaluateTrianglePmPatch(float u, float v, const Vector3 cp[18], Vector3 * pos, Vector3 * du, Vector3 * dv);


} // nv namespace

#endif // NV_MATH_BEZIER_H
