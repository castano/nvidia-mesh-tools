// This code is in the public domain -- castanyo@yahoo.es

#ifndef NV_MESH_BONE_H
#define NV_MESH_BONE_H

#include <nvmath/Matrix.h>
#include <nvmath/Quaternion.h>


namespace nv
{
	/// Build bone matrix from quatertion and offset.
	inline Matrix boneMatrix(Quaternion::Arg q, Vector3::Arg offset)
	{
		// calculate coefficients
		float x2 = 2 * q.x();
		float y2 = 2 * q.y();
		float z2 = 2 * q.z();
		
		float xx, xy, xz, yy, yz, zz, wx, wy, wz;
		xx = q.x() * x2;   xy = q.x() * y2;   xz = q.x() * z2;
		yy = q.y() * y2;   yz = q.y() * z2;   zz = q.z() * z2;
		wx = q.w() * x2;   wy = q.w() * y2;   wz = q.w() * z2;
		
		Matrix m;
		m.data(0) = 1.0f - (yy + zz); 	
		m.data(1) = xy - wz;
		m.data(2) = xz + wy;		
		m.data(3) = 0.0f;
 		
		m.data(4) = xy + wz;		
		m.data(5) = 1.0f - (xx + zz);
		m.data(6) = yz - wx;		
		m.data(7) = 0.0f;
		
		m.data(8) = xz - wy;		
		m.data(9) = yz + wx;
		m.data(10) = 1.0f - (xx + yy);		
		m.data(11) = 0.0f;
		
		m.data(12) = offset.x();
		m.data(13) = offset.y();
		m.data(14) = offset.z();			
		m.data(15) = 1.0f;
		
		return m;
	}
}


#endif // NV_MESH_BONE_H
