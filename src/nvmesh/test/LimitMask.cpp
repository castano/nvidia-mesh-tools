
#include <iostream>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

/*
void test()
{
	Vector3 v(1, 1, 1);
	
	for(int i = 0; i < 2; i++)
	{
		float v0 = 0.5 * (v.x() + v.y());
		float v2 = 0.5 * (v.y() + v.z());
		float v1 = 0.25 * (v0 + v2) + 0.5 * v.y();
		v.set(v0, v1, v2);
	}
	
	printf("v = %f %f %f\n", v.x(), v.y(), v.z());
}
*/

// Compute the limit mask of the cubic B-Spline.
int main(void)
{
	symbol sv0("v0");
	symbol sv1("v1");
	symbol sv2("v2");
	ex v0 = sv0;
	ex v1 = sv1;
	ex v2 = sv2;
	
	for(int i = 0; i < 10; i++)
	{
		v0 = (v0 + v1) / 2;
		v2 = (v1 + v2) / 2;
		v1 = (v0 + v2) / 4 + v1 / 2;
	}
	
//	cout << "v0 = " << evalf(v0) << endl;
//	cout << "v1 = " << evalf(v1) << endl;
//	cout << "v2 = " << evalf(v2) << endl;
	
	cout << evalf(v1.subs(sv0 == 1).subs(sv1 == 0).subs(sv2 == 0)) << endl;
	cout << evalf(v1.subs(sv0 == 0).subs(sv1 == 1).subs(sv2 == 0)) << endl;
	cout << evalf(v1.subs(sv0 == 0).subs(sv1 == 0).subs(sv2 == 1)) << endl;

	matrix A(4,4);
	A = -1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 0, 3, 0,
		1, 4, 1, 0;
	A = A.mul_scalar(ex(1)/ex(6));

	matrix B(4,4);
	B = -1, 3, -3, 1,
	  	3, -6, 3, 0,
		-3, 3, 0, 0,
		1, 0, 0, 0;

	cout << eval(B.inverse().mul(A)) << endl;


	return 0;
}


