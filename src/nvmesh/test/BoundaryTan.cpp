
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

void derivative()
{
	symbol t("t");

	symbol v0("v0");
	symbol v1("v1");
	symbol v2("v2");
	symbol v3("v3");

	ex T = (1-t);
	ex v = v0 * T*T*T + v1 * 3 * T*T*t + v2 * 3 * T*t*t + v3 * t*t*t;

	ex dv = eval(v.diff(t));

	cout << "dv/dt = " << dv << endl;
	cout << endl;

	cout << "dv/dt(0) = " << dv.subs(t == 0) << endl;
	cout << endl;
}


// test linear independance
void test2()
{
	symbol t("t");

	symbol u0("u0");
	symbol u1("u1");
	symbol u2("u2");

	symbol v0("v0");
	symbol v1("v1");
	symbol v2("v2");
	symbol v3("v3");

	symbol w0("w0");
	symbol w1("w1");
	symbol w2("w2");
	symbol w3("w3");

	symbol i_sym("i");	// x,y,z
	idx i(i_sym, 3);

	symbol x("x");
	symbol y("y");
	symbol z("z");

	ex T = (1-t);

	ex u = indexed(u0, i) * T*T + indexed(u1, i) * 2 * T*t + indexed(u2, i) * t * t;
	ex v = indexed(v0, i) * T*T*T + indexed(v1, i) * 3 * T*T*t + indexed(v2, i) * 3 * T*t*t + indexed(v3, i) * t*t*t;
	ex w = indexed(w0, i) * T*T*T + indexed(w1, i) * 3 * T*T*t + indexed(w2, i) * 3 * T*t*t + indexed(w3, i) * t*t*t;

	cout << "u = " << u << endl;
	cout << "v = " << v << endl;
	cout << "w = " << w << endl;
	cout << endl;

	// u,v,w must be linearly dependant!
	matrix A(3, 3); 
	A = u.subs(i == x), u.subs(i == y), u.subs(i == z),
		v.subs(i == x), v.subs(i == y), v.subs(i == z),
		w.subs(i == x), w.subs(i == y), w.subs(i == z);

	cout << "A = " << A << endl;
	cout << endl;

	cout << "|A| = " << eval(determinant(A)) << endl;
	cout << endl;

	cout << "|A|0 = " << eval(determinant(A).subs(t == 0)) << endl;
	cout << "|A|1 = " << eval(determinant(A).subs(t == 1)) << endl;
}

// Compute tangent boundaries.
int main(void)
{
	derivative();

/*
	symbol t("t");

	symbol u0("u0");
	symbol u1("u1");
	symbol u2("u2");

	symbol v0("v0");
	symbol v1("v1");
	symbol v2("v2");
	symbol v3("v3");

	symbol w0("w0");
	symbol w1("w1");
	symbol w2("w2");
	symbol w3("w3");

	symbol i_sym("i");	// x,y,z
	idx i(i_sym, 3);

	symbol x("x");
	symbol y("y");
	symbol z("z");
*/
/*
	matrix u0(3,1); u0 = symbol("u0x"), symbol("u0y"), symbol("u0z");
	matrix u1(3,1); u1 = symbol("u1x"), symbol("u1y"), symbol("u1z");
	matrix u2(3,1); u2 = symbol("u2x"), symbol("u2y"), symbol("u2z");

	matrix v0(3,1); v0 = symbol("v0x"), symbol("v0y"), symbol("v0z");
	matrix v1(3,1); v1 = symbol("v1x"), symbol("v1y"), symbol("v1z");
	matrix v2(3,1); v2 = symbol("v2x"), symbol("v2y"), symbol("v2z");
	matrix v3(3,1); v3 = symbol("v3x"), symbol("v3y"), symbol("v3z");

	matrix w0(3,1); w0 = symbol("w0x"), symbol("w0y"), symbol("w0z");
	matrix w1(3,1); w1 = symbol("w1x"), symbol("w1y"), symbol("w1z");
	matrix w2(3,1); w2 = symbol("w2x"), symbol("w2y"), symbol("w2z");
	matrix w3(3,1); w3 = symbol("w3x"), symbol("w3y"), symbol("w3z");
*/
/*
	ex T = (1-t);

	ex u = indexed(u0, i) * T*T + indexed(u1, i) * 2 * T*t + indexed(u2, i) * t * t;
	ex v = indexed(v0, i) * T*T*T + indexed(v1, i) * 3 * T*T*t + indexed(v2, i) * 3 * T*t*t + indexed(v3, i) * t*t*t;
	ex w = indexed(w0, i) * T*T*T + indexed(w1, i) * 3 * T*T*t + indexed(w2, i) * 3 * T*t*t + indexed(w3, i) * t*t*t;

	cout << "u = " << u << endl;
	cout << "v = " << v << endl;
	cout << "w = " << w << endl;
	cout << endl;

	// u,v,w must be linearly dependant!
	matrix A(3, 3); 
	A = u.subs(i == x), u.subs(i == y), u.subs(i == z),
		v.subs(i == x), v.subs(i == y), v.subs(i == z),
		w.subs(i == x), w.subs(i == y), w.subs(i == z);

	cout << "A = " << A << endl;
	cout << endl;

	cout << "|A| = " << eval(determinant(A)) << endl;
	cout << endl;

	cout << "|A|0 = " << eval(determinant(A).subs(t == 0)) << endl;
	cout << "|A|1 = " << eval(determinant(A).subs(t == 1)) << endl;
*/
/*
	ex u = u0 * T*T + u1 * 2 * T*t + u2 * t * t;
	ex v = v0 * T*T*T + v1 * 3 * T*T*t + v2 * 3 * T*t*t + v3 * t*t*t;
	ex w = w0 * T*T*T + w1 * 3 * T*T*t + w2 * 3 * T*t*t + w3 * t*t*t;

	cout << "u = " << u << endl;
	cout << "v = " << v << endl;
	cout << "w = " << w << endl;
	cout << endl;

	// u,v,w must be linearly dependant!
	matrix A(3, 3); 
	A = idx(u, 0, 0), 
		v, w;

	cout << "A = " << A << endl;
	cout << endl;

	cout << "A = " << determinant(A) << endl;
*/

/*
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
*/

	return 0;
}


