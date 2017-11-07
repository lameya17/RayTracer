#ifndef TRIANGLE_H
#define TRIANGLE_H
//#include<d:\courses\spring 2015\include\cmatrix>

//typedef techsoft::matrix<double> myMatrix;

#include "Object.h"
#include "Vector.h"
#include "Illumination.h"
#include "math.h"


class Triangle : public Object
{
	private:

		Vector vertex1;
		Vector vertex2;
		Vector vertex3;

		Color color1;
		Color color2;

	public:


		Vector normal;
		
		Triangle(Vector, Vector, Vector, Illumination, double, double, double, double, double, double, double);

		//myMatrix worldToCamera(myMatrix);

		virtual Color getColor(int a)
		{
			if (a == 1)
				return color1;
			else
				return color2;
		}

		Vector getNormalAt(Vector pos)
		{
			return normal;
		}


		double Intersection(Ray r)
		{
			double epsilon = 0.00001;
			Vector edge1;
			Vector edge2;
			Vector P, Q, T;
			float det, inv_det, u, v;
			float t;	// intersection point

			// two edges sharing vertex 1
			edge1 = vertex2.subtract(vertex1);
			edge2 = vertex3.subtract(vertex1);

			// determinant
			P = r.getDirection().crossProduct(edge2);
			det = edge1.dotProduct(P);

			if (det > -epsilon && det < epsilon)	// no intersection
				return -1;

			inv_det = 1.0f / det;
			// distance from vertex 1 to ray origin
			T = r.getOrigin().subtract(vertex1);

			// calculate u
			u = T.dotProduct(P) * inv_det;


			
			// intersection outside traingle
			if (u < 0.0f || u > 1.0f)
				return -1;

			Q = T.crossProduct(edge1);

			// calculating v
			v = (float)r.getDirection().dotProduct(Q) * inv_det;

			//printf("lala");

			if (v < 0.0f || u + v > 1.0f)	// intersection outside traingle
				return -1;
			//printf("lala");

			t = (float)edge2.dotProduct(Q) * inv_det;
			//printf("%d", t);

			//printf("lala");



			if (t > epsilon)
			{
				//printf("%d", t);
				return t;
			}

			return -1;
		}
};

/*
myMatrix Triangle::worldToCamera(myMatrix m)
{
	double mat[] = { center.getVectX(), center.getVectY(), center.getVectZ(), 1 };
	myMatrix m1(4, 1, mat);
	myMatrix m2 = m*m1;
	return m2;
}
*/

Triangle::Triangle(Vector v1, Vector v2, Vector v3, Illumination i, double d, double s, double a, double kee, double ref, double trans, double ri)
{

	Vector edge1 = v2.crossProduct(v1);
	Vector edge2 = v3.crossProduct(v1);
	normal = edge1.crossProduct(edge2);
	normal = normal.normalize();
	vertex1 = v1;
	vertex2 = v2;
	vertex3 = v3;
	color1 = i.c1;
	color2 = i.c2;
	Kd = d;
	Ks = s;
	Ka = a;
	Ke = kee;
	Kr = ref;
	Kt = trans;
	transI = ri;
}

#endif