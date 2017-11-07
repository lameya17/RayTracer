#ifndef SPHERE_H
#define SPHERE_H
//#include<d:\courses\spring 2015\include\cmatrix>

//typedef techsoft::matrix<double> myMatrix;

#include "Object.h"
#include "Vector.h"
#include "Illumination.h"
#include "math.h"

class Sphere : public Object
{
	private:
		Vector center;
		double radius;
		Color color;
		
	public:
		
		Sphere();
		Sphere(Vector, double, Illumination, double, double, double, double, double, double, double);

		//myMatrix worldToCamera(myMatrix);
		void setTransI(double t)
		{
			transI = t;
		}

		virtual Color getColor(int a)
		{
			return color;
		}

		double getTransI()
		{
			return transI;
		}

		Vector getNormalAt(Vector pos)
		{
			return pos.subtract(this->center);
		}

		double Intersection(Ray r)
		{
			double a = 1;//r.getDirection().dotProduct(r.getDirection());
			double b = 2.0 * r.getDirection().dotProduct(r.getOrigin().subtract(this->center));
			//double b = 2.0 * r.getDirection().getVectX() * (r.getOrigin().getVectX() - this->center.getVectX()) + 
			Vector t = r.getOrigin().subtract(this->center);
			//double c = pow(t.dotProduct(t), 2) - (this->radius * this->radius);
			double c = t.dotProduct(t) - (this->radius * this->radius);
			double disco;
			if (b*b >= 4 * a*c){
				disco = sqrt(b*b - 4.0*a*c);
				//printf("%f",disco);
			}
			else
			{
				//printf("%f", b*b - 4 * a*c);
				disco = -1;
				return -1;
			}
			
			double s1 = (-b + disco) / 2.0*a;
			
			double s2 = (-b - disco) / 2.0*a;

			if (s2 > 0.0001)
				return s2;
			if (s1 > 0.0001)
				return s1;
			//printf("%d", s1);
			//printf("%d", s2);

			//if (s1 < 0.0001 || s2 < 0.0001)
				return -1;
			/*

			if (s1 > s2)
				return s2;
			else
				return s1;
*/
		}
};

/*
myMatrix Sphere::worldToCamera(myMatrix m)
{
	double mat[] = { center.getVectX(), center.getVectY(), center.getVectZ(), 1 };
	myMatrix m1(4, 1, mat);
	myMatrix m2 = m*m1;
	return m2;
}
*/

Sphere::Sphere()
{
	center = Vector(0, 0, 0);
	radius = 0.0;
	color = Color(0.0f, 0.0f, 0.0f);
}

Sphere::Sphere(Vector v, double r, Illumination i, double d, double s, double a, double Kee, double ref, double trans, double ri)
{
	Kd = d;
	Ks = s;
	Ka = a;
	Ke = Kee;
	center = v;
	radius = r;
	color = i.c1;
	Kr = ref;
	Kt = trans;
	transI = ri;
}


#endif