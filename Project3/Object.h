#ifndef OBJECT_H
#define OBJECT_H

#include"Color.h"
#include"Vector.h"
#include"Ray.h"

class Object
{
		
	public:

		double transI;
		double Kd;
		double Ka;
		double Ks;
		double Ke;

		double Kr;	//reflection coefficient
		double Kt;	//transmission coefficient
		Object();

		virtual Color getColor(int a)
		{
			return Color(0.0, 0.0, 0.0);
		}

		virtual double Intersection(Ray r) = 0;
		virtual Vector getNormalAt(Vector pos) = 0;

};

Object::Object()
{

}
#endif