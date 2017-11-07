#ifndef VECTOR_H
#define VECTOR_H
#include <math.h>

class Vector
{

	private: 
	double x;
	double y;
	double z;


	public:
		
		Vector();

		Vector(double, double, double);

		double getVectX(){return x;}

		double getVectY(){return y;}

		double getVectZ(){return z;}


		double dotProduct(Vector v)
		{
			return this->x*v.getVectX() + this->y*v.getVectY() + this->z*v.getVectZ();
		}

		Vector crossProduct(Vector v)
		{
			return Vector(y*v.getVectZ() - z*v.getVectY(), z*v.getVectX() - x*v.getVectZ(), x*v.getVectY() - y*v.getVectX());
		}

		Vector normalize()
		{
			double magnitude = sqrt((x*x)+(y*y)+(z*z));
			return Vector(x / magnitude, y / magnitude, z / magnitude);
		}

		Vector subtract(Vector a)
		{
			return Vector(this->getVectX() - a.getVectX(), this->getVectY() - a.getVectY(), this->getVectZ() - a.getVectZ());
		}

		/*
		void transform()
		{
		
		}
		*/
};



	Vector:: Vector()
	{
		x = 0;
		y = 0;
		z = 0;
	}


	Vector:: Vector(double i, double j, double k)
	{
		x = i;
		y = j;
		z = k;
	}

#endif
