#ifndef RAY_H
#define RAY_H

#include "Vector.h"

class Ray
{
	private:
		Vector origin;
		Vector direction;
		double ri;

	public:
		Ray();
		Ray(Vector, Vector, double);

		Vector getOrigin(){ return origin;}
		Vector getDirection(){ return direction; }
		double getRI(){ return ri; }

};

Ray::Ray()
{
	origin = Vector(0, 0, 0);
	direction = Vector(1, 0, 0);
}

Ray::Ray(Vector a, Vector b, double ref)
{
	origin = a;
	direction = b;
	ri = ref;
}

#endif;