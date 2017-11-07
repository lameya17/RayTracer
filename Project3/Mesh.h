#ifndef MESH_H
#define MESH_H

#include "Vector.h"

class Mesh
{
public:
	Vector V;
	

	Mesh(double a, double b, double c)
	{
		V = Vector(a, b, c);
	}
};

#endif