#ifndef ILLUMINATION_H
#define ILLUMINATION_H

#include"Color.h"

class Illumination
{
private:
	

public:
	Color c1;
	Color c2;

	Illumination(Color a)
	{
		c1 = a;
	}

	Illumination(Color a, Color b)
	{
		c1 = a;
		c2 = b;
	}

};

#endif