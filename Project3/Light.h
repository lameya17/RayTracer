#ifndef LIGHT_H
#define LIGHT_H

#include "Vector.h"
#include "Color.h"

class Light
{
	private:
	
	Vector position;
	Color color;

	public:
	
	Light();
	Light(Vector, Color);
	Vector getPosition();
	Color getColor();

};

Light::Light()
{
	position = Vector(0.0, 0.0, 0.0);
	color = Color(0.0f, 0.0f, 0.0f);
}

Vector Light::getPosition()
{
	return position;
}

Color Light::getColor()
{
	return color;
}

Light::Light(Vector pos, Color c)
{
	position = pos;
	color = c;
}

#endif