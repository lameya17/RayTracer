#ifndef CAMERA_H
#define CAMERA_H

#include "Vector.h"

class Camera
{
	private :

		Vector position;
		Vector lookat;
		Vector Up;
	
	public :

		Camera();
		Camera(Vector, Vector, Vector);

		Vector getCameraPosition() { return position; }
		Vector getCameraLookAt() { return lookat; }
		Vector getCameraUp() { return Up; }

};

Camera::Camera()
{
	position = Vector(0, 0, 0);
	lookat = Vector(0, 0, 1);
	Up = Vector(0, 1, 0);
}

#endif