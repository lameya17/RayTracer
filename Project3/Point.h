#ifndef POINT_H
#define POINT_H

class Point
{
	private:
		double x;
		double y;
		double z;

	public:
		Point();

		Point(double, double, double);

		double getX(){ return x; }

		double getY(){ return y; }

		double getZ(){ return z; }
};

	Point::Point()
	{
		x = 0;
		y = 0;
		z = 0;
	}

	Point::Point(double i, double j, double k)
	{
		x = i;
		y = j;
		z = k;
	}

#endif