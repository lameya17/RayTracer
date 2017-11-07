#include "Camera.h"
#include "Sphere.h"
#include "Triangle.h"
#include <iostream>
#include <conio.h>
#include <stdio.h>
#include "Light.h"


#include<cstdlib>
#include<ctime>
#include "Mesh.h"

double screenHeight;
double screenWidth;
double La = 0.9;
double air = 1.0;
Color background = Color(0.27f, 0.54f, 0.83f);
Object *obj[4];
int objCount = 4;
int fileSize = (int)screenHeight*(int)screenWidth;
Color *pixels = new Color[fileSize];
//double *Lum = new double[fileSize];
Light *yagami[2];
int lightCount = 1;
bool insideSphere = false;
double r_index;

double randGen()
{
	double hi = -2.1;
	double lo = -3.0;
	double randNum = 0.0;

	randNum = (double)rand() / (double)RAND_MAX;
	return lo + randNum * (hi - lo);

	return randNum;
}

void savebmp(const char *filename, int w, int h, int dpi, Color *data)
{
	double sf;
	double R; 
	double G; 
	double B;
	double delta = 0.001;
	double Lwa = 0.0;
	double Ldmax = 100;
	double Lmax = 100;
	int whichTR = 0;		// 0 for ward, 1 for reinhard, 2 for adaptive logarithmic
	FILE *f;
	int k = w*h;
	int s = 4 * k;
	int filesize = 54 + s;

	//const int size = 640 * 480;
	double L;

	double factor = 39.375;
	int m = static_cast<int>(factor);

	int ppm = dpi*m;

	unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
	unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };

	bmpfileheader[2] = (unsigned char)(filesize);
	bmpfileheader[3] = (unsigned char)(filesize >> 8);
	bmpfileheader[4] = (unsigned char)(filesize >> 16);
	bmpfileheader[5] = (unsigned char)(filesize >> 24);

	bmpinfoheader[4] = (unsigned char)(w);
	bmpinfoheader[5] = (unsigned char)(w >> 8);
	bmpinfoheader[6] = (unsigned char)(w >> 16);
	bmpinfoheader[7] = (unsigned char)(w >> 24);

	bmpinfoheader[8] = (unsigned char)(h);
	bmpinfoheader[9] = (unsigned char)(h >> 8);
	bmpinfoheader[10] = (unsigned char)(h >> 16);
	bmpinfoheader[11] = (unsigned char)(h >> 24);

	bmpinfoheader[21] = (unsigned char)(s);
	bmpinfoheader[22] = (unsigned char)(s >> 8);
	bmpinfoheader[23] = (unsigned char)(s >> 16);
	bmpinfoheader[24] = (unsigned char)(s >> 24);

	bmpinfoheader[25] = (unsigned char)(ppm);
	bmpinfoheader[26] = (unsigned char)(ppm >> 8);
	bmpinfoheader[27] = (unsigned char)(ppm >> 16);
	bmpinfoheader[28] = (unsigned char)(ppm >> 24);

	bmpinfoheader[29] = (unsigned char)(ppm);
	bmpinfoheader[30] = (unsigned char)(ppm >> 8);
	bmpinfoheader[31] = (unsigned char)(ppm >> 16);
	bmpinfoheader[32] = (unsigned char)(ppm >> 24);

	f = fopen(filename, "wb");

	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);

	for (int i = 0; i < k; i++)
	{
		Color rgb = data[i];
		double temp = data[i].getR();

		if (temp < data[i].getG())
			temp = data[i].getG();
		if (temp < data[i].getB())
			temp = data[i].getB();
		double maximum = temp;

		/*if (maximum > 1.0)
		{
			data[i].setR(data[i].getR() / maximum);
			data[i].setG(data[i].getG() / maximum);
			data[i].setB(data[i].getB() / maximum);
		}*/


		R = data[i].getR() * Lmax;
		G = data[i].getG() * Lmax;
		B = data[i].getB() * Lmax;

		//if (i = 1579)
	//		std::cout << "pizdec";
		L = 0.27*R + 0.67*G + 0.06*B;
		
		Lwa = Lwa + log(delta + L);
	}

	Lwa = exp(Lwa/k);
	sf = pow((1.219 + pow(Ldmax / 2, 0.4)) / (1.219 + pow(Lwa, 0.4)), 2.5) / Ldmax;

	for (int i = 0; i < k; i++)
	{
		// Overall luminance
		//L(x, y) = 0.27R(x, y) + 0.67G(x, y) + 0.06B(x, y)
		R = data[i].getR() * Lmax;
		G = data[i].getG() * Lmax;
		B = data[i].getB() * Lmax;

		
		if (whichTR == 0)			// Ward
		{		
			double red = sf * R;
			double green = sf * G;
			double blue = sf * B;

			double t = red;

			if (t < green)
				t = green;
			if (t < blue)
				t = blue;
			double max = t;

			if (max > 1.0)
			{

				red = red / max;
				green = green / max;
				blue = blue / max;
			}
			
			red = red * 255;
			blue = blue * 255;
			green = green * 255;

			
			
			unsigned char color[3] = { (int)floor(blue), (int)floor(green), (int)floor(red) };

			fwrite(color, 1, 3, f);

		}
		else if (whichTR == 1)		// Reinhard
		{
			double a = 0.18;
			double Rs = R * a / Lwa;
			double Gs = G * a / Lwa;
			double Bs = B * a / Lwa;

			double Rr = Rs / (1 + Rs);
			double Gr = Gs / (1 + Gs);
			double Br = Bs / (1 + Bs);



			double red = Rr * Ldmax;
			double green = Gr * Ldmax;
			double blue = Br * Ldmax;

			//double red = Rr * 255;
			//double green = Gr * 255;
			//double blue = Br * 255;
			double t = red;

			if (t < green)
				t = green;
			if (t < blue)
				t = blue;
			double max = t;

			if (max > 1.0)
			{

				red = red / max;
				green = green / max;
				blue = blue / max;
			}

			red = red * 255;
			blue = blue * 255;
			green = green * 255;

			unsigned char color[3] = { (int)floor(blue), (int)floor(green), (int)floor(red) };

			fwrite(color, 1, 3, f);
		}
		else if (whichTR == 2)		// logarithmic(advanced)
		{
			L = 0.2999 * R + 0.587 * G + 0.114 * B;
			double Lwmax = Lmax / Lwa;
			L = L / Lwa;
			double b = 0.85;
			double biasexpo = log(b) / log(0.5);
			double term2 = pow(L / Lwmax, biasexpo) * 0.8;
			double term1 = 1 / log10(Lwmax + 1);
			double Ld = Ldmax * 0.01 * log(L + 1) / log(2 + term2);

			double red = Ld * R;
			double green = Ld * G;
			double blue = Ld * B;

			double t = red;

			if (t < green)
				t = green;
			if (t < blue)
				t = blue;
			double max = t;

			if (max > 1.0)
			{

				red = red / max;
				green = green / max;
				blue = blue / max;
			}

			red = red * 255;
			blue = blue * 255;
			green = green * 255;



			unsigned char color[3] = { (int)floor(blue), (int)floor(green), (int)floor(red) };

			fwrite(color, 1, 3, f);
		}
	}

	fclose(f);
}


Ray TransmitRay(Ray revI, Vector Normal, bool insideSphere, Vector insp, double ni, double nr)
{
	double ri = ni/nr;
	//Vector revdir = Vector(-1*incident.getDirection().getVectX(), -1*incident.getDirection().getVectY(), -1*incident.getDirection().getVectZ());
	//Ray revI = Ray(insp, revdir, incident.getRI());

	double temp1 = 1 - (ri*ri)*(1 - pow(Normal.dotProduct(revI.getDirection()), 2));
	double temp2 = ri * (Normal.dotProduct(revI.getDirection())) - sqrt(temp1);

	Vector v1 = Vector(temp2*Normal.getVectX(), temp2*Normal.getVectY(), temp2*Normal.getVectZ());
	Vector v2 = Vector(ri*revI.getDirection().getVectX(), ri*revI.getDirection().getVectY(), ri*revI.getDirection().getVectZ());

	if (insideSphere)
		ri = air;
	else
		ri = ni;
	return Ray(insp, v1.subtract(v2), ri);
}


//Color background = Color(0.0f, 0.0f, 0.0f);



Color illuminate(Ray r, int depth)
{
	//r_index = r.getRI();
	Ray transmit;
	Ray reflect;
	Color retcolor;
	double lowest = 1000000;
	int lowestIndex = -1;
	//int count = 0;
	for (int count = 0; count < objCount; count++)
	{
		double insp = obj[count]->Intersection(r);
		if (insp == -1)
		{
			continue;
		}
		else
		{
			if (lowest > insp)
			{
				lowestIndex = count;
				lowest = insp;
			}
		}
	}
	if (lowestIndex == -1)
	{
		//set background color
		retcolor.setR(background.getR());
		retcolor.setG(background.getG());
		retcolor.setB(background.getB());

		return retcolor;
	}
	else
	{
		bool tileFlag = false;
		double shadowScale = 1;
		bool shadows = false;
		Color diffuse = Color(0.0f, 0.0f, 0.0f);
		Color ambient = Color(0.0f, 0.0f, 0.0f);
		Color specular = Color(0.0f, 0.0f, 0.0f);
		Vector intersection = Vector(r.getOrigin().getVectX() + r.getDirection().getVectX() * lowest,
			r.getOrigin().getVectY() + r.getDirection().getVectY() * lowest,
			r.getOrigin().getVectZ() + r.getDirection().getVectZ() * lowest);

		for (int l = 0; l < lightCount; l++)
		{
			
			shadows = false;
			Vector N1 = obj[lowestIndex]->getNormalAt(intersection);
			Vector N = obj[lowestIndex]->getNormalAt(intersection);
			N = N.normalize();

			//Vector RayDir = intersection.subtract(r.getOrigin());
			Ray shadowRay = Ray(intersection, yagami[l]->getPosition().subtract(intersection).normalize(), obj[lowestIndex]->transI);

			
			Vector V = (intersection.subtract(yagami[l]->getPosition())).normalize();

			// S - 2 S.n/|n^2| n
			Vector R = shadowRay.getDirection().subtract(
				Vector(2 * N.dotProduct(shadowRay.getDirection()) * N.getVectX(),
				2 * N.dotProduct(shadowRay.getDirection()) * N.getVectY(),
				2 * N.dotProduct(shadowRay.getDirection()) * N.getVectZ()));


			double temp = 2 * (r.getDirection().dotProduct(N));
			R = R.normalize();
			// double temp = 2 * (r.getDirection().dotProduct(N));
			// reflect  = incident ray - (normal * temp);
			Vector ref = r.getDirection().subtract(Vector(N.getVectX() * temp, N.getVectY() * temp, N.getVectZ() * temp));
			reflect = Ray(intersection, ref, r.getRI());

			Vector revdir = Vector(-1*r.getDirection().getVectX(), -1*r.getDirection().getVectY(), -1*r.getDirection().getVectZ());
			Ray revI = Ray(intersection, revdir, r.getRI());
			Vector Ni = N;
			double ratio = revI.getRI() / obj[lowestIndex]->transI;;
			double ni = revI.getRI();
			double nr = obj[lowestIndex]->transI;

			if (revI.getDirection().dotProduct(N) > 0)
			{
				insideSphere = false;
				ratio = revI.getRI() / obj[lowestIndex]->transI;
			}
			else if (lowestIndex == 0)
			{
				insideSphere = true;
				Ni = Vector(-1 * N1.getVectX(), -1 * N1.getVectY(), -1 * N1.getVectZ());
				ratio = obj[lowestIndex]->transI/air;
				ni = obj[lowestIndex]->transI;
				nr = air;
			}
			N1 = N1.normalize();
			Ni = Ni.normalize();
			

			double temp1 = 1 - (ratio*ratio)*(1 - pow(Ni.dotProduct(r.getDirection()), 2));
			if (temp1 < 0)
				transmit = reflect;
			else
				transmit = TransmitRay(revI, Ni, insideSphere, intersection, ni, nr);

			if (insideSphere == true)
			{
			//	insideSphere = false;
			}
			
			// shadows
			for (int i = 0; i < objCount; i++)
			{
				if (obj[i]->Intersection(shadowRay) > 0.0001)
				{
					if (obj[i]->Kt > 0)
					{
						shadowScale = shadowScale - obj[i]->Kt;
					}
					
					shadowScale += 1;
					shadows = true;
					break;
				}
			}



			if (lowestIndex > 1)
			{
				int tile = (int)floor(intersection.getVectX() / 1.35) + (int)floor(intersection.getVectZ() / 1.35);

				if (tile % 2 == 0)
				{
					tileFlag = false;
					if (N.dotProduct(shadowRay.getDirection()) > 0.0001)
					{

						diffuse.setR(diffuse.getR() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getR() * obj[lowestIndex]->getColor(1).getR()));
						diffuse.setG(diffuse.getG() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getG() * obj[lowestIndex]->getColor(1).getG()));
						diffuse.setB(diffuse.getB() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getB() * obj[lowestIndex]->getColor(1).getB()));

					}
					if (V.dotProduct(R) > 0.0001)
					{
						//Specular = Ks * (DOT(V, 2 * DOT(N, L) * N - L))^roughness * Od * Ld

						specular.setR(specular.getR() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getR())));
						specular.setG(specular.getG() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getG())));
						specular.setB(specular.getB() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getB())));

					}
					continue;
				}
				else
				{
					tileFlag = true;
					if (N.dotProduct(shadowRay.getDirection()) > 0.0001)
					{

						diffuse.setR(diffuse.getR() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getR() * obj[lowestIndex]->getColor(2).getR()));
						diffuse.setG(diffuse.getG() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getG() * obj[lowestIndex]->getColor(2).getG()));
						diffuse.setB(diffuse.getB() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getB() * obj[lowestIndex]->getColor(2).getB()));

					}
					if (V.dotProduct(R) > 0.0001)
					{
						//Specular = Ks * (DOT(V, 2 * DOT(N, L) * N - L))^roughness * Od * Ld

						specular.setR(specular.getR() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getR())));
						specular.setG(specular.getG() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getG())));
						specular.setB(specular.getB() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getB())));

					}
					continue;
				}
			}

			if (N.dotProduct(shadowRay.getDirection()) > 0.0001)
			{

				diffuse.setR(diffuse.getR() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getR() * obj[lowestIndex]->getColor(1).getR()));
				diffuse.setG(diffuse.getG() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getG() * obj[lowestIndex]->getColor(1).getG()));
				diffuse.setB(diffuse.getB() + (N.dotProduct(shadowRay.getDirection())*yagami[l]->getColor().getB() * obj[lowestIndex]->getColor(1).getB()));
			}


			if (V.dotProduct(R) > 0.0001)
			{
				specular.setR(specular.getR() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getR())));
				specular.setG(specular.getG() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getG())));
				specular.setB(specular.getB() + (pow(V.dotProduct(R), obj[lowestIndex]->Ke) * (yagami[l]->getColor().getB())));
			}

		}
		if (tileFlag == false)
		{
			ambient.setR(obj[lowestIndex]->Ka * obj[lowestIndex]->getColor(1).getR() * La); // shadowScale);
			ambient.setG(obj[lowestIndex]->Ka * obj[lowestIndex]->getColor(1).getG() * La); // shadowScale);
			ambient.setB(obj[lowestIndex]->Ka * obj[lowestIndex]->getColor(1).getB() * La); // shadowScale);

		}
		else
		{

			ambient.setR(obj[lowestIndex]->Ka * obj[lowestIndex]->getColor(2).getR() * La / shadowScale);
			ambient.setG(obj[lowestIndex]->Ka * obj[lowestIndex]->getColor(2).getG() * La / shadowScale);
			ambient.setB(obj[lowestIndex]->Ka * obj[lowestIndex]->getColor(2).getB() * La / shadowScale);

		}
		// final specular claculation
		specular.setR(obj[lowestIndex]->Ks * specular.getR());
		specular.setG(obj[lowestIndex]->Ks * specular.getG());
		specular.setB(obj[lowestIndex]->Ks * specular.getB());

		// diffuse final calculation
		diffuse.setR(obj[lowestIndex]->Kd * diffuse.getR());
		diffuse.setG(obj[lowestIndex]->Kd * diffuse.getG());
		diffuse.setB(obj[lowestIndex]->Kd * diffuse.getB());
		

		retcolor.setR((diffuse.getR() + ambient.getR() + specular.getR()) / shadowScale);
		retcolor.setG((diffuse.getG() + ambient.getG() + specular.getG()) / shadowScale);
		retcolor.setB((diffuse.getB() + ambient.getB() + specular.getB()) / shadowScale);
		if (insideSphere)
		{
			//retcolor.setR(0.0f);
			//retcolor.setG(0.0f);
			//retcolor.setB(0.0f);
		}

		if (insideSphere)
			insideSphere = false;
		//return retcolor;


		if (depth < 6)
		{
			if (obj[lowestIndex]->Kr > 0)
			{
				retcolor.setR(retcolor.getR() + obj[lowestIndex]->Kr * illuminate(reflect, depth + 1).getR());
				retcolor.setG(retcolor.getG() + obj[lowestIndex]->Kr * illuminate(reflect, depth + 1).getG());
				retcolor.setB(retcolor.getB() + obj[lowestIndex]->Kr * illuminate(reflect, depth + 1).getB());
			}
			if (obj[lowestIndex]->Kt > 0)
			{
				retcolor.setR(retcolor.getR() + obj[lowestIndex]->Kt * illuminate(transmit, depth + 1).getR());
				retcolor.setG(retcolor.getG() + obj[lowestIndex]->Kt * illuminate(transmit, depth + 1).getG());
				retcolor.setB(retcolor.getB() + obj[lowestIndex]->Kt * illuminate(transmit, depth + 1).getB());
			}
		}
		return retcolor;

	}
}


int main()
{
	srand(time(NULL));

	Camera camera = Camera();
	Vector cameraU, cameraV, cameraN;
	Vector center;
	screenWidth = 640.0;
	screenHeight = 480.0;
	double aspectratio = screenWidth / screenHeight;
	yagami[0] = new Light(Vector(2.0, 2.0, 0.0), Color(1.0f, 1.0f, 1.0f));
	yagami[1] = new Light(Vector(0.0, 0.0, 0.0), Color(1.0f, 1.0f, 1.0f));

	/*
	//Mesh *mesh[9];
	int x_rows = 31;
	int z_rows = 31;

	Mesh *mesh[31][31];
	int meshCount = 0;

	double X = -5;
	double xval = 0.5;
	double Z = 2.0;
	double zval = 0.5;
	
	
	
	for(int z = 0; z < z_rows; z++)
	{
		X = -5;
		for(int x = 0; x <x_rows; x++)
		{
			mesh[z][x] = new Mesh(X, randGen(), Z);	
			X += xval;
		}
		Z += zval;
	}
	
	int tmp = 0;
	int x = 0;
	int z = 0;

	
	while (z < z_rows - 1)
	{
		std ::cout << z;
		x = 0;
		while (x < x_rows - 1)
		{
			obj[tmp++] = new Triangle(mesh[z][x]->V, mesh[z][x+1]->V, mesh[z+1][x+1]->V, Illumination(Color(0.2f, 0.3f, 1.0f), Color(0.2f, 0.3f, 1.0f)), 0.6, 0.2, 0.6, 5, 0.0, 0.0);
			obj[tmp++] = new Triangle(mesh[z][x]->V, mesh[z+1][x+1]->V, mesh[z+1][x]->V, Illumination(Color(0.2f, 0.3f, 1.0f), Color(0.2f, 0.3f, 1.0f)), 0.6, 0.2, 0.6, 5, 0.0, 0.0);
			x++;
		}
		z++;
	}
	*/	

	
	
	obj[0] = new Sphere(Vector(-0.8, 0.22, 5.1), 1.4, Color(0.2f, 0.2f, 0.2f), 0.3, 0.5, 0.3, 5, 0.8, 0.0, 1.0);
	obj[1] = new Sphere(Vector( 1.5, -0.20, 8.5), 2.0, Color(0.0f, 0.0f, 0.4f), 0.3, 0.6, 0.3, 5, 0.0, 0.3, 1.0);


	obj[2] = new Triangle(Vector(4.0, -2.69, 2.68), Vector(4.0, -2.69, 16.61), Vector(-4.0, -2.69, 16.61), Illumination(Color(0.4f, 0.0f, 0.0f), Color(0.5f, 0.5f, 0.1f)), 0.6, 0.4, 0.9, 5, 0.0, 0.0, 1.0);
	obj[3] = new Triangle(Vector(4.0, -2.69, 2.68), Vector(-4.0, -2.69, 16.61), Vector(-4.0, -2.69, 2.68), Illumination(Color(0.4f, 0.0f, 0.0f), Color(0.5f, 0.5f, 0.1f)), 0.6, 0.4, 0.9, 5, 0.0, 0.0, 1.0);
	


	cameraN = camera.getCameraPosition().subtract(camera.getCameraLookAt());
	cameraN = cameraN.normalize();
	cameraU = cameraN.crossProduct(camera.getCameraUp());
	cameraU = cameraU.normalize();
	cameraV = cameraU.crossProduct(cameraN);

	double focalLength = 1;

	// + for maya && - for unity
	center = Vector(camera.getCameraPosition().getVectX() - focalLength * cameraN.getVectX(), 
					camera.getCameraPosition().getVectY() - focalLength * cameraN.getVectY(),
					camera.getCameraPosition().getVectZ() - focalLength * cameraN.getVectZ());

	
	double worldWidth = 2;
	double worldHeight =  worldWidth / aspectratio;	// W = H* aspect ratio

	double pixelDx = worldWidth/screenWidth;
	double pixelDy = worldHeight/screenHeight;

	//starting pixel
	double startX = center.getVectX() - (worldWidth * cameraU.getVectX() + (worldHeight * cameraV.getVectX())) / 2.0;
	double startY = center.getVectY() - (worldWidth * cameraU.getVectY() + (worldHeight * cameraV.getVectY())) / 2.0;
	double startZ = center.getVectZ() - (worldWidth * cameraU.getVectZ() + (worldHeight * cameraV.getVectZ())) / 2.0;

	Vector start = Vector(startX, startY, startZ);
	int fileSize = (int)screenHeight*(int)screenWidth;

	Color *pixels = new Color[fileSize];
	int position = 0;
	Color retcolor;

	for (int i = 0; i < screenHeight; i++)
	{
		for (int j = 0; j < screenWidth; j++)
		{
			Vector temp = Vector(start.getVectX() + cameraU.getVectX() * (j + 0.5) * pixelDx + cameraV.getVectX() * (i + 0.5) * pixelDy,
								 start.getVectY() + cameraU.getVectY() * (j + 0.5) * pixelDx + cameraV.getVectY() * (i + 0.5) * pixelDy,
								 start.getVectZ() + cameraU.getVectZ() * (j + 0.5) * pixelDx + cameraV.getVectZ() * (i + 0.5) * pixelDy);
			Vector d = temp.subtract(camera.getCameraPosition());
			d = d.normalize();
			insideSphere = false;
			Ray r(camera.getCameraPosition(),d, air);

			retcolor = illuminate(r, 0);

			
			

					pixels[position].setR(retcolor.getR());
					pixels[position].setG(retcolor.getG());
					pixels[position].setB(retcolor.getB());

					//Lum[position] = 0.0;
			//}

			position++;
		}
	}

	savebmp("scene.bmp", screenWidth, screenHeight, 72, pixels);
	
	return 0;
}