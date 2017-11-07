#ifndef COLOR_H
#define COLOR_H

class Color
{
	private:
		float r;
		float g;
		float b;

	public:
		
		Color();
		Color(float, float, float);
		float getR(){ return r; }
		float getG(){ return g; }
		float getB(){ return b; }
		
		void setR(float red){ this->r = red; }
		void setG(float green){ this->g = green; }
		void setB(float blue){ this->b = blue; }
		/*
		Color add(Color c)
		{
			Color a;
			a.setR(this->r + c->r);
			a.setG(this->g + c->g);
			a.setB(this->b + c->b);
		}
		*/
};

Color::Color()
{
	r = 0.0f;
	g = 0.0f;
	b = 0.0f;
}

Color::Color(float i, float j, float k)
{
	r = i;
	g = j;
	b = k;
}


#endif