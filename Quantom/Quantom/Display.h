#pragma once

#include <SFML/Graphics.hpp>
#include "Raytracer.cuh"

class Display
{
public:
	Display() {}
	Display(Simulation* simulation);
	void render(Simulation* simulation);
	void animate(Trajectory* traj);

	void terminate();
	sf::RenderWindow* window;
	Raytracer* raytracer;

private:
	void draw(uint8_t* image);
	void sleep(int ms);
	uint8_t* enhance(uint8_t* im, int from_size);	// doubles image resolution

	int xyToIndex(int x, int y, int size_x) {
		return (x + y * size_x) * 4;
	}




	sf::Texture texture;
	sf::Sprite sprite;

};

