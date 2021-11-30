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






	sf::Texture texture;
	sf::Sprite sprite;

};

