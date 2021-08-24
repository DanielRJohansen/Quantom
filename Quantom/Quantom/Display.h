#pragma once

#include <SFML/Graphics.hpp>
#include "Raytracer.cuh"

class Display
{
public:
	Display() {}
	Display(Simulation* simulation);
	void render();

	void terminate();
	sf::RenderWindow* window;
	Raytracer* raytracer;

private:
	

};

