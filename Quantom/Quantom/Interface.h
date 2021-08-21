#pragma once

#include <SFML/Graphics.hpp>


class Interface
{
public:
	Interface(sf::RenderWindow* window);
	void handleEvents();


	bool quit = false;


private:
	sf::RenderWindow* window;
};

