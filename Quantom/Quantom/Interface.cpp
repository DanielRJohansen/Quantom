#include "Interface.h"



Interface::Interface(sf::RenderWindow* window) : window(window) {
	
}

void Interface::handleEvents() {


	sf::Event event;
	while (window->pollEvent(event))
	{
		if (event.type == sf::Event::Closed) {

			quit = true;
		}
	}
}