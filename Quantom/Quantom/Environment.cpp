#include "Environment.h"



Environment::Environment() {
	display = new Display();
	interface = new Interface(display->window);

	run();
}

void Environment::run() {


	Molecule h2o;

	while (display->window->isOpen()) {






		display->render();

		interface->handleEvents();	
		if (interface->quit)
			display->terminate();
	}
}