#include "Environment.h"



Environment::Environment() {
	engine = new Engine(&simulation);
	simulation.moveToDevice();	// Must be done before initiating raytracer!
	display = new Display(&simulation);
	interface = new Interface(display->window);

	
	run();
}

void Environment::run() {


	Molecule h2o;
	printf("Simulation started\n");
	while (display->window->isOpen()) {






		display->render(&simulation);

		interface->handleEvents();	
		if (interface->quit)
			display->terminate();
	}
}