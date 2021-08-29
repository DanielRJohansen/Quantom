#include "Environment.h"



Environment::Environment() {
	engine = new Engine(&simulation);
	printf("Block 0 bodies: %d\n", simulation.box->blocks[0].n_bodies);

	simulation.moveToDevice();	// Must be done before initiating raytracer!

	printf("Block 0 bodies: %d\n", simulation.box->blocks[0].n_bodies);
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