#include "Environment.h"



Environment::Environment() {
	MoleculeLibrary* mol_library = new MoleculeLibrary;

	simulation = new Simulation(mol_library);

	engine = new Engine;
	simulation = engine->prepSimulation(simulation);
	

	display = new Display(simulation);
	interface = new Interface(display->window);


	
}

void Environment::run() {


	Molecule h2o;
	printf("Simulation started\n\n");
	while (display->window->isOpen()) {



		engine->step();


		display->render(simulation);

		interface->handleEvents();	
		if (interface->quit)
			display->terminate();
	}
}