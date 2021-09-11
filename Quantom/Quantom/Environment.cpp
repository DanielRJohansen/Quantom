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
		//auto start = std::chrono::high_resolution_clock::now();



		engine->step();
		//exit(1);

		display->render(simulation);

		interface->handleEvents();	
		if (interface->quit)
			display->terminate();



		//auto stop = std::chrono::high_resolution_clock::now();
		//int duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
		
		
		
		
		


	}
}