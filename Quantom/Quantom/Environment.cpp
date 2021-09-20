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
	int steps = 0;
	engine->countBodies();


	

	while (display->window->isOpen()) {
		auto t0 = std::chrono::high_resolution_clock::now();


		for (int i = 0; i < simulation->steps_per_render; i++) {
			engine->step();
			printf("\r\tStep #%05d", simulation->step);
		}
		float duration = (float) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
		printf("\tAverage step time: %.1f ms.\t",  duration/simulation->steps_per_render);


		if (!(simulation->step % 50)) {



			display->render(simulation);

			interface->handleEvents();
			if (interface->quit)
				display->terminate();
		}
		

		if (steps++ == -1)
			break;

		if (simulation->finished)
			break;
		
		if (simulation->step == simulation->n_steps)
			break;

		//if (simulation->step == 650)
			//break;

		/*
		int duration;
		do {
			auto stop = std::chrono::high_resolution_clock::now();
			duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
		} while (duration < 50);
		*/
		
		
		
		
		


	}

	printf("\n\n\n########################## SIMULATION FINISHED ##########################\n");
	engine->countBodies();
}