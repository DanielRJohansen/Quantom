#include "Environment.h"



Environment::Environment() {
	MoleculeLibrary* mol_library = new MoleculeLibrary;

	simulation = new Simulation(mol_library);



	engine = new Engine;
	simulation = engine->prepSimulation(simulation);
	if (!verifySimulationParameters()) {
		exit(0);
	}


	display = new Display(simulation);
	interface = new Interface(display->window);


	
}

bool Environment::verifySimulationParameters() {
	if (simulation->box->blocks_per_dim < 3) {
		printf("\n\n\t\tError: Simulation does not contain a minimum of 27 blocks required for PBC\n");
		return false;
	}
	return true;
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
			//exit(1);
			if (simulation->finished)
				break;
		}
		float duration = (float) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
		int remaining_seconds = (int) (1.f/1000 * duration/simulation->steps_per_render * (simulation->n_steps - simulation->step));
		printf("\tAverage step time: %.1fms (%d/%d/%d)  Remaining: %d s\t",  duration/simulation->steps_per_render, engine->timings.x/simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z / simulation->steps_per_render, remaining_seconds);
		engine->timings = Int3(0, 0, 0);

		if (!(simulation->step % simulation->steps_per_render)) {
			display->render(simulation);

			interface->handleEvents();
			if (interface->quit)
				display->terminate();
		}
		

		if (steps++ == -1)
			break;

		if (simulation->finished)
			break;
		
		if (simulation->step >= simulation->n_steps)
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