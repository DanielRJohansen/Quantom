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

bool Environment::verifySimulationParameters() {	// Not yet implemented
	return true;
}

void Environment::run() {


	Molecule h2o;
	printf("Simulation started\n\n");
	int steps = 0;


	

	while (display->window->isOpen()) {
		auto t0 = std::chrono::high_resolution_clock::now();


		for (int i = 0; i < simulation->steps_per_render; i++) {
			engine->step();
			printf("\r\tStep #%05d", simulation->step);
			//simulation->finished = true;
			//exit(1);
			if (simulation->finished)
				break;
		}
		float duration = (float) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
		int remaining_seconds = (int) (1.f/1000 * duration/simulation->steps_per_render * (simulation->n_steps - simulation->step));
		printf("\tAverage step time: %.1fms (%d/%d/%d)  Remaining: %ds",  duration/simulation->steps_per_render, engine->timings.x/simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z / simulation->steps_per_render, remaining_seconds);
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


		
		
		
		


	}

	printf("\n\n\n########################## SIMULATION FINISHED ##########################\n");
	printOut(simulation->box->outdata);
}


void Environment::printOut(float* data) {
	std::ofstream myfile("D:\\Quantom\\log.csv");
	for (int i = 0; i < 10; i++) {
		myfile << "Data" << std::to_string(i) << ";";
	}
	myfile << "\n";

	int n_datapoints = 10000;

	for (int i = 0; i < n_datapoints; i++) {
		for (int j = 0; j < 10; j++) {
			myfile << data[i + j * n_datapoints] << ";";
		}
		myfile << "\n";

	}
		
	myfile.close();
}