#include "Environment.h"



Environment::Environment() {
	simulation = new Simulation();

	engine = new Engine;

	
	Compound main_molecule = compoundbuilder.buildMolecule("D:\\Quantom\\m40.pdb", "D:\\Quantom\\dpc.itp");


	simulation = engine->prepSimulation(simulation, &main_molecule);
	printf("Engine ready\n");
	if (!verifySimulationParameters()) {
		exit(0);
	}

	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error before Display Initiation\n");
		exit(1);
	}

	display = new Display(simulation);
	interface = new Interface(display->window);
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during Display Initiation\n");
		exit(1);
	}

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

			if (simulation->finished)
				break;
		}



		printf("\r\tStep #%06d", simulation->box->step);
		double duration = (double) std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count();
		int remaining_seconds = (int) (1.f/1000 * duration/simulation->steps_per_render * (simulation->n_steps - simulation->box->step));
		printf("\tAvg. step time: %.1fms (%d/%d/%d) \tRemaining: %04ds",  duration/simulation->steps_per_render, engine->timings.x/simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z / simulation->steps_per_render, remaining_seconds);
		engine->timings = Int3(0, 0, 0);

		if (!(simulation->box->step % simulation->steps_per_render)) {
			display->render(simulation);

			interface->handleEvents();
			if (interface->quit)
				display->terminate();
		}


		if (steps++ == -1)
			break;

		if (simulation->finished)
			break;
		
		if (simulation->box->step >= simulation->n_steps)
			break;


		
		
		
		


	}

	//exit(1);


	printf("\n\n\n########################## SIMULATION FINISHED ##########################\n\n\n\n");
	printTrajectory(simulation);
	printWaterforce(simulation);
	analyzer.analyzeEnergy(simulation);

	printOut(simulation->box->outdata, simulation->n_steps);
	//printDataBuffer(simulation->box);
}

void Environment::renderTrajectory(string trj_path)
{
	Trajectory* trj = new Trajectory(trj_path);
	for (int i = 0; i < trj->n_particles; i++) {
		trj->particle_type[i] = 0;
	}
	trj->particle_type[0] = 1;

	display->animate(trj);
}

void Environment::makeVirtualTrajectory(string trj_path, string waterforce_path) {
	Trajectory* trj = new Trajectory(trj_path);
	Trajectory* force_buffer = new Trajectory(waterforce_path);
	int n_steps = trj->n_steps;

	printf(" part: %d\n", trj->n_particles);


	Float3* particle_position = new Float3[n_steps];
	for (int step = 0; step < n_steps; step++)
		particle_position[step] = trj->positions[0 + step * trj->n_particles];
	Float3* forces = force_buffer->positions;
	

	VirtualPathMaker VPM;
	Float3* vp_path = VPM.makeVirtualPath(particle_position, forces, n_steps);

	std::ofstream myfile("D:\\Quantom\\virtrj.csv");
	for (int step = 0; step < n_steps; step++) {

		for (int k = 0; k < 3; k++) {
			myfile << particle_position[step].at(k) << ";";
		}
		for (int k = 0; k < 3; k++) {
			myfile << vp_path[step].at(k) << ";";
		}

		myfile << "\n";
	}
	myfile.close();

}


















void Environment::printOut(double* data, int n_steps) {
	std::ofstream myfile("D:\\Quantom\\log.csv");
	for (int i = 0; i < 10; i++) {
		myfile << "Data" << std::to_string(i+1) << ";";
	}
	myfile << "\n";


	for (int i = 0; i < n_steps; i++) {
		for (int j = 0; j < 10; j++) {
			myfile << data[j + i * 10] << ";";
		}
		myfile << "\n";
	}
		
	myfile.close();
}

void Environment::printTrajectory(Simulation* simulation) {
	std::ofstream myfile("D:\\Quantom\\trajectory.csv");

	int n_points = simulation->box->total_particles * simulation->n_steps;
	Float3* traj_host = new Float3[n_points];
	cudaMemcpy(traj_host, simulation->box->trajectory, sizeof(Float3) * n_points, cudaMemcpyDeviceToHost);

	for (int i = 0; i < simulation->box->step; i++) {
		for (int j = 0; j < simulation->box->total_particles; j++) {
			for (int k = 0; k < 3; k++) {
				myfile << traj_host[j + i * simulation->box->total_particles].at(k) << ";";
			}			
		}
		myfile << "\n";
	}
	myfile.close();
}

void Environment::printWaterforce(Simulation* simulation)
{
	std::ofstream myfile("D:\\Quantom\\waterforce.csv");

	for (int i = 0; i < simulation->n_steps; i++) {
		myfile << simulation->box->outdata[7 + i * 10] << ";";
		myfile << simulation->box->outdata[8 + i * 10] << ";";
		myfile << simulation->box->outdata[9 + i * 10] << ";";
		myfile << "\n";
	}

	myfile.close();
}










/*
void Environment::printDataBuffer(Box* box) {
	std::ofstream myfile("D:\\Quantom\\energy.csv");

	double* host_data = engine->getDatabuffer();

	printf("Printing %d columns per step\n", box->n_compounds * 3 * 2);

	for (int i = 0; i < box->step; i++) {
		for (int j = 0; j < box->n_compounds * 3; j++) {
			myfile << host_data[0 + j*2 + i * box->n_compounds * 3 * 2] << ';';
			myfile << host_data[1 + j*2 + i * box->n_compounds * 3 * 2] << ';';
		}
		myfile << "\n";
	}
	myfile.close();

	delete host_data;
}
*/