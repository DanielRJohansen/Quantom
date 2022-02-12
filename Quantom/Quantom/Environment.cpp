#include "Environment.h"



Environment::Environment() {
	simulation = new Simulation();
	if (!verifySimulationParameters()) {
		exit(0);
	}

	//engine = new Engine;

	
	//Compound mol_dpc = compoundbuilder.buildMolecule("D:\\Quantom\\m40.pdb", "D:\\Quantom\\dpc.itp", 1);
	//Compound mol_4pcw10 = compoundbuilder.buildMolecule("D:\\Quantom\\filaggrin\\4pcw_first10.pdb", "D:\\Quantom\\filaggrin\\topol.top", 7);
	//Compound mol_4pcw10 = compoundbuilder.buildMolecule("D:\\Quantom\\filaggrin\\conf.gro", "D:\\Quantom\\filaggrin\\topol.top", 4);
	Compound mol_6lzm_10 = compoundbuilder.buildMolecule("D:\\Quantom\\t4lys\\conf.gro", "D:\\Quantom\\t4lys\\topol.top", 10);

	boxbuilder.buildBox(simulation);
	boxbuilder.addSingleMolecule(simulation, &mol_6lzm_10);
	//boxbuilder.addScatteredMolecules(simulation, &mol_dpc, N_LIPID_COPIES);
	boxbuilder.finishBox(simulation);



	simulation->moveToDevice();	// Only moves the Box to the device

	engine = new Engine(simulation);



	display = new Display(simulation);
	interface = new Interface(display->window);
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during Display Initiation\n");
		exit(1);
	}

}

bool Environment::verifySimulationParameters() {	// Not yet implemented
	assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES);
	assert(THREADS_PER_SOLVENTBLOCK >= N_SOLVATE_MOLECULES);
	assert(BOX_LEN > 3.f);
	assert(BOX_LEN >= CUTOFF + 0.5f);
	assert(simulation->n_compounds <= 1);	// Otherwise data_GAN goes haywire
	return true;
}






void Environment::run() {
	printf("Simulation started\n\n");
	time0 = std::chrono::high_resolution_clock::now();

	

	while (display->window->isOpen()) {

		
		engine->hostMaster();
		engine->deviceMaster();
		
		
		


		handleStatus(simulation);
		handleDisplay(simulation);



		if (handleTermination(simulation)) {
			break;
		}
	}
	printf("\n\n\n########################## SIMULATION FINISHED ##########################\n\n\n\n");
	postRunEvents();
}

void Environment::postRunEvents() {
	printFloat3Matrix(simulation->box->data_GAN, Int3(simulation->getStep(), MAX_COMPOUND_PARTICLES, 6), "D:\\Quantom\\Training\\sim_out.csv");
	//printTrajectory(simulation);
	//printWaterforce(simulation);
	//analyzer.analyzeEnergy(simulation);

	//printOut(simulation->box->outdata, simulation->n_steps);
	//printDataBuffer(simulation->box);
}

void Environment::handleStatus(Simulation* simulation) {
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		printf("\r\tStep #%06d", simulation->box->step);
		double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count();
		int remaining_seconds = (int)(1.f / 1000 * duration / simulation->steps_per_render * (simulation->n_steps - simulation->box->step));
		printf("\tAvg. step time: %.1fms (%05d/%05d/%05d) \tRemaining: %04ds", duration / simulation->steps_per_render, engine->timings.x / simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z, remaining_seconds);
		engine->timings = Int3(0, 0, 0);

		time0 = std::chrono::high_resolution_clock::now();
	}
}

void Environment::handleDisplay(Simulation* simulation) {
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		display->render(simulation);

		interface->handleEvents();
		if (interface->quit)
			display->terminate();
	}
}

bool Environment::handleTermination(Simulation* simulation)
{
	if (simulation->finished)
		return true;
	if (simulation->getStep() >= simulation->n_steps)
		return true;
	if (simulation->box->critical_error_encountered)
		return true;

	return false;
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

void Environment::printFloat3Matrix(Float3* data_matrix, Int3 dim, string filename)
{
	printf("Printing data to file \n");
	cout << filename << endl;

	std::ofstream myfile(filename);

	int n_points = simulation->box->total_particles * simulation->n_steps;
	Float3* traj_host = new Float3[n_points];
	cudaMemcpy(traj_host, simulation->box->trajectory, sizeof(Float3) * n_points, cudaMemcpyDeviceToHost);


	for (int i = 0; i < dim.x; i++) {
		if (i % 1000 == 999)
			printf("Exporting line %d\r", i);


		for (int ii = 0; ii < dim.y; ii++) {
			for (int iii = 0; iii < dim.z; iii++) {
				int index = i * dim.y * dim.z + ii * dim.z + iii;
				Float3 data = data_matrix[index];
				for (int iiii = 0; iiii < 3; iiii++)
					myfile << data.at(iiii) << ';';
			}
		}
		myfile << "\n";
	}
	printf("\n");
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