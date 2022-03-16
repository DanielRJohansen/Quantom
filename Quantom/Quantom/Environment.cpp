#include "Environment.h"



Environment::Environment() {

	display = new DisplayV2();


	simulation = new Simulation();
	verifySimulationParameters();




	//engine = new Engine;

	
	//Compound mol_dpc = compoundbuilder.buildMolecule("D:\\Quantom\\m40.pdb", "D:\\Quantom\\dpc.itp", 1);
	//Compound mol_4pcw10 = compoundbuilder.buildMolecule("D:\\Quantom\\filaggrin\\4pcw_first10.pdb", "D:\\Quantom\\filaggrin\\topol.top", 7);
	//Compound mol_4pcw10 = compoundbuilder.buildMolecule("D:\\Quantom\\filaggrin\\conf.gro", "D:\\Quantom\\filaggrin\\topol.top", 4);
	//Compound mol_6lzm_10 = compoundbuilder.buildMolecule(MOL_FOLDER + "conf.gro", MOL_FOLDER + "topol.top", 10);
	
	
	Compound temp = compoundbuilder.buildMolecule(MOL_FOLDER + "conf.gro", MOL_FOLDER + "topol.top", 14);
	//Molecule mol_6lzm_10 = compoundbuilder.buildMolecule(MOL_FOLDER + "conf.gro", MOL_FOLDER + "topol.top", 14, 0);
	Molecule mol_6lzm_10;
	mol_6lzm_10.n_atoms_total = temp.n_particles;
	mol_6lzm_10.compounds[0] = temp;

	//printf("here %d", temp.n_particles);
	
	boxbuilder.buildBox(simulation);
	//boxbuilder.addSingleMolecule(simulation, &mol_6lzm_10);
	boxbuilder.addSingleMolecule(simulation, &mol_6lzm_10);
	//boxbuilder.addScatteredMolecules(simulation, &mol_dpc, N_LIPID_COPIES);
	delete[] mol_6lzm_10.compounds;
	boxbuilder.finishBox(simulation);

	simulation->moveToDevice();	// Only moves the Box to the device
	verifyBox();

	engine = new Engine(simulation);




}

void Environment::verifySimulationParameters() {	// Not yet implemented
	assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES);
	//assert(THREADS_PER_SOLVENTBLOCK >= N_SOLVATE_MOLECULES);
	assert(BOX_LEN > 3.f);
	//assert(BOX_LEN >= CUTOFF + 0.5f);
	assert(simulation->n_compounds <= 1);	// Otherwise data_GAN goes haywire

	assert(simulation->n_steps % STEPS_PER_LOGTRANSFER == 0);
	assert(simulation->n_steps % STEPS_PER_THERMOSTAT == 0);
	assert(simulation->n_steps % STEPS_PER_TRAINDATATRANSFER == 0);

	assert(STEPS_PER_THERMOSTAT % STEPS_PER_LOGTRANSFER == 0);		// Change to trajtransfer later
	//assert(STEPS_PER_THERMOSTAT >= STEPS_PER_LOGTRANSFER);

}

void Environment::verifyBox() {

	printf("FSDA OJØÆMVPSEWM %d\n", simulation->n_compounds);
	for (int c = 0; c < simulation->n_compounds; c++) {
		printf("Compound radius: %f\n", simulation->compounds_host[c].confining_particle_sphere);
		if ((simulation->compounds_host[c].confining_particle_sphere * 1.1) > BOX_LEN_HALF) {
			printf("Compound %d too large for simulation-box\n");
			exit(1);
		}
	}
}







void Environment::run() {
	printf("Simulation started\n\n");

#ifdef ENABLE_CHRONO
	time0 = std::chrono::high_resolution_clock::now();
#endif

	//while (display->window->isOpen()) {
	while (display->checkWindowStatus()) {
		
		engine->deviceMaster();		// Device first, otherwise offloading data always needs the last datapoint!
		engine->hostMaster();

		handleStatus(simulation);
		handleDisplay(simulation);

		if (handleTermination(simulation)) {
			break;
		}
	}
	printf("\n\n\n########################## SIMULATION FINISHED ##########################\n\n\n\n");

	if (simulation->finished || simulation->box->critical_error_encountered) {
		postRunEvents();
	}
}

void Environment::postRunEvents() {
#ifndef __linux__
	simulation->out_dir += "\\Steps_" + to_string(simulation->getStep());
	int check = _mkdir(&(simulation->out_dir[0]));

	Analyzer::AnalyzedPackage analyzed_package = analyzer.analyzeEnergy(simulation);
	dumpToFile(analyzed_package.energy_data, analyzed_package.n_energy_values, simulation->out_dir + "\\energy.bin");
	dumpToFile(analyzed_package.temperature_data, analyzed_package.n_temperature_values, simulation->out_dir + "\\temperature.bin");


	dumpToFile(simulation->traindata_buffer,
		simulation->getStep() * 6 * MAX_COMPOUND_PARTICLES,
		simulation->out_dir + "\\sim_out.bin");

	/*dumpToFile(simulation->box->data_GAN,
		simulation->getStep() * MAX_COMPOUND_PARTICLES * 6,
		simulation->out_dir + "\\sim_out.bin"
	);*/

	string data_processing_command = "C:\\Users\\Daniel\\git_repo\\Quantom\\LIMA_services\\x64\\Debug\\LIMA_services.exe "
		+ simulation->out_dir + " "
		+ to_string(simulation->getStep())
		+ " 0";											// do_shuffle

	cout << data_processing_command << "\n\n";
	system(&data_processing_command[0]);		
#else
	Analyzer::AnalyzedPackage analyzed_package = analyzer.analyzeEnergy(simulation);
	dumpToFile(analyzed_package.energy_data, analyzed_package.n_energy_values, simulation->out_dir + "\\energy.bin");
	dumpToFile(analyzed_package.temperature_data, analyzed_package.n_temperature_values, simulation->out_dir + "\\temperature.bin");
#endif
}

void Environment::handleStatus(Simulation* simulation) {
#ifdef ENABLE_CHRONO
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		printf("\r\tStep #%06d", simulation->box->step);
		double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count();
		int remaining_minutes = (int)(1.f / 1000 * duration / simulation->steps_per_render * (simulation->n_steps - simulation->box->step) / 60);
		printf("\tAvg. step time: %.1fms (%05d/%05d/%05d) \tRemaining: %04d min", duration / simulation->steps_per_render, engine->timings.x / simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z/simulation->steps_per_render, remaining_minutes);
		engine->timings = Int3(0, 0, 0);


		time0 = std::chrono::high_resolution_clock::now();

	}
#else
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		printf("\r\tStep #%06d", simulation->box->step);
	}
#endif
}

void Environment::handleDisplay(Simulation* simulation) {
	//display->render(simulation);
	
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		display->render(simulation);

		//interface->handleEvents();
		//if (interface->quit)
			//display->terminate();
	}
}

bool Environment::handleTermination(Simulation* simulation)
{
	if (simulation->finished)
		return true;
	if (simulation->getStep() >= simulation->n_steps) {
		simulation->finished = true;
		return true;
	}		
	if (simulation->box->critical_error_encountered)
		return true;

	return false;
}

void Environment::renderTrajectory(string trj_path)
{
	/*
	Trajectory* trj = new Trajectory(trj_path);
	for (int i = 0; i < trj->n_particles; i++) {
		trj->particle_type[i] = 0;
	}
	trj->particle_type[0] = 1;

	display->animate(trj);
	*/
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





template <typename T>
void Environment::dumpToFile(T* data, int n_datapoints, string file_path_s) {	
	char* file_path;
	file_path = &file_path_s[0];
	cout << "Printing to file " << file_path << endl;

	FILE* file;
#ifndef __linux__
	fopen_s(&file, file_path, "wb");
#else
	file = fopen(file_path, "wb");
#endif

	fwrite(data, sizeof(T), n_datapoints, file);
	fclose(file);
}

