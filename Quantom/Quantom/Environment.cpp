#include "Environment.h"



Environment::Environment() {
	display = new DisplayV2();


	simulation = new Simulation();
	verifySimulationParameters();


	ForceFieldMaker* forcefieldmaker = new ForceFieldMaker();
	compoundbuilder = new CompoundBuilder(forcefieldmaker);



	Molecule mol_6lzm_10 = compoundbuilder->buildMolecule(MOL_FOLDER + "conf.gro", MOL_FOLDER + "topol.top", MAX_RESIDUES_TO_LOAD, 0);

	printf("bridges after %d\n", mol_6lzm_10.compound_bridge_bundle->n_bridges);
	//printf("here %d", temp.n_particles);
	printf("mol bridge 0 particles %d\n", mol_6lzm_10.compound_bridge_bundle->compound_bridges[0].n_particles);
	boxbuilder.buildBox(simulation);
	//boxbuilder.addSingleMolecule(simulation, &mol_6lzm_10);
	boxbuilder.addSingleMolecule(simulation, &mol_6lzm_10);
	boxbuilder.solvateBox(simulation);	// Always do after placing compounds

	//boxbuilder.addScatteredMolecules(simulation, &mol_dpc, N_LIPID_COPIES);
	delete[] mol_6lzm_10.compounds;
	boxbuilder.finishBox(simulation);


	for (int c = 0; c < simulation->n_compounds; c++) {
		break;
		Compound* compound = &simulation->box->compounds[c];
		for (int i = 0; i < compound->n_particles; i++) {
			printf("%d %d %d:   ", c, i, compound->particle_global_ids[i]);
			compound->prev_positions[i].print();

			if (compound->particle_global_ids[i] == 200) {
				for (int ii = 0; ii < compound->lj_ignore_list[i].n_ignores; ii++) {
					//printf("%d %d %d\n", compound->lj_ignore_list[i].local_ids[ii], compound->lj_ignore_list[i].compound_ids[ii], 0);// compound->lj_ignore_list[i].)
				}
				//exit(0);
			}

		}
	}


	simulation->moveToDevice();	// Only moves the Box to the device
	verifyBox();

	//engine = new Engine(simulation, forcefieldmaker->getForcefield());
	engine = new Engine(simulation, forcefieldmaker->getNBForcefield());



}

void Environment::verifySimulationParameters() {	// Not yet implemented
	assert(THREADS_PER_COMPOUNDBLOCK >= MAX_COMPOUND_PARTICLES);
	//assert(THREADS_PER_SOLVENTBLOCK >= N_SOLVATE_MOLECULES);
	assert(BOX_LEN > 3.f);
	//assert(BOX_LEN >= CUTOFF + 0.5f);
	//assert(simulation->n_compounds <= 1);	// Otherwise data_GAN goes haywire

	assert(simulation->n_steps % STEPS_PER_LOGTRANSFER == 0);
	assert(simulation->n_steps % STEPS_PER_THERMOSTAT == 0);
	assert(simulation->n_steps % STEPS_PER_TRAINDATATRANSFER == 0);

	assert(STEPS_PER_THERMOSTAT % STEPS_PER_LOGTRANSFER == 0);		// Change to trajtransfer later
	//assert(STEPS_PER_THERMOSTAT >= STEPS_PER_LOGTRANSFER);
	assert(THREADS_PER_SOLVENTBLOCK >= MAX_COMPOUND_PARTICLES);


	printf("Simulation parameters verified\n");
}

void Environment::verifyBox() {
	for (int c = 0; c < simulation->n_compounds; c++) {
		printf("Compound radius: %f\t center: %f %f %f\n", simulation->compounds_host[c].confining_particle_sphere, simulation->compounds_host[c].center_of_mass.x, simulation->compounds_host[c].center_of_mass.y, simulation->compounds_host[c].center_of_mass.z);
		if ((simulation->compounds_host[c].confining_particle_sphere * 1.1) > BOX_LEN_HALF) {
			printf("Compound %d too large for simulation-box\n");
			exit(1);
		}
	}
}







void Environment::run() {
	printf("Simulation started\n\n");

	time0 = std::chrono::high_resolution_clock::now();

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

	dumpToFile(simulation->box->data_GAN,
		simulation->getStep() * simulation->n_compounds * MAX_COMPOUND_PARTICLES * 6,
		simulation->out_dir + "\\sim_out.bin"
	);
	return;
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
	if (!(simulation->getStep() % simulation->steps_per_render)) {
		printf("\r\tStep #%06d", simulation->box->step);
		double duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - time0).count();
		int remaining_minutes = (int)(1.f / 1000 * duration / simulation->steps_per_render * (simulation->n_steps - simulation->box->step) / 60);
		printf("\tAvg. step time: %.1fms (%05d/%05d/%05d) \tRemaining: %04d min", duration / simulation->steps_per_render, engine->timings.x / simulation->steps_per_render, engine->timings.y / simulation->steps_per_render, engine->timings.z/simulation->steps_per_render, remaining_minutes);
		engine->timings = Int3(0, 0, 0);


		time0 = std::chrono::high_resolution_clock::now();

	}
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

