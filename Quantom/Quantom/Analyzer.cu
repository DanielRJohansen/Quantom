#include "Analyzer.cuh"






template<typename T>
void __device__ distributedSummation(T* arrayptr, int array_len) {				// Places the result at pos 0 of input_array
	T temp;			// This is a lazy soluation, but maybe it is also fast? Definitely simple..
	for (int i = 1; i < array_len; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
		if ((threadIdx.x + i) < array_len) {
			temp = arrayptr[threadIdx.x] + arrayptr[threadIdx.x + i];// *0.5f;	// easier to just divide by sum of solvents at host
		}
		__syncthreads();
		arrayptr[threadIdx.x] = temp;
		__syncthreads();
	}
}

double __device__ calcKineticEnergy(Float3* pos1, Float3* pos2, double mass, double dt) {
	LIMAENG::applyHyperpos(pos1, pos2);


	double vel = (*pos1 - *pos2).len() * 0.5f / dt;
	double kinE = 0.5 * mass * vel * vel;
	return kinE;
}

void __global__ monitorCompoundEnergyKernel(Box* box, Float3* traj_buffer, double* potE_buffer, Float3* data_out) {		// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
	__shared__ Float3 energy[MAX_COMPOUND_PARTICLES];
	__shared__ Compound compound;


	const int step = blockIdx.x + 1;
	const int compound_index = blockIdx.y;
	energy[threadIdx.x] = Float3(0.f);

	if (threadIdx.x == 0) {
		//printf("index: %d\n", compound_index + (step - 1) * N_MONITORBLOCKS_PER_STEP);
		data_out[compound_index + (step - 1) * box->n_compounds] = Float3(0, 0, 0);
		//mass = box->compounds[compound_index].particles[0]
		compound = box->compounds[compound_index];
	}

	__syncthreads();

	if (threadIdx.x >= compound.n_particles) {
		return;
	}
	__syncthreads();


	double potE = potE_buffer[threadIdx.x + compound_index * MAX_COMPOUND_PARTICLES + step * box->total_particles_upperbound];

	Float3 pos_tsub1 = traj_buffer[threadIdx.x + compound_index * MAX_COMPOUND_PARTICLES + (step - 1) * box->total_particles_upperbound];
	Float3 pos_tadd1 = traj_buffer[threadIdx.x + compound_index * MAX_COMPOUND_PARTICLES + (step + 1) * box->total_particles_upperbound];
	//LIMAENG::applyHyperpos(&pos_tadd1, &pos_tsub1);
	//testspace::testerfn(4);
	double kinE = calcKineticEnergy(&pos_tadd1, &pos_tsub1, compound.particles[threadIdx.x].mass, box->dt);
	/*
	applyHyperposA(&pos_tadd1, &pos_tsub1);
	double vel = (pos_tadd1 - pos_tsub1).len() * 0.5f / box->dt;
	double kinE = 0.5 * compound.particles[threadIdx.x].mass * vel * vel;
	*/


	double totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();


	/*if (energy[threadIdx.x].len() > 40000) {
		printf("Step %d\n", step);
		energy[threadIdx.x].print('c');
	}*/

	distributedSummation(energy, MAX_COMPOUND_PARTICLES);
	
	__syncthreads();
	if (threadIdx.x == 0) {
		data_out[compound_index + (step - 1) * box->n_compounds] = energy[0];
	}
}





void __global__ monitorSolventEnergyKernel(Box* box, Float3* traj_buffer, double* potE_buffer, Float3* data_out) {
	__shared__ Float3 energy[256];

	//int solvent_index = threadIdx.x + blockIdx.y * THREADS_PER_MONITORBLOCK;
	int solvent_index = threadIdx.x;
	int step = blockIdx.x + 1;
	int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;


	energy[threadIdx.x] = Float3(0.f);

	if (threadIdx.x == 0) {
		//data_out[blockIdx.y + (step - 1) * N_MONITORBLOCKS_PER_STEP] = Float3(0, 0, 0);
	}
	if (solvent_index >= box->n_solvents) {	// Shouldn't be necessary right now..
		return;
	}



	Float3 pos_tsub1 = traj_buffer[compounds_offset + solvent_index + (step - 1) * box->total_particles_upperbound];
	Float3 pos_tadd1 = traj_buffer[compounds_offset + solvent_index + (step + 1) * box->total_particles_upperbound];
	double kinE = calcKineticEnergy(&pos_tadd1, &pos_tsub1, SOLVENT_MASS, box->dt);

	double potE = potE_buffer[compounds_offset + solvent_index + step * box->total_particles_upperbound];

	if (potE > 200'000) {
		//printf("step %04d solvate %04d pot %f\n", step, solvent_index, potE);
	}

	double totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();


		

	distributedSummation(energy, 256);
	if (threadIdx.x == 0) {
		//energy[0].print('b');
		data_out[step-1] = energy[0];
		//data_out[blockIdx.y + (step - 1) * 256] = Float3(1,2,3);
		//data_out[blockIdx.y + step * 256].print('g');
	}	
}



void Analyzer::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol // calculate energies separately for compounds and solvents. weigh averages based on amount of each
	printf("Get step %d\n", simulation->getStep());
	int analysable_steps = simulation->getStep() - 3;
	if (analysable_steps < 1)
		return;

	Float3 a(0.4);
	Float3 b(2.f);
	//LIMAENG::applyHyperpos(&a, &b);
	//testspace::testerfn(4);

	Float3* average_energy = new Float3[analysable_steps];
	
	for (int i = 0; i < simulation->getStep(); i++) {
		//printf("Step %d potE %f\n", i, simulation->potE_buffer[i]);
	}
	//exit(0);


																		// TODO: I think maybe i am missing 1 datapoint here? Something about only loading 99 steps in??
	cudaMalloc(&traj_buffer_device, sizeof(Float3) * simulation->total_particles_upperbound * simulation->getStep());
	cudaMalloc(&potE_buffer_device, sizeof(double) * simulation->total_particles_upperbound * simulation->getStep());
	cudaMemcpy(traj_buffer_device, simulation->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * simulation->getStep(), cudaMemcpyHostToDevice);
	cudaMemcpy(potE_buffer_device, simulation->potE_buffer, sizeof(double) * simulation->total_particles_upperbound * simulation->getStep(), cudaMemcpyHostToDevice);


	Float3* average_solvent_energy = analyzeSolvateEnergy(simulation, analysable_steps);
	Float3* average_compound_energy = analyzeCompoundEnergy(simulation, analysable_steps);	//avg energy PER PARTICLE in compound

	cudaFree(traj_buffer_device);
	cudaFree(potE_buffer_device);


	
	for (int i = 0; i < analysable_steps; i++) {
		//printf("\n");
		//average_compound_energy[i].print('c');
		//average_solvent_energy[i].print('s');
		average_energy[i] = (average_solvent_energy[i] * simulation->box->n_solvents * 1 + average_compound_energy[i] * simulation->box->compounds[0].n_particles) * (1.f/ (simulation->box->n_solvents + simulation->box->compounds[0].n_particles));
		//average_energy[i].print('a');
	}
	

	printEnergies(average_energy, analysable_steps, simulation);


	delete [] average_solvent_energy, average_compound_energy, average_energy;


	printf("\n########## Finished analyzing energies ##########\n\n");

}

Float3* Analyzer::analyzeSolvateEnergy(Simulation* simulation, int n_steps) {
	dim3 block_dim(n_steps, 1, 1);
	Float3* average_solvent_energy = new Float3[n_steps];
	//Float3* host_data = new Float3[1 * n_steps];
	Float3* data_out;
	cudaMalloc(&data_out, sizeof(Float3) * 1 * n_steps);


	//monitorSolventEnergyKernel <<< block_dim, simulation->n_solvents >>> (simulation->box, traj_buffer_device, potE_buffer_device, data_out);
	monitorSolventEnergyKernel << < block_dim, 256 >> > (simulation->box, traj_buffer_device, potE_buffer_device, data_out);
	cudaDeviceSynchronize();
	//cudaMemcpy(host_data, device_data, sizeof(Float3) * 1 * (n_steps), cudaMemcpyDeviceToHost);
	cudaMemcpy(average_solvent_energy, data_out, sizeof(Float3) * 1 * (n_steps), cudaMemcpyDeviceToHost);

	for (int step = 0; step < n_steps; step++) {
		/*for (int i = 0; i < 256; i++) {
			if (host_data[i + step * 256].x > 10000)
				printf("Block: %d energy: %f\n", i, host_data[i + step * 256].x);
			average_solvent_energy[step] += host_data[i + step * 256];
		}*/
		average_solvent_energy[step] *= (1.f / simulation->n_solvents);
	}


	cudaFree(data_out);
	//delete[] host_data;

	return average_solvent_energy;
}


Float3* Analyzer::analyzeCompoundEnergy(Simulation* simulation, int n_steps) {
	int n_datapoints = simulation->n_compounds * n_steps;
	printf("n steps: %d\n", n_steps);
	dim3 block_dim(n_steps, simulation->box->n_compounds, 1);
	Float3* average_compound_energy = new Float3[n_steps];
	Float3* host_data = new Float3[n_datapoints];
	Float3* data_out;
	cudaMalloc(&data_out, sizeof(Float3) * n_datapoints);




	monitorCompoundEnergyKernel << < block_dim, MAX_COMPOUND_PARTICLES >> > (simulation->box, traj_buffer_device, potE_buffer_device, data_out);
	cudaDeviceSynchronize();
	cudaMemcpy(host_data, data_out, sizeof(Float3) * n_datapoints, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();


	printf("This energy analysis is only valid if just 1 compound exists!\n");
	for (int step = 0; step < n_steps; step++) {
		for (int i = 0; i < simulation->box->n_compounds; i++) {
			//if (host_data[i + step * 256].x > 10000)
				//printf("Block: %d energy: %f\n", i, host_data[i + step * 256].x);
			average_compound_energy[step] += host_data[i + step * simulation->box->n_compounds];
		}
		average_compound_energy[step] *= (1.f / (simulation->box->compounds[0].n_particles));
	}


	cudaFree(data_out);
	delete[] host_data;

	return average_compound_energy;
}

void Analyzer::printEnergies(Float3* energy_data, int analysable_steps, Simulation* simulation) {
	string file_path_s = "D:\\Quantom\\energies_steps_" + to_string(analysable_steps) + ".bin";
	char* file_path;
	file_path = &file_path_s[0];
	cout << "Printing to file " << file_path << endl;

	FILE* file;
	fopen_s(&file, file_path, "wb");
	fwrite(energy_data, sizeof(Float3), analysable_steps, file);
	fclose(file);



	file_path_s = "D:\\Quantom\\temperatures_steps_" + to_string(analysable_steps) + ".bin";
	file_path = &file_path_s[0];
	cout << "Printing to file " << file_path << endl;

	fopen_s(&file, file_path, "wb");
	fwrite(simulation->temperature_buffer, sizeof(float), simulation->getStep()/STEPS_PER_THERMOSTAT, file);
	fclose(file);
}





/*

	std::ofstream myfile("D:\\Quantom\\energies_steps_" + to_string(analysable_steps) + ".csv");



	for (int i = 0; i < analysable_steps; i++) {
		for (int j = 0; j < 3; j++) {
			myfile << energy_data[i].at(j) << ';';
			//myfile << energy_data[j + i * 3] << ";";
		}
		myfile << "\n";
	}
	myfile.close();
*/


















