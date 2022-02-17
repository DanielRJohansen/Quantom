#include "Analyzer.cuh"






void __device__ applyHyperposA(Float3* static_particle, Float3* movable_particle) {
	for (int i = 0; i < 3; i++) {
		*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
		*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
}


void __global__ monitorCompoundEnergyKernel(Box* box, Float3* traj_buffer, double* potE_buffer, Float3* data_out) {		// everything here breaks if not all compounds are identical in particle count and particle mass!!!!!!!
	__shared__ Float3 energy[MAX_COMPOUND_PARTICLES];
	__shared__ Compound compound;


	const int step = blockIdx.x + 1;
	const int compound_index = blockIdx.y;
	energy[threadIdx.x] = Float3(0, 0, 0);

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
	applyHyperposA(&pos_tadd1, &pos_tsub1);
	

	double vel = (pos_tadd1 - pos_tsub1).len() * 0.5f / box->dt;
	double kinE = 0.5 * compound.particles[threadIdx.x].mass * vel * vel;

	double totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();



	for (int i = 1; i < MAX_COMPOUND_PARTICLES; i *= 2) {	// Distributed averaging							// Make a generic and SAFER function for this, PLEASE OK??
		Float3 temp;			// This is a lazy soluation, but maybe it is also fast? Definitely simple..
		if ((threadIdx.x + i) < MAX_COMPOUND_PARTICLES) {
			//energy[threadIdx.x] = (energy[threadIdx.x] + energy[threadIdx.x + i]);// *0.5f;	// easier to just divide by sum of solvents at host
			temp = energy[threadIdx.x] + energy[threadIdx.x + i];// *0.5f;	// easier to just divide by sum of solvents at host
		}
		__syncthreads();
		energy[threadIdx.x] = temp;
		__syncthreads();
	}


	__syncthreads();
	if (threadIdx.x == 0) {
		data_out[compound_index + (step - 1) * box->n_compounds] = energy[0];
	}
}




/*
void __global__ monitorSolventEnergyKernel(Box* box, Float3* data_out) {
	__shared__ Float3 energy[THREADS_PER_MONITORBLOCK];
	int solvent_index = threadIdx.x + blockIdx.y * THREADS_PER_MONITORBLOCK;
	int step = blockIdx.x + 1;
	int compounds_offset = box->n_compounds * PARTICLES_PER_COMPOUND;


	energy[threadIdx.x] = Float3(0, 0, 0);
	if (solvent_index >= box->n_solvents) {
		return;
	}


	if (threadIdx.x == 0) {
		data_out[blockIdx.y + (step - 1) * N_MONITORBLOCKS_PER_STEP] = Float3(0, 0, 0);
	}

	Float3 pos_tsub1 = box->trajectory[compounds_offset + solvent_index + (step - 1) * box->total_particles_upperbound];
	Float3 pos_tadd1 = box->trajectory[compounds_offset + solvent_index + (step + 1) * box->total_particles_upperbound];
	applyHyperposA(&pos_tadd1, &pos_tsub1);


	double potE = box->potE_buffer[compounds_offset + solvent_index + step * box->total_particles_upperbound];

	
	double vel = (pos_tadd1 - pos_tsub1).len() * 0.5f / box->dt;
	if (vel > 10'000) {
		printf("step %04d solvate %04dvel:  %f\n", step, solvent_index, vel);
		//pos_tadd1.print('a');
		//printf("analyzer index %d\n", compounds_offset + solvent_index + (step - 1) * box->total_particles_upperbound);
		//pos_tsub1.print('s');
	}
	if (potE > 200'000) {
		printf("step %04d solvate %04d pot %f\n", step, solvent_index, potE);
	}
	
		
	double kinE = 0.5 * SOLVENT_MASS * vel * vel;


	double totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);
	__syncthreads();

	for (int i = 1; i < THREADS_PER_MONITORBLOCK; i *= 2) {	// Distributed averaging
		if ((threadIdx.x + i) < THREADS_PER_MONITORBLOCK) {
			energy[threadIdx.x] = (energy[threadIdx.x] + energy[threadIdx.x + i]);// *0.5f;	// easier to just divide by sum of solvents at host
		}
		__syncthreads();
	}

	if (threadIdx.x == 0) {
		//energy[0].print('b');
		data_out[blockIdx.y + (step-1) * N_MONITORBLOCKS_PER_STEP] = energy[0];
		//data_out[blockIdx.y + (step - 1) * 256] = Float3(1,2,3);
		//data_out[blockIdx.y + step * 256].print('g');
	}
	
}
*/


void Analyzer::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol // calculate energies separately for compounds and solvents. weigh averages based on amount of each
	int analysable_steps = simulation->getStep() - 2;
	if (analysable_steps < 1)
		return;

	Float3* average_energy = new Float3[analysable_steps];
	

	//Float3* average_solvent_energy = analyzeSolvateEnergy(simulation, analysable_steps);
	Float3* average_compound_energy = analyzeCompoundEnergy(simulation, analysable_steps);	//avg energy PER PARTICLE in compound

	

	for (int i = 0; i < analysable_steps; i++) {
		average_compound_energy->print('e');
		//average_energy[i] = (average_solvent_energy[i] * simulation->box->n_solvents * 1 + average_compound_energy[i] * simulation->box->n_compounds * PARTICLES_PER_COMPOUND) * (1.f/ simulation->box->total_particles_upperbound);
	}


	printEnergies(average_energy, analysable_steps);


	//delete [] average_solvent_energy, average_compound_energy, average_energy;


	printf("\n########## Finished analyzing energies ##########\n\n");

}
/*
Float3* Analyzer::analyzeSolvateEnergy(Simulation* simulation, int n_steps)
{
	// Solvent energies first //
	dim3 block_dim(n_steps, 256, 1);
	Float3* average_solvent_energy = new Float3[n_steps];
	Float3* host_data = new Float3[256 * n_steps];
	Float3* device_data;
	cudaMalloc(&device_data, sizeof(Float3) * 256 * n_steps);


	monitorSolventEnergyKernel <<< block_dim, 256 >>> (simulation->box, device_data);
	cudaDeviceSynchronize();
	cudaMemcpy(host_data, device_data, sizeof(Float3) * 256 * (n_steps), cudaMemcpyDeviceToHost);

	for (int step = 0; step < n_steps; step++) {
		for (int i = 0; i < 256; i++) {
			if (host_data[i + step * 256].x > 10000)
				printf("Block: %d energy: %f\n", i, host_data[i + step * 256].x);
			average_solvent_energy[step] += host_data[i + step * 256];
		}
		average_solvent_energy[step] *= (1.f / simulation->box->n_solvents);
	}


	cudaFree(device_data);
	delete[] host_data;

	return average_solvent_energy;
}
*/
Float3* Analyzer::analyzeCompoundEnergy(Simulation* simulation, int n_steps) {
	int n_datapoints = simulation->n_compounds * n_steps;
	printf("n steps: %d\n", n_steps);
	dim3 block_dim(n_steps, simulation->box->n_compounds, 1);
	Float3* average_compound_energy = new Float3[n_steps];
	Float3* host_data = new Float3[n_datapoints];
	Float3* device_data;
	cudaMalloc(&device_data, sizeof(Float3) * n_datapoints);


	Float3* traj_buffer_device;
	double* potE_buffer_device;
	cudaMalloc(&traj_buffer_device, sizeof(Float3) * simulation->total_particles_upperbound * simulation->getStep());
	cudaMalloc(&potE_buffer_device, sizeof(double) * simulation->total_particles_upperbound * simulation->getStep());
	cudaMemcpy(traj_buffer_device, simulation->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * simulation->getStep(), cudaMemcpyHostToDevice);
	cudaMemcpy(potE_buffer_device, simulation->potE_buffer, sizeof(double) * simulation->total_particles_upperbound * simulation->getStep(), cudaMemcpyHostToDevice);

	monitorCompoundEnergyKernel << < block_dim, MAX_COMPOUND_PARTICLES >> > (simulation->box, traj_buffer_device, potE_buffer_device, device_data);
	cudaDeviceSynchronize();
	cudaMemcpy(host_data, device_data, sizeof(Float3) * n_datapoints, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();


	printf("What %d\n", simulation->box->n_compounds);
	for (int step = 0; step < n_steps; step++) {
		for (int i = 0; i < simulation->box->n_compounds; i++) {
			//if (host_data[i + step * 256].x > 10000)
				//printf("Block: %d energy: %f\n", i, host_data[i + step * 256].x);
			average_compound_energy[step] += host_data[i + step * simulation->box->n_compounds];
		}
		//average_compound_energy[step] *= (1.f / (simulation->box->n_compounds * PARTICLES_PER_COMPOUND));
	}


	cudaFree(device_data);
	delete[] host_data;

	return average_compound_energy;
}


void Analyzer::printEnergies(Float3* energy_data, int analysable_steps) {
	std::ofstream myfile("D:\\Quantom\\energies.csv");


	for (int i = 0; i < analysable_steps; i++) {
		for (int j = 0; j < 3; j++) {
			myfile << energy_data[i].at(j) << ';';
			//myfile << energy_data[j + i * 3] << ";";
		}
		myfile << "\n";
	}
	myfile.close();
}
























