#include "Analyzer.cuh"





void __global__ monitorCompoundEnergyKernel(Box* box, float* data_out, int step) {
	
	__shared__ Compound compound;
	__shared__ Float3 energy[3];

	//printf("%d\n", blockIdx.x);

	if (threadIdx.x == 0) {
		compound = box->compounds[blockIdx.x];
	}
	__syncthreads();

	float potE = box->potE_buffer[0 + threadIdx.x + blockIdx.x * 3 + step * box->n_compounds * 3];	

	Float3 postsub1 = box->trajectory[threadIdx.x + blockIdx.x * blockDim.x + (step - 1) * box->n_compounds * blockDim.x];
	Float3 postplus1 = box->trajectory[threadIdx.x + blockIdx.x * blockDim.x + (step + 1) * box->n_compounds * blockDim.x];

	for (int i = 0; i < 3; i++) {
		//*postsub1.placeAt(i) += BOX_LEN * ((postplus1.x - postsub1.x) > BOX_LEN_HALF);
		//*postsub1.placeAt(i) -= BOX_LEN * ((postplus1.x - postsub1.x) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
	float vel = (postplus1 - postsub1).len() * 0.5f / box->dt;
	float kinE = 0.5 * compound.particles[threadIdx.x].mass * vel * vel;

	float totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);


	__syncthreads();
	if (threadIdx.x == 0) {
		Float3 sum = Float3(0, 0, 0);
		for (int i = 0; i < compound.n_particles; i++)
			sum = sum + energy[i] * (1.f / compound.n_particles);

		for (int i = 0; i < 3; i++) {
			data_out[i + blockIdx.x * 3] = sum.at(i);
		}		
	}
}


void __global__ monitorSolventEnergyKernel(Box* box, Float3* data_out) {

	__shared__ Float3 energy[THREADS_PER_MONITORBLOCK];
	int solvent_index = threadIdx.x + blockIdx.y * THREADS_PER_MONITORBLOCK;
	int step = blockIdx.x + 1;
	int compounds_offset = box->n_compounds * PARTICLES_PER_COMPOUND;



	if (solvent_index >= box->n_solvents) {
		return;
	}
	if (threadIdx.x == 0) {
		data_out[blockIdx.y + (step - 1) * N_MONITORBLOCKS_PER_STEP] = Float3(0, 0, 0);
	}

	Float3 pos_tsub1 = box->trajectory[compounds_offset + solvent_index + (step - 1) * box->total_particles];
	Float3 pos_tadd1 = box->trajectory[compounds_offset + solvent_index + (step + 1) * box->total_particles];


	for (int i = 0; i < 3; i++) {
		*pos_tsub1.placeAt(i) += BOX_LEN * ((pos_tadd1.at(i) - pos_tsub1.at(i)) > BOX_LEN_HALF);
		*pos_tsub1.placeAt(i) -= BOX_LEN * ((pos_tadd1.at(i) - pos_tsub1.at(i)) < -BOX_LEN_HALF);
	}


	float potE = box->potE_buffer[compounds_offset + solvent_index + step * box->total_particles];

	
	float vel = (pos_tadd1 - pos_tsub1).len() * 0.5f / box->dt;
	if (solvent_index == 1) {
		if (vel > 2000) {
			printf("step %04d solvate %04d %f\n", step, solvent_index, vel);
			//pos_tadd1.print('a');
			printf("analyzer index %d\n", compounds_offset + solvent_index + (step - 1) * box->total_particles);
			//pos_tsub1.print('s');
		}
		
	}
		
	float kinE = 0.5 * SOLVENT_MASS * vel * vel;


	float totalE = potE + kinE;

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



void Analyzer::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol
	// calculate energies separately for compounds and solvents. weigh averages based on amount of each


	int analysable_steps = simulation->n_steps - 2;
	// Solvent energies first //
	dim3 block_dim(analysable_steps, 256, 1);
	Float3* average_solvent_energy = new Float3[analysable_steps];
	Float3* host_data = new Float3[256 * analysable_steps];
	Float3* device_data;
	cudaMalloc(&device_data, sizeof(Float3) * 256 * analysable_steps);


	monitorSolventEnergyKernel << < block_dim , 256 >> > (simulation->box, device_data);
	cudaDeviceSynchronize();
	cudaMemcpy(host_data, device_data, sizeof(Float3) * 256 * (analysable_steps), cudaMemcpyDeviceToHost);

	for (int step = 0; step < analysable_steps; step++) {
		for (int i = 0; i < 256; i++) {
			average_solvent_energy[step] += host_data[i + step * 256];
		}
		average_solvent_energy[step] *= (1.f / simulation->box->n_solvents);
	}
	printEnergies(average_solvent_energy, analysable_steps);

	/*
	for (int i = 1; i < (simulation->n_steps - 1); i++) {
		if (!((i + 2) % 100))
			printf("Analyzing step %d\r", i + 2);
		monitorCompoundEnergyKernel << <n_compounds, 3 >> > (simulation->box, &device_data[data_per_step * i], i);
	}

	cudaMemcpy(host_data, device_data, sizeof(float) * n_datapoints, cudaMemcpyDeviceToHost);

	float* out_data = new float[3 * (simulation->n_steps - 2)];
	for (int step = 1; step < (simulation->n_steps - 1); step++) {
		for (int energy_type = 0; energy_type < 3; energy_type++) {
			out_data[energy_type + (step - 1) * 3] = 0;
			for (int m = 0; m < n_compounds; m++) {
				out_data[energy_type + (step - 1) * 3] += host_data[energy_type + m * 3 + step * data_per_step] / simulation->box->n_compounds;
			}
		}
	}
	




	//int n_compounds = simulation->box->n_compounds;

	printEnergies(out_data, simulation->n_steps);
*/
	cudaFree(device_data);
	delete [] host_data, average_solvent_energy;


	printf("\n########## Finished analyzing energies ##########\n\n");

}


/*
void Analyzer::analyzeEnergy1(Simulation* simulation) {	// Calculates the avg J/mol

	dim3 block_dim(simulation->n_steps, 256, 0);


	int data_per_step = 3 * simulation->box->n_compounds;
	int n_datapoints = data_per_step * simulation->n_steps;

	printf("N Compounds: %d\n", simulation->box->n_compounds);

	float* host_data = new float[n_datapoints];
	float* device_data;
	cudaMalloc(&device_data, sizeof(float) * n_datapoints);
	cudaDeviceSynchronize();
	int n_compounds = simulation->box->n_compounds;
	for (int i = 1; i < (simulation->n_steps-1); i++) {
		if (!((i+2) % 100))
			printf("Analyzing step %d\r", i+2);
		monitorCompoundEnergyKernel << <n_compounds, 3 >> > (simulation->box, &device_data[data_per_step * i], i);
	}

	cudaMemcpy(host_data, device_data, sizeof(float) * n_datapoints, cudaMemcpyDeviceToHost);

	float* out_data = new float[3 * (simulation->n_steps-2)];
	for (int step = 1; step < (simulation->n_steps-1); step++) {	
		for (int energy_type = 0; energy_type < 3; energy_type++) {
			out_data[energy_type + (step-1) * 3] = 0;
			for (int m = 0; m < n_compounds; m++) {
				out_data[energy_type + (step-1) * 3] += host_data[energy_type + m * 3 + step * data_per_step] / simulation->box->n_compounds;
			}
		}
	}

	printEnergies(out_data, simulation->n_steps);

	cudaFree(device_data);
	delete host_data, out_data;
	

	printf("\n########## Finished analyzing energies ##########\n\n");

}
*/

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
























