#include "Analyzer.cuh"








void __global__ monitorEnergyKernel(Box* box, float* data_out, int step) {
	__shared__ Compound_H2O compound;
	__shared__ CompoundState self_state;
	__shared__ CompoundState neighborstate;	// Buffer for all neighbors state, one at a time.
	__shared__ CompoundNeighborList neighborlist;
	__shared__ Float3 utility_buffer[64];
	__shared__ Float3 hyperpos_offset;

	__shared__ Float3 energy[3];

	if (threadIdx.x == 0) {
		compound = box->compounds[blockIdx.x];
		self_state = box->compound_state_array[blockIdx.x];
		neighborlist = box->compound_neighborlist_array[blockIdx.x];
	}


	Float3 force(0, 0, 0);

	getAngle(Float3(1, 2, 3), Float3(2, 2, 2));
/*

	force = force + computeLJForces(box, &compound, &neighborlist, &self_state, &neighborstate, utility_buffer);
	force = force + computePairbondForces(&compound, &self_state);
	force = force + computeAnglebondForces(&compound, &self_state);
*/
	float potE = force.len();


	float vel = (box->trajectory[threadIdx.x + blockIdx.x * blockDim.x + (step + 1) * box->n_compounds * blockDim.x] - box->trajectory[threadIdx.x + blockIdx.x * blockDim.x + (step - 1) * box->n_compounds * blockDim.x]).len() * 0.5 * 1 / box->dt;
	float kinE = vel * compound.particles[threadIdx.x].mass * 1 / 1000.f;
	float totalE = potE + kinE;

	energy[threadIdx.x] = Float3(potE, kinE, totalE);

	__syncthreads();
	if (threadIdx.x == 0) {
		Float3 sum = energy[0] + energy[1] + energy[2];
		for (int i = 0; i < 3; i++) {
			data_out[i + blockIdx.x * 3] = sum.at(i);
		}		
	}

	
}


void Analyzer::analyzeEnergy(Simulation* simulation) {
	int data_per_step = 3 * simulation->box->n_compounds;
	int n_datapoints = data_per_step * simulation->n_steps;

	
	float* host_data = new float[n_datapoints];
	float* device_data;
	cudaMalloc(&device_data, sizeof(float) * n_datapoints);

	for (int i = 0; i < simulation->n_steps; i++) {
		//monitorEnergyKernel<<<simulation->box->n_compounds, 3>>>(simulation->box, &device_data[data_per_step * i], i);
	}

	cudaMemcpy(host_data, device_data, sizeof(float) * n_datapoints, cudaMemcpyDeviceToHost);

	float* out_data = new float[3 * simulation->n_steps];
	for (int step = 0; step < simulation->n_steps; step++) {
		for (int energy_type = 0; energy_type < 3; energy_type++) {
			out_data[energy_type + step * 3] = 0;
			for (int m = 0; m < simulation->box->n_compounds; m++) {
				out_data[energy_type + step * 3] += host_data[energy_type + m * 3 + step * data_per_step] / simulation->box->n_compounds;
			}
		}
	}

	printEnergies(out_data, simulation->n_steps);

	cudaFree(device_data);
	delete host_data, out_data;
	
}


void Analyzer::printEnergies(float* energy_data, int n_steps) {
	std::ofstream myfile("D:\\Quantom\\energies.csv");


	for (int i = 0; i < n_steps; i++) {
		for (int j = 0; j < 3; j++) {
			myfile << energy_data[j + i * 3] << ";";
		}
		myfile << "\n";
	}

	myfile.close();
}