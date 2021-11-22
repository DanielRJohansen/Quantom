#include "Analyzer.cuh"


//using namespace what;




















void __global__ monitorEnergyKernel(Box* box, float* data_out, int step) {
	
	__shared__ Compound compound;
	__shared__ CompoundState self_state;
	__shared__ CompoundState neighborstate;	// Buffer for all neighbors state, one at a time.
	__shared__ CompoundNeighborList neighborlist;
	__shared__ Float3 utility_buffer[64];
	__shared__ Float3 hyperpos_offset;

	__shared__ Float3 energy[3];

	//printf("%d\n", blockIdx.x);

	if (threadIdx.x == 0) {
		compound = box->compounds[blockIdx.x];
		self_state = box->compound_state_array[blockIdx.x];
		neighborlist = box->compound_neighborlist_array[blockIdx.x];
	}
	__syncthreads();


	//float potE = force.len();


	float potE = box->data_buffer[0 + threadIdx.x + blockIdx.x * 3 + step * box->n_compounds * 3];	

	Float3 postsub1 = box->trajectory[threadIdx.x + blockIdx.x * blockDim.x + (step - 1) * box->n_compounds * blockDim.x];
	Float3 postplus1 = box->trajectory[threadIdx.x + blockIdx.x * blockDim.x + (step + 1) * box->n_compounds * blockDim.x];

	for (int i = 0; i < 3; i++) {
		*postsub1.placeAt(i) += BOX_LEN * ((postplus1.x - postsub1.x) > BOX_LEN_HALF);
		*postsub1.placeAt(i) -= BOX_LEN * ((postplus1.x - postsub1.x) < -BOX_LEN_HALF);
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








void Analyzer::analyzeEnergy(Simulation* simulation) {	// Calculates the avg J/mol



	int data_per_step = 3 * simulation->box->n_compounds;
	int n_datapoints = data_per_step * simulation->n_steps;

	printf("N Compounds: %d\n", simulation->box->n_compounds);

	float* host_data = new float[n_datapoints];
	float* device_data;
	cudaMalloc(&device_data, sizeof(float) * n_datapoints);
	cudaDeviceSynchronize();
	int n_compounds = simulation->box->n_compounds;
	for (int i = 1; i < (simulation->n_steps-1); i++) {
		if (!(i % 100))
			printf("Analyzing step %d\r", i);
		monitorEnergyKernel << <n_compounds, 3 >> > (simulation->box, &device_data[data_per_step * i], i);
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


void Analyzer::printEnergies(float* energy_data, int n_steps) {
	std::ofstream myfile("D:\\Quantom\\energies.csv");


	for (int i = 0; i < n_steps-2; i++) {
		for (int j = 0; j < 3; j++) {
			myfile << energy_data[j + i * 3] << ";";
		}
		myfile << "\n";
	}

	myfile.close();
}





























// ################################################################# //
// COPY OF ENGINE FUNCTIONS, DONT KNOW ANY OTHER WAY AROUND THIS.... //
// ################################################################# //



__device__ float _getAngle(Float3 v1, Float3 v2) {
	return acos((v1.dot(v2)) / (v1.len() * v2.len()));
}


__device__ void _determineMoleculerHyperposOffset(Float3* utility_buffer, Float3* compund_center, Float3* neighbor_center) {		// Called with 3 threads!! I hope i dont have to search for any bugs here! :D
	*utility_buffer[0].placeAt(threadIdx.x) =
		(
			1 * (
				(compund_center->at(threadIdx.x) - neighbor_center->at(threadIdx.x)
					) > BOX_LEN_HALF)
			- 1 * (
				neighbor_center->at(threadIdx.x) - (compund_center->at(threadIdx.x)
					) > BOX_LEN_HALF)
			)
		* BOX_LEN;
	;
}





constexpr float sigma = 0.3923;	//nm, basicllay break 0 point
constexpr float epsilon = 0.5986 * 1'000; // J/mol
__device__ Float3 _calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* data_ptr2) {	// Applying force to p0 only! Returns force in J/mol
	float dist = (*pos0 - *pos1).len();
	if (dist < *data_ptr2)
		*data_ptr2 = dist;

	if (dist < 0.000001) {
		printf("WTF");
	}

	float fraction = sigma / dist;		//nm/nm, so unitless

	float f2 = fraction * fraction;
	float f6 = f2 * f2 * f2;
	float f12 = f6 * f6;

	float LJ_pot = 4 * epsilon * (f12 - f6);
	*data_ptr += LJ_pot;

	Float3 force_unit_vector = (*pos0 - *pos1).norm();

	if (LJ_pot > 50'000) {
		printf("\n\n KILOFORCE! Block %d thread %d\n", blockIdx.x, threadIdx.x);
	}
	return force_unit_vector * LJ_pot;	//J/mol*M	(M is direction)
}

__device__ Float3 _calcLJForce(Float3* pos0, Float3* pos1) {	// Applying force to p0 only! Returns force in J/mol
	float dist = (*pos0 - *pos1).len();
	float fraction = sigma / dist;		//nm/nm, so unitless

	float f2 = fraction * fraction;
	float f6 = f2 * f2 * f2;
	float f12 = f6 * f6;

	float LJ_pot = 4 * epsilon * (f12 - f6);
	Float3 force_unit_vector = (*pos0 - *pos1).norm();

	if (LJ_pot > 50'000) {
		printf("\n\n KILOFORCE! Block %d thread %d\n", blockIdx.x, threadIdx.x);
	}
	return force_unit_vector * LJ_pot;	//J/mol*M	(M is direction)
}



constexpr float kb = 17.5 * 1e+6;		//	J/(mol*nm^2)
__device__ Float3 _calcPairbondForce(Float3* self_pos, Float3* other_pos, float* reference_dist) {
	Float3 v = *self_pos - *other_pos;
	float dif = v.len() - *reference_dist;
	float invert = v.len() > *reference_dist ? -1 : 1;
	return v.norm() * (0.5 * kb * (dif * dif) * invert);
}



constexpr float ktheta = 65 * 1e+3;	// J/mol
__device__ Float3 _calcAngleForce(CompoundState* statebuffer, AngleBond* anglebond) {	// We fix the middle particle and move the other particles so they are closest as possible
	// HOLY FUUUCK this shit is unsafe. Works if atoma are ordered left to right, with angles BELOW 180. Dont know which checks to implement yet

	Float3 v1 = statebuffer->positions[anglebond->atom_indexes[0]] - statebuffer->positions[anglebond->atom_indexes[1]];
	Float3 v2 = statebuffer->positions[anglebond->atom_indexes[2]] - statebuffer->positions[anglebond->atom_indexes[1]];
	Float3 plane_normal = v1.cross(v2);

	Float3 force_direction = anglebond->atom_indexes[0] == threadIdx.x ? plane_normal.cross(v1) : v2.cross(plane_normal);	// Calculates "inward" force orthogonal to the bond direction of the atom of interest to middle atom
	force_direction = force_direction.norm();

	float angle = _getAngle(v1, v2);
	//*dataptr = angle;					// Temp
	float dif = angle - anglebond->reference_angle;
	float force_scalar = 0.5 * ktheta * (dif * dif);

	force_scalar = dif < 0 ? force_scalar * -1 : force_scalar;	// Invert force if we need to push outwards

	return force_direction * force_scalar;
}








__device__ Float3 _computeLJForces(Box* box, Compound* compound, CompoundNeighborList* neighborlist, CompoundState* self_state, CompoundState* neighborstate_buffer, Float3* utility_buffer) {
	Float3 force(0, 0, 0);
	for (int neighbor_index = 0; neighbor_index < neighborlist->n_neighbors; neighbor_index++) {
		// ------------ Load and process neighbor molecule ------------ //
		neighborstate_buffer->n_particles = box->compound_state_array[neighborlist->neighborcompound_indexes[neighbor_index]].n_particles;						// These two lines may be optimizable.
		neighborstate_buffer->positions[threadIdx.x] = box->compound_state_array[neighborlist->neighborcompound_indexes[neighbor_index]].positions[threadIdx.x];
		__syncthreads(); // Necessary?


		// Hyperpositioning the entire molecule
		if (threadIdx.x < 3)
			_determineMoleculerHyperposOffset(&utility_buffer[0], &compound->center_of_mass, &box->compounds[neighborlist->neighborcompound_indexes[neighbor_index]].center_of_mass);
		__syncthreads();
		neighborstate_buffer->positions[threadIdx.x] = neighborstate_buffer->positions[threadIdx.x] + utility_buffer[0];
		__syncthreads();
		// ------------------------------------------------------------ //



		for (int neighbor_particle_index = 0; neighbor_particle_index < neighborstate_buffer->n_particles; neighbor_particle_index++) {
			if (threadIdx.x < compound->n_particles) {
				force = force + _calcLJForce(&self_state->positions[threadIdx.x], &neighborstate_buffer->positions[neighbor_particle_index]);
			}
		}
		__syncthreads();
	}
	return force;
}


__device__ Float3 _computePairbondForces(Compound* compound, CompoundState* self_state) {
	Float3 force(0, 0, 0);
	for (int i = 0; i < compound->n_pairbonds; i++) {
		PairBond* pb = &compound->pairbonds[i];
		if (pb->atom_indexes[0] == threadIdx.x || pb->atom_indexes[1] == threadIdx.x) {
			int other_index = pb->atom_indexes[0] != threadIdx.x ? pb->atom_indexes[0] : pb->atom_indexes[1];
			force = force + _calcPairbondForce(&self_state->positions[threadIdx.x], &self_state->positions[other_index], &pb->reference_dist);
		}
	}
	return force;
}

__device__ Float3 _computeAnglebondForces(Compound* compound, CompoundState* self_state) {
	Float3 force(0, 0, 0);
	for (int i = 0; i < compound->n_anglebonds; i++) {
		AngleBond* ab = &compound->anglebonds[i];
		if (ab->atom_indexes[0] == threadIdx.x || ab->atom_indexes[2] == threadIdx.x) {
			force = force + _calcAngleForce(self_state, ab);
		}
	}
	return force;
}