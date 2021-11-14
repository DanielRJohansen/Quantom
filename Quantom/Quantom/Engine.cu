#include "Engine.cuh"


int thing() {
	return 4;
}

Simulation* Engine::prepSimulation(Simulation* simulation) {
	this->simulation = simulation;
	srand(290128301);
	boxbuilder.build(simulation);
	printf("Boxbuild complete!\n");

	//updateNeighborLists();
	//printf("Neighborlists ready\n");



	simulation->moveToDevice();

	//initKernel << < simulation->box->n_compounds, 64 >> > (simulation->box);
	//cudaDeviceSynchronize();

	return this->simulation;
}

float* Engine::getDatabuffer() {
	float* host_data;
	int n_datapoints = simulation->box->n_compounds * 3 * 2 * 10000;
	host_data = new float[n_datapoints];

	//cudaMemcpy(host_data, box->data_buffer, sizeof(float) * box->n_compounds * 3 * 2 * 10000, cudaMemcpyDeviceToHost);
	for (int i = 0; i < n_datapoints; i++)
		host_data[i] = 0;
	printf("Copying %d floats, %d MB\n", n_datapoints, (int) (n_datapoints * sizeof(float)/1000000.f));
	cudaMemcpy(host_data, simulation->box->data_buffer, sizeof(float) * n_datapoints, cudaMemcpyDeviceToHost);
	
	cudaDeviceSynchronize();
	return host_data;
}




void Engine::updateNeighborLists() {	// Write actual function later;
	/*
	int maxc = 1'000'000; // this is temporary!
	CompoundState* statebuffer_host = new CompoundState[maxc];
	CompoundNeighborList* neighborlists_host = new CompoundNeighborList[maxc];
	cudaMemcpy(statebuffer_host, simulation->box->compound_state_array, sizeof(CompoundState) * maxc, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();


	// This only needs to be done the first time... Or does it????
	for (int i = 0; i < maxc; i++) {
		//neighborlists_host[i].n_neighbors = 0;
	}
		



	cudaMemcpy(simulation->box->compound_neighborlist_array, neighborlists_host, sizeof(CompoundNeighborList) * maxc, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	*/
}





//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//


void Engine::step() {
	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Error before step!");
		exit(1);
	}


	auto t0 = std::chrono::high_resolution_clock::now();
	forceKernel <<< simulation->box->n_compounds, 64 >>> (simulation->box, testval++);
	cudaDeviceSynchronize();
	auto t1 = std::chrono::high_resolution_clock::now();

	cudaMemcpy(simulation->box->compound_state_array, simulation->box->compound_state_array_next, sizeof(CompoundState) * boxbuilder.max_compounds, cudaMemcpyDeviceToDevice);	// Update all positions, after all forces have been calculated
	cudaDeviceSynchronize();
	auto t2 = std::chrono::high_resolution_clock::now();



	simulation->box->step++;

	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
	


	simulation->step++;
}





// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//

__device__ float cudaMax(float a, float b) {
	if (a > b)
		return a;
	return b;
}
__device__ float cudaMin(float a, float b) {
	if (a < b)
		return a;
	return b;
}


__device__ Float3 getHyperPosition(Float3 pos) {			// Dont thjínk this is right anymore??????
	pos = pos + Float3(BOX_LEN, BOX_LEN, BOX_LEN) * 1.5;
	pos = pos.elementwiseModulus(BOX_LEN);
	return pos - Float3(BOX_LEN, BOX_LEN, BOX_LEN) * 0.5;
}

__device__ Float3 getClosestMirrorPos(Float3 pos, Float3 block_center) {	// Block center can also just be whatever particle you want the pos to be close to...
	Float3 dists = block_center - pos;
	Float3 abs_dists = dists.abs();
	abs_dists = abs_dists + Float3(0.000001 * !(abs_dists.x), 0.000001 * !(abs_dists.y), 0.000001 * !(abs_dists.z));
	Float3 directional_hotencoded_vector(
		round(dists.x / abs_dists.x) * (abs_dists.x > BOX_LEN_HALF),
		round(dists.y / abs_dists.y) * (abs_dists.y > BOX_LEN_HALF),
		round(dists.z / abs_dists.z) * (abs_dists.z > BOX_LEN_HALF)
	);
	return pos + directional_hotencoded_vector * Float3(BOX_LEN, BOX_LEN, BOX_LEN);
}


__device__ float getAngle(Float3 v1, Float3 v2) {
	return acos((v1.dot(v2)) / (v1.len() * v2.len()));
}


__device__ void determineMoleculerHyperposOffset(Float3* utility_buffer, Float3* compund_center, Float3* neighbor_center) {		// Called with 3 threads!! I hope i dont have to search for any bugs here! :D
	*utility_buffer[0].placeAt(threadIdx.x) =
		(
			1 * (
					(compund_center->at(threadIdx.x) - neighbor_center->at(threadIdx.x)
				) > BOX_LEN_HALF)
			- 1 *	(
				neighbor_center->at(threadIdx.x) - (compund_center->at(threadIdx.x)
					) > BOX_LEN_HALF)
		)
		* BOX_LEN;
		;
}





constexpr float sigma = 0.3923;	//nm, basicllay break 0 point
constexpr float epsilon = 0.5986 * 1'000; // J/mol
__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* data_ptr2) {	// Applying force to p0 only! Returns force in J/mol
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

__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1) {	// Applying force to p0 only! Returns force in J/mol
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



constexpr float kb = 17.5 * 10e+6;		//	J/(mol*nm^2)
__device__ Float3 calcPairbondForce(Float3* self_pos, Float3* other_pos, float* reference_dist) {
	Float3 v = *self_pos - *other_pos;
	float dif = v.len() - *reference_dist;
	float invert = v.len() > *reference_dist ? -1 : 1;
	return v.norm() * (0.5 * kb * (dif * dif) * invert);
}



constexpr float ktheta = 65 * 10e+3;	// J/mol
__device__ Float3 calcAngleForce(CompoundState* statebuffer, AngleBond* anglebond) {	// We fix the middle particle and move the other particles so they are closest as possible
	// HOLY FUUUCK this shit is unsafe. Works if atoma are ordered left to right, with angles BELOW 180. Dont know which checks to implement yet

	Float3 v1 = statebuffer->positions[anglebond->atom_indexes[0]] - statebuffer->positions[anglebond->atom_indexes[1]];
	Float3 v2 = statebuffer->positions[anglebond->atom_indexes[2]] - statebuffer->positions[anglebond->atom_indexes[1]];
	Float3 plane_normal = v1.cross(v2);

	Float3 force_direction = anglebond->atom_indexes[0] == threadIdx.x ? plane_normal.cross(v1) : v2.cross(plane_normal);	// Calculates "inward" force orthogonal to the bond direction of the atom of interest to middle atom
	force_direction = force_direction.norm();

	float angle = getAngle(v1, v2);
	//*dataptr = angle;					// Temp
	float dif = angle - anglebond->reference_angle;
	float force_scalar = 0.5 * ktheta * (dif * dif);

	force_scalar = dif < 0 ? force_scalar * -1 : force_scalar;	// Invert force if we need to push outwards
	
	return force_direction * force_scalar;
}


__device__ void integratePosition(CompactParticle* particle, Float3* particle_pos, Float3* particle_force, float* dt) {
	Float3 temp = *particle_pos;
	*particle_pos = *particle_pos * 2 - particle->pos_tsub1 + *particle_force * (1000 / particle->mass) * *dt * *dt;	// no *0.5?
	particle->pos_tsub1 = temp;
}








__device__ Float3 computeLJForces(Box * box, Compound_H2O* compound, CompoundNeighborList* neighborlist, CompoundState* self_state, CompoundState* neighborstate_buffer, Float3* utility_buffer) {
	Float3 force(0, 0, 0);
	for (int neighbor_index = 0; neighbor_index < neighborlist->n_neighbors; neighbor_index++) {
		// ------------ Load and process neighbor molecule ------------ //
		neighborstate_buffer->n_particles = box->compound_state_array[neighborlist->neighborcompound_indexes[neighbor_index]].n_particles;						// These two lines may be optimizable.
		neighborstate_buffer->positions[threadIdx.x] = box->compound_state_array[neighborlist->neighborcompound_indexes[neighbor_index]].positions[threadIdx.x];
		__syncthreads(); // Necessary?


		// Hyperpositioning the entire molecule
		if (threadIdx.x < 3)
			determineMoleculerHyperposOffset(&utility_buffer[0], &compound->center_of_mass, &box->compounds[neighborlist->neighborcompound_indexes[neighbor_index]].center_of_mass);
		__syncthreads();
		neighborstate_buffer->positions[threadIdx.x] = neighborstate_buffer->positions[threadIdx.x] + utility_buffer[0];
		__syncthreads();
		// ------------------------------------------------------------ //



		for (int neighbor_particle_index = 0; neighbor_particle_index < neighborstate_buffer->n_particles; neighbor_particle_index++) {
			if (threadIdx.x < compound->n_particles) {
				force = force + calcLJForce(&self_state->positions[threadIdx.x], &neighborstate_buffer->positions[neighbor_particle_index]);
			}
		}
		__syncthreads();
	}
	return force;
}

__device__ Float3 computePairbondForces(Compound_H2O* compound, CompoundState* self_state) {
	Float3 force(0, 0, 0);
	for (int i = 0; i < compound->n_pairbonds; i++) {
		PairBond* pb = &compound->pairbonds[i];
		if (pb->atom_indexes[0] == threadIdx.x || pb->atom_indexes[1] == threadIdx.x) {
			int other_index = pb->atom_indexes[0] != threadIdx.x ? pb->atom_indexes[0] : pb->atom_indexes[1];
			force = force + calcPairbondForce(&self_state->positions[threadIdx.x], &self_state->positions[other_index], &pb->reference_dist);
		}
	}
	return force;
}

__device__ Float3 computeAnglebondForces(Compound_H2O* compound, CompoundState* self_state) {
	Float3 force(0, 0, 0);
	for (int i = 0; i < compound->n_anglebonds; i++) {
		AngleBond* ab = &compound->anglebonds[i];
		if (ab->atom_indexes[0] == threadIdx.x || ab->atom_indexes[2] == threadIdx.x) {
			force = force + calcAngleForce(self_state, ab);
		}
	}
	return force;
}



// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//





__global__ void forceKernel(Box* box, int step_test) {
	__shared__ Compound_H2O compound;
	__shared__ CompoundState self_state;
	__shared__ CompoundState neighborstate_buffer;
	__shared__ CompoundNeighborList neighborlist;
	__shared__ Float3 utility_buffer[64];
	__shared__ Float3 hyperpos_offset;
	
	if (threadIdx.x == 0) {
		compound = box->compounds[blockIdx.x];
		self_state = box->compound_state_array[blockIdx.x];
		neighborlist = box->compound_neighborlist_array[blockIdx.x];
	}
	
	__syncthreads();

	//bool thread_particle_active = (compound.n_particles > threadIdx.x);


	Float3 force(0, 0, 0);
	float data_ptr1 = 404;
	float data_ptr2 = 0;
	float data_ptr3 = 404;



	
	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //
	force = force + computeLJForces(box, &compound, &neighborlist, &self_state, &neighborstate_buffer, utility_buffer) * 0.1;
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //







	// ------------------------------------------------------------ Intramolecular forces ------------------------------------------------------------ //
	//float intramolecular_potential = 0;
	Float3 intramolecular_force = Float3(0,0,0);

	intramolecular_force = intramolecular_force + computePairbondForces(&compound, &self_state);

	intramolecular_force = intramolecular_force + computeAnglebondForces(&compound, &self_state);
	

	force = force + intramolecular_force;
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //





	
	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	if (threadIdx.x < compound.n_particles) {
		integratePosition(&compound.particles[threadIdx.x], &self_state.positions[threadIdx.x], &force, &box->dt);
		box->compounds[blockIdx.x].particles[threadIdx.x].pos_tsub1 = compound.particles[threadIdx.x].pos_tsub1;
	}
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------- //
	




	

	
	// ------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------------------------- //
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	if (threadIdx.x < compound.n_particles)
		utility_buffer[threadIdx.x] = self_state.positions[threadIdx.x];
	__syncthreads();
	
	for (int i = 1; i < 64; i*=2) {
		if (threadIdx.x + i < 64)	
			utility_buffer[threadIdx.x] = utility_buffer[threadIdx.x] + utility_buffer[threadIdx.x + i];
		__syncthreads();
	}

	if (threadIdx.x == 0) {
		if (utility_buffer[0].x < -10) {
			printf("\nBlock %d Thread %d\n", blockIdx.x, threadIdx.x);
			utility_buffer[0].print('u');
		}
			
		utility_buffer[0] = utility_buffer[0] * (1.f / compound.n_particles);
		hyperpos_offset.x = 1 * (utility_buffer[0].x < 0) - 1 * (utility_buffer[0].x > BOX_LEN);
		hyperpos_offset.y = 1 * (utility_buffer[0].y < 0) - 1 * (utility_buffer[0].y > BOX_LEN);
		hyperpos_offset.z = 1 * (utility_buffer[0].z < 0) - 1 * (utility_buffer[0].z > BOX_LEN);
		hyperpos_offset = hyperpos_offset * Float3(BOX_LEN, BOX_LEN, BOX_LEN);
		box->compounds[blockIdx.x].center_of_mass = utility_buffer[0] + hyperpos_offset;
	}
	__syncthreads();

	self_state.positions[threadIdx.x] = self_state.positions[threadIdx.x] + hyperpos_offset;	// Uhh, move hyperpos to utilbuffer?	
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
	





	
	// ------------------------------------ DATA LOG ------------------------------- //

	if (threadIdx.x < 3) {
		//float kin_e = 0.5 * compound.particles[threadIdx.x].mass * compound.particles[threadIdx.x].vel.len() * compound.particles[threadIdx.x].vel.len();
		box->data_buffer[0 + threadIdx.x * 2 + blockIdx.x * 3 * 2 + step_test * box->n_compounds * 3 * 2] = data_ptr2 + intramolecular_force.len();
		//box->data_buffer[1 + threadIdx.x * 2 + blockIdx.x * 3 * 2 + box->step * box->n_compounds * 3 * 2] =  kin_e;
		box->trajectory[threadIdx.x + blockIdx.x * 3 + (box->step) * box->n_compounds * 3] = self_state.positions[threadIdx.x];
	}
	__syncthreads();

	if (blockIdx.x == LOGBLOCK && threadIdx.x == LOGTHREAD && box->step < 10000) {
		
		//box->outdata[0 + box->step * 10] = bondlen;
		box->outdata[1 + box->step * 10] = data_ptr1;	// Angles
		box->outdata[2 + box->step * 10] = intramolecular_force.len();
		box->outdata[3 + box->step * 10] = data_ptr2;	// LJ pot
		box->outdata[4 + box->step * 10] = -1;//compound.particles[threadIdx.x].vel.len();
		box->outdata[5 + box->step * 10] = data_ptr3;	// closest particle
		
	}
	// ----------------------------------------------------------------------------- //


	
	__syncthreads();

	if (threadIdx.x == 0) {
		box->compound_state_array_next[blockIdx.x] = self_state;
	}
		
}






/*	VELOCITY VERLET STORMER
__device__ void integratePosition(CompactParticle* particle, Float3* particle_pos, Float3* particle_force, float* dt) {
	*particle_pos = *particle_pos + (particle->vel + *particle_force * (0.5 / particle->mass) * *dt) * *dt;
}
__device__ void integrateVelocity(CompactParticle* particle, Float3* particle_force, float* dt) {
	particle->vel = particle->vel + (*particle_force + particle->force_prev) * (0.5 / particle->mass) * *dt;
}
*/