#include "Engine.cuh"



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
	if (simulation->box->step==simulation->n_steps)
		return;
	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Error before step!");
		exit(1);
	}


	auto t0 = std::chrono::high_resolution_clock::now();
	forceKernel <<< simulation->box->n_compounds, MAX_COMPOUND_PARTICLES >>> (simulation->box);
	solventForceKernel <<<BLOCKS_PER_SOLVENTKERNEL, THREADS_PER_SOLVENTBLOCK >>> (simulation->box);

	cudaDeviceSynchronize();
	auto t1 = std::chrono::high_resolution_clock::now();


	//IMPORTANT TODOOOOOOOOOOOOOOO!!
	// INSTEAD OF MOVING THE DATA AROUND, FREE THE OLD POINTERS, AND SIMPLY POINT TO THE NEW LOCATION IN VRAM. THEN MAKE ALLOC A NEW ARRAY FOR NEXT STEP'S ARRAY_NEXT !
	// FOR THIS IT WOULD BE OPTIMAL TO ALSO KEEP THE PREVIOUS POSITIONS SEPARATE.

	cudaMemcpy(simulation->box->compound_state_array, simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS, cudaMemcpyDeviceToDevice);	// Update all positions, after all forces have been calculated
	cudaMemcpy(simulation->box->solvents, simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyDeviceToDevice);
	cudaDeviceSynchronize();
	auto t2 = std::chrono::high_resolution_clock::now();



	simulation->box->step++;

	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
	

	//simulation->step++;
}





// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//

__device__ double cudaMax(double a, double b) {
	if (a > b)
		return a;
	return b;
}
__device__ double cudaMin(double a, double b) {
	if (a < b)
		return a;
	return b;
}




void __device__ applyHyperpos(Float3* static_particle, Float3* movable_particle) {
	for (int i = 0; i < 3; i++) {
		*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
		*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
}

__device__ void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
	*current_position += Float3(BOX_LEN, BOX_LEN, BOX_LEN);
	*current_position = current_position->elementwiseModulus(BOX_LEN);
}




constexpr double sigma = 0.3923f;	//nm, basicllay break 0 point
constexpr double epsilon = 0.5986 * 1'000.f; // J/mol
__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, double* data_ptr) {	// Applying force to p0 only! Returns force in J/mol
	double dist = (*pos0 - *pos1).len();
	double fraction = sigma / dist;		//nm/nm, so unitless

	double f2 = fraction * fraction;
	double f6 = f2 * f2 * f2;
	double f12 = f6 * f6;

	//double force = 1.f / (dist) * f6 * (1.f - 2.f * f6);
	//double force = 24.f * (epsilon / dist) * f6 * (1.f - 2.f * f6);
	double force = 24.f * (epsilon / (dist*dist)) * f6 * (1.f - 2.f * f6);

	double LJ_pot = 4.f * epsilon * (f12 - f6);
	Float3 force_unit_vector = (*pos1 - *pos0).norm();


	//printf("Mol %d force %f\n", blockIdx.x, (force_unit_vector * force).x);
	data_ptr[0] += LJ_pot *0.5 * 2;																	// WATCH OUT FOR 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	data_ptr[1] += force;
	data_ptr[2] = cudaMin(data_ptr[2], dist);

	if (dist < 0.18f) {
	//if (abs(LJ_pot) > 7000) {
		printf("\nThread %d step %d dist %f pot %f\n", threadIdx.x, (int) data_ptr[3], (*pos0 - *pos1).len(), LJ_pot);
		//(*pos0 - *pos1).print('v');
		//(force_unit_vector * force).print('f');
		pos0->print('0');
		pos1->print('1');
		//printf("\n\n KILOFORCE! Block %d thread %d\n", blockIdx.x, threadIdx.x);
		//(force_unit_vector * force).print('f');
	}
	return force_unit_vector * force;	//J/mol*M	(M is direction)
}



constexpr double kb = 17.5 * 1e+6;		//	J/(mol*nm^2)
__device__ Float3 calcPairbondForce(Float3* self_pos, Float3* other_pos, double* reference_dist, double* potE) {
	Float3 v = *self_pos - *other_pos;
	double dif = v.len() - *reference_dist;
	double invert = v.len() > *reference_dist ? -1 : 1;
	
	Float3 force = v.norm() * kb * (-dif);

	if (abs(*potE) > 300000) {
		printf("\n");
		printf("threadid %d from: %f to: %f\n", threadIdx.x, *potE, 0.5 * kb * (dif * dif) * 0.5);

		self_pos->print('s');
		other_pos->print('o');
		printf("len %f, force %f\n", v.len(), force.len());
	}

	// Logging stuff
	if (threadIdx.x == LOGTHREAD && blockIdx.x == LOGBLOCK) {
		//self_pos[63 - threadIdx.x].x = v.len();
		//self_pos[63 - threadIdx.x].y = *potE;
	}
	*potE += 0.5 * kb * (dif * dif) * 0.5;// *1 * 10e+4;
	return force;
	//return v.norm() * (0.5 * kb * (dif * dif) * invert) *1.99 * 10e+4;
}



constexpr double ktheta = 65 * 1e+3;	// J/mol
__device__ Float3 calcAngleForce(CompoundState* statebuffer, AngleBond* anglebond, double* potE) {	// We fix the middle particle and move the other particles so they are closest as possible
	// HOLY FUUUCK this shit is unsafe. Works if atoma are ordered left to right, with angles BELOW 180. Dont know which checks to implement yet

	Float3 v1 = statebuffer->positions[anglebond->atom_indexes[0]] - statebuffer->positions[anglebond->atom_indexes[1]];
	Float3 v2 = statebuffer->positions[anglebond->atom_indexes[2]] - statebuffer->positions[anglebond->atom_indexes[1]];
	Float3 plane_normal = v1.cross(v2);

	Float3 force_direction = anglebond->atom_indexes[0] == threadIdx.x ? plane_normal.cross(v1) : v2.cross(plane_normal);	// Calculates "inward" force orthogonal to the bond direction of the atom of interest to middle atom
	force_direction = force_direction.norm();

	
		


	double angle = Float3::getAngle(v1, v2);
	potE[1] = angle;					// Temp
	double dif = angle - anglebond->reference_angle;
	dif = dif / 2.f / PI * 360.f;
	*potE += 0.5 * ktheta * dif * dif * 0.5;
	double force_scalar = ktheta * (dif);


	if (threadIdx.x == 1 && 0) {
		printf("\n");
		v1.print('1');
		v2.print('2');
		plane_normal.print('n');
		force_direction.print('d');
		printf("dif: %f\n", dif);
		printf("\nAngle %f ref %f\n", angle, anglebond->reference_angle);
		(force_direction * force_scalar).print('f');
	}
	

	
	//force_scalar = dif < 0 ? force_scalar * -1 : force_scalar;	// Invert force if we need to push outwards
	
	return force_direction * force_scalar;
}




// compoundkernel for other compounds
/*
__device__ Float3 computeLJForces(Box * box, Compound* compound, CompoundNeighborList* neighborlist, CompoundState* self_state, CompoundState* neighborstate_buffer, Float3* utility_buffer, double* data_ptr) {
	Float3 force(0, 0, 0);
	for (int neighbor_index = 0; neighbor_index < neighborlist->n_neighbors; neighbor_index++) {
		// ------------ Load and process neighbor molecule ------------ //
		neighborstate_buffer->n_particles = box->compound_state_array[neighborlist->neighborcompound_indexes[neighbor_index]].n_particles;						// These two lines may be optimizable.
		neighborstate_buffer->positions[threadIdx.x] = box->compound_state_array[neighborlist->neighborcompound_indexes[neighbor_index]].positions[threadIdx.x];
		__syncthreads(); // Necessary?


		// Hyperpositioning the entire molecule
		determineMoleculerHyperposOffset(&utility_buffer[0], &compound->center_of_mass, &box->compounds[neighborlist->neighborcompound_indexes[neighbor_index]].center_of_mass);
		neighborstate_buffer->positions[threadIdx.x] = neighborstate_buffer->positions[threadIdx.x] + utility_buffer[0];
		__syncthreads();		
		// ------------------------------------------------------------ //



		for (int neighbor_particle_index = 0; neighbor_particle_index < neighborstate_buffer->n_particles; neighbor_particle_index++) {
			if (threadIdx.x < compound->n_particles) {
				force = force + calcLJForce(&self_state->positions[threadIdx.x], &neighborstate_buffer->positions[neighbor_particle_index], data_ptr);

				if (threadIdx.x == LOGTHREAD && blockIdx.x == LOGBLOCK) {
					double len = (self_state->positions[threadIdx.x] - neighborstate_buffer->positions[neighbor_particle_index]).len();
					self_state->positions[60].x = cudaMin(self_state->positions[60].x, len);
				}
			}
		}
		__syncthreads();
	}
	//force = force * 24.f * epsilon;
	return force;
}
*/

// compoundkernel for solvents
__device__ Float3 computeLJForces(Box* box, Compound* compound, SolventNeighborList* neighborlist, CompoundState* self_state, Float3* utility_buffer, double* data_ptr) {
	Float3 force(0, 0, 0);

	// ------------ Load and process neighbor molecule ------------ //
	for (int i = 0; i < 256; i += MAX_COMPOUND_PARTICLES) {
		int solvent_index = threadIdx.x + i;
		if ((solvent_index) < neighborlist->n_neighbors) {				// Dont think we have enough threads for this...
			utility_buffer[solvent_index] = box->solvents[neighborlist->neighborsolvent_indexes[solvent_index]].pos;
			applyHyperpos(&self_state->positions[threadIdx.x], &utility_buffer[solvent_index]);
			//utility_buffer[solvent_index] = utility_buffer[solvent_index] + getSolventHyperposOffset(&self_state->positions[threadIdx.x], &utility_buffer[solvent_index]);
		}
	}
	__syncthreads(); 
	// ---------------------------------------------- ------------ //

	if (threadIdx.x < compound->n_particles) {

		for (int i = 0; i < neighborlist->n_neighbors; i++) {
			force = force + calcLJForce(&self_state->positions[threadIdx.x], &utility_buffer[i], data_ptr);			
			//force.print('f');
		}
	}

	//force = force * 24.f * epsilon;
	return force;
}

// solventkernel for solvents and compounds
__device__ Float3 computeLJForces(Box* box, Float3* self_pos, CompoundNeighborList* compoundneighbor_list, SolventNeighborList* solventneighbor_list,  double* data_ptr, Float3* utility_buffer) {
	Float3 force(0, 0, 0);

	// ------------ Load and process neighbor molecule ------------ //
	for (int i = 0; i < compoundneighbor_list->n_neighbors; i++) {
		int compound_index = compoundneighbor_list->neighborcompound_indexes[i];
		CompoundState* compoundstate = &box->compound_state_array[compound_index];

		if (threadIdx.x < compoundstate->n_particles) {
			utility_buffer[threadIdx.x] = compoundstate->positions[threadIdx.x];
			applyHyperpos(self_pos, &utility_buffer[threadIdx.x]);
		}
		__syncthreads();

		for (int j = 0; j < compoundstate->n_particles; j++) {
			force += calcLJForce(self_pos, &utility_buffer[j], data_ptr);

			/*
			float dist = (utility_buffer[j] - *self_pos).len();
			if (dist < 0.7) {
				self_pos->print('s');
				utility_buffer[j].print('o');
				printf("dist: %f force: %f\n", dist, force.len());
			}
				*/
		}
		__syncthreads();
	}

	for (int i = 0; i < solventneighbor_list->n_neighbors; i++) {
		Float3 other_pos = box->solvents[solventneighbor_list->neighborsolvent_indexes[i]].pos;
		applyHyperpos(self_pos, &other_pos);
		force += calcLJForce(self_pos, &other_pos, data_ptr);
	}
	return force;
}


__device__ Float3 computePairbondForces(Compound* compound, CompoundState* self_state, double* potE) {
	Float3 force(0, 0, 0);
	for (int i = 0; i < compound->n_pairbonds; i++) {
		PairBond* pb = &compound->pairbonds[i];
		if (pb->atom_indexes[0] == threadIdx.x || pb->atom_indexes[1] == threadIdx.x) {
			int other_index = pb->atom_indexes[0] != threadIdx.x ? pb->atom_indexes[0] : pb->atom_indexes[1];
			force = force + calcPairbondForce(&self_state->positions[threadIdx.x], &self_state->positions[other_index], &pb->reference_dist, potE);
		}
	}
	return force;
}

__device__ Float3 computeAnglebondForces(Compound* compound, CompoundState* self_state, double* potE) {
	Float3 force(0, 0, 0);
	for (int i = 0; i < compound->n_anglebonds; i++) {
		AngleBond* ab = &compound->anglebonds[i];
		if (ab->atom_indexes[0] == threadIdx.x || ab->atom_indexes[2] == threadIdx.x) {
			force = force + calcAngleForce(self_state, ab, potE);
		}
	}
	return force;
}


__device__ void integratePosition(Float3* pos, Float3* pos_tsub1, Float3* force, double mass, double dt) {
	Float3 temp = *pos;
	*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;	// no *0.5?
	*pos_tsub1 = temp;
	//printf("force: %f\n", force->len());
}




// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//





__global__ void forceKernel(Box* box) {
	__shared__ Compound compound;
	__shared__ CompoundState self_state;
	__shared__ CompoundState neighborstate_buffer;
	__shared__ CompoundNeighborList neighborlist;
	//__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];
	__shared__ Float3 utility_buffer[256];							// waaaaay to big, but needed for PCB of the many solvent neighbors!
	__shared__ Float3 hyperpos_offset;
	
	if (threadIdx.x == 0) {
		compound = box->compounds[blockIdx.x];
		self_state.n_particles = box->compound_state_array[blockIdx.x].n_particles;
		neighborlist = box->compound_neighborlist_array[blockIdx.x];
	}
	__syncthreads();
	if (threadIdx.x < compound.n_particles)
		self_state.positions[threadIdx.x] = box->compound_state_array[blockIdx.x].positions[threadIdx.x];
	else
		self_state.positions[threadIdx.x] = Float3(0, 0, 0);

	//printf("n neighbors")



	double data_ptr[4];
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999;




	Float3 force(0, 0, 0);
	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //
	//force = force + computeLJForces(box, &compound, &neighborlist, &self_state, &neighborstate_buffer, utility_buffer, data_ptr);
	Float3 waterforce = computeLJForces(box, &compound, &box->solvent_neighborlist_array[blockIdx.x], &self_state, utility_buffer, data_ptr);
	force += waterforce;
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //
	if (blockIdx.x == LOGBLOCK && threadIdx.x == LOGTHREAD && LOGTYPE == 1) {
		//waterforce.print('\n');
	}


	// ------------------------------------------------------------ Intramolecular forces ------------------------------------------------------------ //
	//force = force + computePairbondForces(&compound, &self_state, &data_ptr[2]);

	//force = force + computeAnglebondForces(&compound, &self_state, &data_ptr[2]);
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //





	
	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	if (threadIdx.x < compound.n_particles) {
		integratePosition(&self_state.positions[threadIdx.x], &compound.particles[threadIdx.x].pos_tsub1, &force, compound.particles[threadIdx.x].mass, box->dt);
		box->compounds[blockIdx.x].particles[threadIdx.x].pos_tsub1 = compound.particles[threadIdx.x].pos_tsub1;
	}
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------- //
	




	// ------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------------------------- // 
	{
		/*
		//for (int i = threadIdx.x; i < 256; i+= MAX_COMPOUND_PARTICLES)
			//utility_buffer[i] = Float3(0, 0, 0);
		utility_buffer[threadIdx.x] = Float3(0, 0, 0);
		if (threadIdx.x < compound.n_particles)
			utility_buffer[threadIdx.x] = self_state.positions[threadIdx.x];
		__syncthreads();

		for (int i = 1; i < 64; i *= 2) {
			if (threadIdx.x + i < 64)
				utility_buffer[threadIdx.x] = utility_buffer[threadIdx.x] + utility_buffer[threadIdx.x + i];
			__syncthreads();
		}

		if (threadIdx.x == 0) {
			if (utility_buffer[0].x < -10) {
				printf("\nBlock %d Thread %d\n", blockIdx.x, threadIdx.x);
				utility_buffer[0].print('u');
				force.print('f');
				self_state.positions[0].print('0');
				self_state.positions[1].print('1');
				self_state.positions[2].print('2');
				compound.particles[0].pos_tsub1.print('a');
				compound.particles[1].pos_tsub1.print('b');
				compound.particles[2].pos_tsub1.print('c');
			}

			utility_buffer[0] = utility_buffer[0] * (1.f / compound.n_particles);
			hyperpos_offset.x = 1 * (utility_buffer[0].x < 0) - 1 * (utility_buffer[0].x > BOX_LEN);
			hyperpos_offset.y = 1 * (utility_buffer[0].y < 0) - 1 * (utility_buffer[0].y > BOX_LEN);
			hyperpos_offset.z = 1 * (utility_buffer[0].z < 0) - 1 * (utility_buffer[0].z > BOX_LEN);
			hyperpos_offset = hyperpos_offset * Float3(BOX_LEN, BOX_LEN, BOX_LEN);
			box->compounds[blockIdx.x].center_of_mass = utility_buffer[0] + hyperpos_offset;
		}
		__syncthreads();
		applyPBC()
		//self_state.positions[threadIdx.x] = self_state.positions[threadIdx.x] + hyperpos_offset;	// Uhh, move hyperpos to utilbuffer?
		*/
		applyPBC(&self_state.positions[threadIdx.x]);
		if (threadIdx.x == 0) {
			box->compounds[blockIdx.x].center_of_mass = self_state.positions[0];	//// yeeeeeeeeeeeeeeeeeeaaah, dunno.
		}
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
	


	
	// ------------------------------------ DATA LOG ------------------------------- //
	{
		if (threadIdx.x < 3) {	// TODO: UPDATE 3 TO PARTICLES_PER_COMPOUND
			box->potE_buffer[threadIdx.x + blockIdx.x * PARTICLES_PER_COMPOUND + box->step * box->total_particles] = data_ptr[0] + data_ptr[2];
			box->trajectory[threadIdx.x + blockIdx.x * PARTICLES_PER_COMPOUND + (box->step) * box->total_particles] = self_state.positions[threadIdx.x];
		}
		__syncthreads();

		if (blockIdx.x == LOGBLOCK && threadIdx.x == LOGTHREAD && LOGTYPE == 1) {
			//box->outdata[2 + box->step * 10] = potE;	//pairbond force
			//box->outdata[1 + box->step * 10] = data_ptr[3];	// Angles
			//box->outdata[2 + box->step * 10] = data_ptr[2];
			box->outdata[3 + box->step * 10] = data_ptr[0];	// LJ pot


			Float3 pos_temp = compound.particles[threadIdx.x].pos_tsub1;
			applyHyperpos(&self_state.positions[threadIdx.x], &pos_temp);
			box->outdata[4 + box->step * 10] = (self_state.positions[threadIdx.x] - pos_temp).len() / box->dt;

			box->outdata[5 + box->step * 10] = data_ptr[2];// closest particle
			box->outdata[6 + box->step * 10] = data_ptr[1];// force.len();


			box->outdata[7 + box->step * 10] = waterforce.x;
			box->outdata[8 + box->step * 10] = waterforce.y;
			box->outdata[9 + box->step * 10] = waterforce.z;
		}
	}
	
	// ----------------------------------------------------------------------------- //


	

	if (threadIdx.x == 0) {
		box->compound_state_array_next[blockIdx.x] = self_state;
	}	
}




__global__ void solventForceKernel(Box* box) {
	__shared__ Float3 utility_buffer[THREADS_PER_SOLVENTBLOCK];


	int solvent_index = threadIdx.x + blockIdx.x * blockDim.x;
	if (solvent_index >= box->n_solvents) { return; }			// No use, we need all threads later for pbc with compounds.. find another solution

	double data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;
	Float3 force(0,0,0);


	Solvent solvent = box->solvents[solvent_index];
	data_ptr[3] = box->step;
	force += computeLJForces(box, &solvent.pos, &box->compound_neighborlist_array[MAX_COMPOUNDS + solvent_index], &box->solvent_neighborlist_array[MAX_COMPOUNDS + solvent_index], data_ptr, utility_buffer);

	integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, SOLVENT_MASS, box->dt);

	
	applyPBC(&solvent.pos);	// forcePositionToInsideBox





	// ------------------------------------ DATA LOG ------------------------------- //
	{
		int compounds_offset = box->n_compounds * PARTICLES_PER_COMPOUND;
		box->potE_buffer[compounds_offset + solvent_index + (box->step) * box->total_particles] = data_ptr[0];
		//printf("pot: %f\n", box->potE_buffer[compounds_offset + solvent_index + (box->step) * box->total_particles]);
		box->trajectory[compounds_offset + solvent_index + (box->step) * box->total_particles] = solvent.pos;
		if (solvent_index == LOGTHREAD && LOGTYPE == 0) {
			//printf("\nindex: %d\n", compounds_offset + solvent_index + (box->step) * box->total_particles);			
			//solvent.pos.print('p');
			box->outdata[3 + box->step * 10] = data_ptr[0];	// LJ pot
			box->outdata[4 + box->step * 10] = (solvent.pos-solvent.pos_tsub1).len()/box->dt;	// approx. vel
			box->outdata[5 + box->step * 10] = data_ptr[2];	// closest particle	
			//printf("%f\n", data_ptr[2]);
		}
	}

	// ----------------------------------------------------------------------------- //



	box->solvents_next[solvent_index] = solvent;
}


























/*	VELOCITY VERLET STORMER
__device__ void integratePosition(CompactParticle* particle, Float3* particle_pos, Float3* particle_force, double* dt) {
	*particle_pos = *particle_pos + (particle->vel + *particle_force * (0.5 / particle->mass) * *dt) * *dt;
}
__device__ void integrateVelocity(CompactParticle* particle, Float3* particle_force, double* dt) {
	particle->vel = particle->vel + (*particle_force + particle->force_prev) * (0.5 / particle->mass) * *dt;
}
*/