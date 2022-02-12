#include "Engine.cuh"


Engine::Engine() {}
Engine::Engine(Simulation* simulation) {
	genericErrorCheck("Error before engine initialization.\n");

	this->simulation = simulation;

	nlist_data_collection = new NListDataCollection(simulation);





	int kernel_shared_mem = sizeof(Compound) + sizeof(CompoundState) + sizeof(NeighborList) + sizeof(Float3) * NEIGHBORLIST_MAX_SOLVENTS;
	printf("Forcekernel shared mem. size: %d B\n", kernel_shared_mem);


	printf("Engine ready\n");
}





void Engine::deviceMaster() {
	genericErrorCheck("Error before step!");
	step();
}

void Engine::hostMaster() {
	if (updated_neighborlists_ready) {
		cudaMemcpy(simulation->box->compound_neighborlists, nlist_data_collection->compound_neighborlists, sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyHostToDevice);
		cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);
		updated_neighborlists_ready = 0;
	}


	if (simulation->getStep() - prev_nlist_update_step >= STEPS_PER_NLIST_UPDATE) {
		offLoadPositionData(simulation);
		// Lots of waiting time here...
		cudaDeviceSynchronize();

		
		std::thread nlist_worker(Engine::updateNeighborLists, simulation, nlist_data_collection, &updated_neighborlists_ready, &timings.z);
		nlist_worker.detach();

		prev_nlist_update_step = simulation->getStep();
	}
}









//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//


void Engine::offLoadPositionData(Simulation* simulation) {
	//cudaMemcpyAsync(compoundstatearray_host, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->box->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(nlist_data_collection->compoundstates, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(nlist_data_collection->solvents, simulation->box->solvents, sizeof(Solvent) * simulation->n_solvents, cudaMemcpyDeviceToHost);
}

void __device__ __host__ applyHyperpos(Float3* static_particle, Float3* movable_particle);
void Engine::cullDistantNeighbors(NListDataCollection* nlist_data_collection) {	// Calling with nlist avoids writing function for both solvent and compound
	for (int i = 0; i < nlist_data_collection->n_compounds; i++) {
		int id_self = i;
		NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[i];
		Float3 pos_self = nlist_data_collection->compound_key_positions[i];

		for (int j = 0; j < nlist_self->n_compound_neighbors; j++) {		// Cull compound-compound
			int id_neighbor = nlist_self->neighborcompound_ids[j];

			if (id_self < id_neighbor) {
				Float3 pos_neighbor = nlist_data_collection->compound_key_positions[id_neighbor];
				NeighborList* nlist_neighbor = &nlist_data_collection->compound_neighborlists[id_neighbor];
				applyHyperpos(&pos_self, &pos_neighbor);
				double dist = (pos_neighbor - pos_self).len();
				if (dist > CUTOFF) {
					nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::COMPOUND);
					nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
					j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
				}
			}
		}


		for (int j = 0; j < nlist_self->n_solvent_neighbors; j++) {			// Cull compound-solvent
			int id_neighbor = nlist_self->neighborsolvent_ids[j];


			Float3 pos_neighbor = nlist_data_collection->solvent_positions[id_neighbor];
			NeighborList* nlist_neighbor = &nlist_data_collection->solvent_neighborlists[id_neighbor];

			applyHyperpos(&pos_self, &pos_neighbor);
			double dist = (pos_neighbor - pos_self).len();
			if (dist > CUTOFF) {
				nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
				j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
			}
		}
	}

	for (int i = 0; i < nlist_data_collection->n_solvents; i++) {																// Cull solvent-solvent
		int id_self = i;
		NeighborList* nlist_self = &nlist_data_collection->solvent_neighborlists[i];
		Float3 pos_self = nlist_data_collection->solvent_positions[i];


		for (int j = 0; j < nlist_self->n_compound_neighbors; j++) {			/// NOT FINISHED HERE
			int id_neighbor = nlist_self->neighborsolvent_ids[j];

			if (id_self < id_neighbor) {
				Float3 pos_neighbor = nlist_data_collection->solvent_positions[id_neighbor];
				NeighborList* nlist_neighbor = &nlist_data_collection->solvent_neighborlists[id_neighbor];

				applyHyperpos(&pos_self, &pos_neighbor);
				double dist = (pos_neighbor - pos_self).len();
				if (dist > CUTOFF) {
					nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
				}
			}
		}
	}
}

void Engine::updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, volatile bool* finished, int* timing) {	// This is a thread worker-function, so it can't own the object, thus i pass a ref to the engine object..
	auto t0 = std::chrono::high_resolution_clock::now();

	nlist_data_collection->compressPositionData();

	// First do culling of neighbors that has left CUTOFF
	cullDistantNeighbors(nlist_data_collection);



	// First add compound->solvent, compound->compound
	for (int id_self = 0; id_self < simulation->n_compounds; id_self++) {										
		NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
		HashTable hashtable_compoundneighbors(nlist_self->neighborcompound_ids, (int)nlist_self->n_compound_neighbors, NEIGHBORLIST_MAX_COMPOUNDS * 2);
		HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, (int)nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
		Float3 pos_self = nlist_data_collection->compound_key_positions[id_self];



		for (int i = 0; i < nlist_self->n_compound_neighbors; i++) {											// Use neighbor compounds to find new solvents
			NeighborList* nlist_neighbor = &nlist_data_collection->compound_neighborlists[nlist_self->neighborcompound_ids[i]];

			for (int j = 0; j < nlist_neighbor->n_solvent_neighbors; j++) {
				int id_candidate = nlist_neighbor->neighborsolvent_ids[j];
				Float3 pos_candidate = nlist_data_collection->solvent_positions[id_candidate];
				applyHyperpos(&pos_self, &pos_candidate);
				double dist = (pos_self - pos_candidate).len();
				if (dist < CUTOFF) {
					NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_candidate];

					if (hashtable_solventneighbors.insert(id_candidate)) {
						nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::SOLVENT);
						nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
					}
				}
			}
		}


		for (int id_candidate = id_self + 1; id_candidate < simulation->n_compounds; id_candidate++) {	// For finding new nearby compounds, it is faster and simpler to just check all compounds, since there is so few
			Float3 pos_candidate = nlist_data_collection->compound_key_positions[id_candidate];
			applyHyperpos(&pos_self, &pos_candidate);
			double dist = (pos_self, pos_candidate).len();			
			if (dist < CUTOFF) {
				NeighborList* nlist_candidate = &nlist_data_collection->compound_neighborlists[id_candidate];
				if (hashtable_compoundneighbors.insert(id_candidate)) {
					nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::COMPOUND);
					nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
				}
			}
		}
	}

	// Finally add all solvent->solvent candidates
	for (int id_self = 0; id_self < simulation->n_solvents; id_self++) {
		NeighborList* nlist_self = &nlist_data_collection->solvent_neighborlists[id_self];
		HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, (int)nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
		Float3 pos_self = nlist_data_collection->solvent_positions[id_self];

		for (int i = 0; i < nlist_self->n_compound_neighbors; i++) {											// Use neighbor compounds to find new solvents
			NeighborList* nlist_neighbor = &nlist_data_collection->compound_neighborlists[nlist_self->neighborcompound_ids[i]];

			for (int j = 0; j < nlist_neighbor->n_solvent_neighbors; j++) {
				int id_candidate = nlist_neighbor->neighborsolvent_ids[j];
				if (id_candidate <= id_self) {
					continue;
				}


				Float3 pos_candidate = nlist_data_collection->solvent_positions[id_candidate];
				applyHyperpos(&pos_self, &pos_candidate);
				double dist = (pos_self - pos_candidate).len();
				if (dist < CUTOFF) {
					NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_candidate];

					if (hashtable_solventneighbors.insert(id_candidate)) {
						nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::SOLVENT);
						nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					}
				}
			}
		}
	}




	auto t1 = std::chrono::high_resolution_clock::now();
	*timing = (int) std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

	*finished = 1;
}





//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//


void Engine::step() {
	auto t0 = std::chrono::high_resolution_clock::now();
	forceKernel <<< simulation->box->n_compounds, THREADS_PER_COMPOUNDBLOCK >>> (simulation->box);
	solventForceKernel <<< BLOCKS_PER_SOLVENTKERNEL, THREADS_PER_SOLVENTBLOCK >>> (simulation->box);

	cudaDeviceSynchronize();
	auto t1 = std::chrono::high_resolution_clock::now();



	CompoundState* temp = simulation->box->compound_state_array;
	simulation->box->compound_state_array = simulation->box->compound_state_array_next;
	simulation->box->compound_state_array_next = temp;
	
	Solvent* temp_s = simulation->box->solvents;
	simulation->box->solvents = simulation->box->solvents_next;
	simulation->box->solvents_next = temp_s;

	//cudaMemcpy(simulation->box->compound_state_array, simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS, cudaMemcpyDeviceToDevice);	// Update all positions, after all forces have been calculated
	
	
	
	//cudaMemcpy(simulation->box->solvents, simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyDeviceToDevice);
	cudaDeviceSynchronize();
	auto t2 = std::chrono::high_resolution_clock::now();


	simulation->incStep();


	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
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




void __device__ __host__ applyHyperpos(Float3* static_particle, Float3* movable_particle) {
	//Float3 tmp = *movable_particle;
	for (int i = 0; i < 3; i++) {
		*movable_particle->placeAt(i) += BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) > BOX_LEN_HALF);
		*movable_particle->placeAt(i) -= BOX_LEN * ((static_particle->at(i) - movable_particle->at(i)) < -BOX_LEN_HALF);	// use at not X!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	}
	/*
	if ((*movable_particle - tmp).len() > 0.1 && !(tmp == Float3(0,0,0))) {
		static_particle->print('s');
		tmp.print('f');
		movable_particle->print('t');
	}*/
}

__device__ void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
	for (int dim = 0; dim < 3; dim++) {
		*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0);
		*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
	}
	
	/**current_position += Float3(BOX_LEN, BOX_LEN, BOX_LEN);
	current_position->x -= BOX_LEN * (current_position->x > BOX_LEN);
	current_position->y -= BOX_LEN * (current_position->y > BOX_LEN);
	current_position->z -= BOX_LEN * (current_position->z > BOX_LEN);
	//*current_position = current_position->elementwiseModulus(BOX_LEN);*/
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
	double force = 24.f * (epsilon / (dist)) * f6 * (1.f - 2.f * f6);

	double LJ_pot = 4.f * epsilon * (f12 - f6);
	Float3 force_unit_vector = (*pos1 - *pos0).norm();


	//printf("Mol %d force %f\n", blockIdx.x, (force_unit_vector * force).x);
	data_ptr[0] += LJ_pot *0.5 * 2;																	// WATCH OUT FOR 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	data_ptr[1] += force;
	data_ptr[2] = cudaMin(data_ptr[2], dist);

	Float3 force_vec = force_unit_vector * force;
	if (force_vec.x != force_vec.x) {
		printf("Force: %f\n", force);
		force_unit_vector.print('u');
	}

	//|| (threadIdx.x == 0 && blockIdx.x == 1)
	if (dist < 0.1f || abs(force) > 100e+6 ) {
	//if (abs(LJ_pot) > 7000) {
		printf("\nThread %d Block %d step %d dist %f force %f\n", threadIdx.x, blockIdx.x, (int) data_ptr[3], (*pos0 - *pos1).len(), force);
		//(*pos0 - *pos1).print('v');
		//(force_unit_vector * force).print('f');
		//pos0->print('0');
		//pos1->print('1');
		//printf("\n\n KILOFORCE! Block %d thread %d\n", blockIdx.x, threadIdx.x);
		//(force_unit_vector * force).print('f');
	}
	return force_unit_vector * force;	//J/mol*M	(M is direction)
}



constexpr double kb = 17.5 * 1e+6;		//	J/(mol*nm^2)
/*
__device__ Float3 calcPairbondForce(Float3* self_pos, Float3* other_pos, double* reference_dist, double* potE) {
	Float3 v = *self_pos - *other_pos;
	double dif = v.len() - *reference_dist;
	double invert = v.len() > *reference_dist ? -1 : 1;
	
	Float3 force = v.norm() * kb * (-dif);


	*potE += 0.5 * kb * (dif * dif) * 0.5;// *1 * 10e+4;

	if (abs(*potE) > 300000) {
		printf("\n");
		printf("threadid %d from: %f to: %f\n", threadIdx.x, *potE - 0.5 * kb * (dif * dif) * 0.5, *potE);

		self_pos->print('s');
		other_pos->print('o');
		printf("len %f\tdif: %f\t force %f\n", v.len(), dif, force.len());
	}

	// Logging stuff
	if (threadIdx.x == LOGTHREAD && blockIdx.x == LOGBLOCK) {
		//self_pos[63 - threadIdx.x].x = v.len();
		//self_pos[63 - threadIdx.x].y = *potE;
	}
	return force;
	//return v.norm() * (0.5 * kb * (dif * dif) * invert) *1.99 * 10e+4;
}
*/

__device__ void calcPairbondForces(Float3* pos_a, Float3* pos_b, double* reference_dist, Float3* results, double* potE) {
	Float3 vec1 = *pos_a - *pos_b;
	double error = vec1.len() - *reference_dist;

	*potE += 0.5 * kb * (error * error) * 0.5;// *1 * 10e+4;

	vec1 = vec1.norm();
	double force_scalar = -kb * error;
	results[0] = vec1 * force_scalar;
	results[1] = vec1 * force_scalar * -1;
}


constexpr double ktheta = 65 * 1e+3;	// J/mol
/*
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


	//if (threadIdx.x == 1 ) {
	if (force_scalar > 1000) {
		//printf("\n");
		//v1.print('1');
		//v2.print('2');
		//plane_normal.print('n');
		//force_direction.print('d');
		printf("dif: %f\tforce: %f\n", dif, force_scalar);
		//printf("\nAngle %f ref %f\n", angle, anglebond->reference_angle);
		//(force_direction * force_scalar).print('f');
	}

	//force_scalar = dif < 0 ? force_scalar * -1 : force_scalar;	// Invert force if we need to push outwards
	
	return force_direction * force_scalar;
}
*/
__device__ void calcAnglebondForces(Float3* pos_left, Float3* pos_middle, Float3* pos_right, double* reference_angle, Float3* results, double* potE) {
	Float3 v1 = *pos_left - *pos_middle;
	Float3 v2 = *pos_right - *pos_middle;
	Float3 normal1 = v1.cross(v2);
	Float3 normal2 = v2.cross(v1);

	Float3 inward_force_direction1 = (v1.cross(normal1)).norm();
	Float3 inward_force_direction2 = (v2.cross(normal2)).norm();

	double angle = Float3::getAngle(v1, v2);
	double error = angle - *reference_angle;

	//error = (error / (2.f * PI)) * 360.f;

	*potE += 0.5 * ktheta * error * error * 0.5;
	double force_scalar = ktheta * (error);

	results[0] = inward_force_direction1 * force_scalar;
	results[2] = inward_force_direction2 * force_scalar;
	results[1] = (results[0] + results[2]) * -1;

	//if (abs(error) > 0.0001)
		//printf("Error %f force %f \n", error, force_scalar);
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
/*
__device__ Float3 computeLJForces(Box* box, Compound* compound, SolventNeighborList* neighborlist, CompoundState* self_state, Float3* utility_buffer, double* data_ptr) {
	Float3 force(0, 0, 0);

	// ------------ Load and process neighbor molecule ------------ //
	for (int i = 0; i < 256; i += MAX_COMPOUND_PARTICLES) {
		int solvent_index = threadIdx.x + i;
		if ((solvent_index) < neighborlist->n_neighbors) {				// Dont think we have enough threads for this...
			utility_buffer[solvent_index] = box->solvents[neighborlist->neighborsolvent_indexes[solvent_index]].pos;
			//applyHyperpos(&self_state->positions[threadIdx.x], &utility_buffer[solvent_index]);
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
				//*
		}
		__syncthreads();
	}

	for (int i = 0; i < solventneighbor_list->n_neighbors; i++) {
		Float3 other_pos = box->solvents[solventneighbor_list->neighborsolvent_indexes[i]].pos;
		//applyHyperpos(self_pos, &other_pos);
		force += calcLJForce(self_pos, &other_pos, data_ptr);
	}
	return force;
}
*/
__device__ Float3 computeLJForces(Float3* self_pos, int n_particles, Float3* positions, double* data_ptr) {	// Assumes all positions are 
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		force += calcLJForce(self_pos, &positions[i], data_ptr);
	}
	return force;
}



__device__ Float3 computeSolventToSolventLJForces(Float3* self_pos, int self_index, int n_particles, Float3* positions, double* data_ptr) {	// Specific to solvent kernel
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		if (i != self_index) {
			Float3 hyperpos = positions[i];			// copy, DONT ref it as all threads will cann applyHyperpos
			applyHyperpos(self_pos, &hyperpos);
			force += calcLJForce(self_pos, &hyperpos, data_ptr);
		}		
	}
	return force;
}


__device__ Float3 computePairbondForces(Compound* compound, CompoundState* compound_state, Float3* utility_buffer, double* potE) {	// only works if n threads >= n bonds
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int i = 0; (i * blockDim.x) < compound->n_pairbonds; i++) {					
		PairBond* pb = nullptr;
		Float3 forces[2];
		int bond_index = threadIdx.x + i * blockDim.x;

		if (bond_index < compound->n_pairbonds) {
			pb = &compound->pairbonds[bond_index];

			calcPairbondForces(
				&compound_state->positions[pb->atom_indexes[0]],
				&compound_state->positions[pb->atom_indexes[1]],
				&pb->reference_dist,
				forces, potE
			);
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && pb != nullptr) {
				utility_buffer[pb->atom_indexes[0]] += forces[0];
				utility_buffer[pb->atom_indexes[1]] += forces[1];
			}
			__syncthreads();
		}
	}

	return utility_buffer[threadIdx.x];
}

__device__ Float3 computeAnglebondForces(Compound* compound, CompoundState* compound_state, Float3* utility_buffer, double* potE) {
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int i = 0; (i * blockDim.x) < compound->n_anglebonds; i++) {
		AngleBond* ab = nullptr;
		Float3 forces[3];
		int bond_index = threadIdx.x + i * blockDim.x;

		if (bond_index < compound->n_pairbonds) {
			ab = &compound->anglebonds[bond_index];

			calcAnglebondForces(
				&compound_state->positions[ab->atom_indexes[0]],
				&compound_state->positions[ab->atom_indexes[1]],
				&compound_state->positions[ab->atom_indexes[2]],
				&ab->reference_angle,
				forces, potE
			);
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && ab != nullptr) {
				utility_buffer[ab->atom_indexes[0]] += forces[0];
				utility_buffer[ab->atom_indexes[1]] += forces[1];
				utility_buffer[ab->atom_indexes[2]] += forces[2];
			}
			__syncthreads();
		}
	}

	return utility_buffer[threadIdx.x];
}


__device__ void integratePosition(Float3* pos, Float3* pos_tsub1, Float3* force, double mass, double dt) {
	Float3 temp = *pos;
	*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;	// no *0.5?
	*pos_tsub1 = temp;
	if (force->len() > 200e+6) {
		printf("\nThread %d blockId %d\tForce %f \tFrom %f %f %f\tTo %f %f %f\n", threadIdx.x, blockIdx.x, force->len(), pos_tsub1->x, pos_tsub1->y, pos_tsub1->z, pos->x, pos->y, pos->z);
		
	}
	//printf("force: %f\n", force->len());
}




// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//





__global__ void forceKernel(Box* box) {
	__shared__ Compound compound;
	__shared__ CompoundState compound_state;
	__shared__ NeighborList neighborlist;
	__shared__ Float3 utility_buffer[NEIGHBORLIST_MAX_SOLVENTS];							// waaaaay to biggg
	



	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
		compound_state.setMeta(compound.n_particles);
		neighborlist.loadMeta(&box->compound_neighborlists[blockIdx.x]);
	}
	__syncthreads();
	compound.loadData(&box->compounds[blockIdx.x]);
	compound_state.loadData(&box->compound_state_array[blockIdx.x]);
	neighborlist.loadData(&box->compound_neighborlists[blockIdx.x]);



	//if (!blockIdx.x && !threadIdx.x)
		//compound_state.positions[0].print('p');

	double data_ptr[4];
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999;
	data_ptr[3] = box->step + 1;

	Float3 force_bond(0.f);
	Float3 force_angle(0.f);
	Float3 force_LJ_com(0.f);
	Float3 force_LJ_sol(0.f);

	// ------------------------------------------------------------ Intramolecular Operations ------------------------------------------------------------ //
	applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);
	force_bond =  computePairbondForces(&compound, &compound_state, utility_buffer, &data_ptr[2]);
	if (force_bond.x != force_bond.x ) {
		force_bond.print('p');
		box->critical_error_encountered = 1;
	}
	force_angle = computeAnglebondForces(&compound, &compound_state, utility_buffer, &data_ptr[2]);
	if (force_angle.x != force_angle.x) {
		force_angle.print('a');
		box->critical_error_encountered = 1;
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //


	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {
		if (i == blockIdx.x)	// Only needed untill we have proper neighbor lists
			continue;
		int n_particles = box->compound_state_array[i].n_particles;

		if (threadIdx.x < n_particles) {	
			utility_buffer[threadIdx.x] = box->compound_state_array[i].positions[threadIdx.x];
			applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}
		__syncthreads();
		force_LJ_com += computeLJForces(&compound_state.positions[threadIdx.x], n_particles, utility_buffer, data_ptr);
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //




	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
	for (int i = threadIdx.x; i < neighborlist.n_solvent_neighbors; i += blockDim.x) {
		utility_buffer[i] = box->solvents[neighborlist.neighborsolvent_ids[i]].pos;
		applyHyperpos(&compound_state.positions[0], &utility_buffer[i]);
	}
	__syncthreads();
	if (threadIdx.x < compound.n_particles) {
		force_LJ_sol += computeLJForces(&compound_state.positions[threadIdx.x], neighborlist.n_solvent_neighbors, utility_buffer, data_ptr);
	}		
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //

	if (force_LJ_sol.x != force_LJ_sol.x) {
		force_LJ_com.print('1');
		force_LJ_sol.print('2');
		box->critical_error_encountered = 1;
	}


	/*if (force2 > 50e+6) {
		printf("Bond %f LJ %f total %f\n", force.len(), force2.len(), (force + force2).len());
	}
	force += force2;*/
	Float3 force = force_bond + force_angle + force_LJ_com + force_LJ_sol;

	if (force.x != force.x) {
		compound_state.positions[threadIdx.x].print('p');
	}



	
	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	if (threadIdx.x < compound.n_particles) {
		integratePosition(&compound_state.positions[threadIdx.x], &compound.particles[threadIdx.x].pos_tsub1, &force, compound.particles[threadIdx.x].mass, box->dt);
		box->compounds[blockIdx.x].particles[threadIdx.x].pos_tsub1 = compound.particles[threadIdx.x].pos_tsub1;
	}
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------- //
	




	// ------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------------------------- // 
	if (threadIdx.x == 0) {
		applyPBC(&compound_state.positions[threadIdx.x]);
		//if (blockIdx.x == 0 && ((compound_state.positions[0]-Float3(0.f)).len() > BOX_LEN))
		//compound_state.positions[0].print('s');
	}
	applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);	// So all particles follows p0
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
	


	
	// ------------------------------------ DATA LOG ------------------------------- //
	{
		/*
		if (threadIdx.x < 3) {	// TODO: UPDATE 3 TO PARTICLES_PER_COMPOUND
			box->potE_buffer[threadIdx.x + blockIdx.x * PARTICLES_PER_COMPOUND + box->step * box->total_particles] = data_ptr[0] + data_ptr[2];
			box->trajectory[threadIdx.x + blockIdx.x * PARTICLES_PER_COMPOUND + (box->step) * box->total_particles] = compound_state.positions[threadIdx.x];
		}
		__syncthreads();

		if (blockIdx.x == LOGBLOCK && threadIdx.x == LOGTHREAD && LOGTYPE == 1) {
			box->outdata[3 + box->step * 10] = data_ptr[0];	// LJ pot


			Float3 pos_temp = compound.particles[threadIdx.x].pos_tsub1;
			applyHyperpos(&compound_state.positions[threadIdx.x], &pos_temp);
			box->outdata[4 + box->step * 10] = (compound_state.positions[threadIdx.x] - pos_temp).len() / box->dt;

			box->outdata[5 + box->step * 10] = data_ptr[2];// closest particle
			box->outdata[6 + box->step * 10] = data_ptr[1];// force.len();
		}
		*/
		box->data_GAN[0 + threadIdx.x * 6 + box->step * MAX_COMPOUND_PARTICLES * 6] = compound_state.positions[threadIdx.x];
		box->data_GAN[1 + threadIdx.x * 6 + box->step * MAX_COMPOUND_PARTICLES * 6] = force_bond + force_angle;
		box->data_GAN[2 + threadIdx.x * 6 + box->step * MAX_COMPOUND_PARTICLES * 6] = force_LJ_com;
		box->data_GAN[3 + threadIdx.x * 6 + box->step * MAX_COMPOUND_PARTICLES * 6] = force_LJ_sol;
		if (threadIdx.x == 0) {
			//printf("Neighbors: %d\n", neighborlist.n_solvent_neighbors);
			//force_LJ_sol.print('s');
		}
			
	}
	
	// ----------------------------------------------------------------------------- //

	
	
	box->compound_state_array_next[blockIdx.x].loadData(&compound_state);
}












__global__ void solventForceKernel(Box* box) {
	__shared__ Float3 solvent_positions[THREADS_PER_SOLVENTBLOCK];
	__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];




	double data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;
	Float3 force(0,0,0);


	Solvent solvent = box->solvents[threadIdx.x];
	solvent_positions[threadIdx.x] = solvent.pos;




	// For compound LJ forces, maybe apply hyperpos to solvent, based on particle0 ? No that wont work, compound may be on both sides of box...
	//force += computeLJForces(box, &solvent.pos, &box->compound_neighborlist_array[MAX_COMPOUNDS + solvent_index], &box->solvent_neighborlist_array[MAX_COMPOUNDS + solvent_index], data_ptr, utility_buffer);

	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {
		int n_particles = box->compound_state_array[i].n_particles;

		// First all threads help loading the molecule
		if (threadIdx.x < n_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[i].positions[threadIdx.x];
		}
		__syncthreads();

		//Only do if mol is neighbor to the particular solvent
		applyHyperpos(&utility_buffer[0], &solvent_positions[threadIdx.x]);
		force += computeLJForces(&solvent_positions[threadIdx.x], n_particles, utility_buffer, data_ptr);

	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //
	if (threadIdx.x < box->n_solvents) {
		//applyHyperpos
		force += computeSolventToSolventLJForces(&solvent.pos, threadIdx.x, box->n_solvents, solvent_positions, data_ptr);
	}
	



	if (threadIdx.x < box->n_solvents) {
		integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, SOLVENT_MASS, box->dt);
	}

	if (threadIdx.x < box->n_solvents) {
		applyPBC(&solvent.pos);	// forcePositionToInsideBox	// This break the integration, as that doesn't accound for PBC in the vel conservation
	}





	// ------------------------------------ DATA LOG ------------------------------- //
	{/*
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
		}*/
	}

	// ----------------------------------------------------------------------------- //


	if (threadIdx.x < box->n_solvents) {
		box->solvents_next[threadIdx.x] = solvent;
	}
}


























/*	VELOCITY VERLET STORMER
__device__ void integratePosition(CompactParticle* particle, Float3* particle_pos, Float3* particle_force, double* dt) {
	*particle_pos = *particle_pos + (particle->vel + *particle_force * (0.5 / particle->mass) * *dt) * *dt;
}
__device__ void integrateVelocity(CompactParticle* particle, Float3* particle_force, double* dt) {
	particle->vel = particle->vel + (*particle_force + particle->force_prev) * (0.5 / particle->mass) * *dt;
}
*/