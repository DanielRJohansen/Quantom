#include "Engine.cuh"


__constant__ ForceField forcefield_device;
//__constant__ ForceField forcefield_nb_device;


Engine::Engine() {}
Engine::Engine(Simulation* simulation, ForceField forcefield_host) {
	LIMAENG::genericErrorCheck("Error before engine initialization.\n");

	this->simulation = simulation;

	nlist_data_collection = new NListDataCollection(simulation);





	int Ckernel_shared_mem = sizeof(Compound) + sizeof(CompoundState) + sizeof(NeighborList) + sizeof(Float3) * NEIGHBORLIST_MAX_SOLVENTS + sizeof(uint8_t) * NEIGHBORLIST_MAX_SOLVENTS;	
	int Skernel_shared_mem = sizeof(Float3) * MAX_COMPOUND_PARTICLES + sizeof(uint8_t) * MAX_COMPOUND_PARTICLES;
	printf("Compoundkernel shared mem. size: %d B\n", Ckernel_shared_mem);
	printf("Solventkernel shared mem. size: %d B\n", Skernel_shared_mem);



	//ForceField forcefield_host = FFM.getForcefield();
	//ForceField forcefield_host = FFM.getNBForcefield();
	cudaMemcpyToSymbol(forcefield_device, &forcefield_host, sizeof(ForceField), 0, cudaMemcpyHostToDevice);

	//ForceField forcefield_nb_host = FFM.getNBForcefield();
	//cudaMemcpyToSymbol(forcefield_nb_device, &forcefield_nb_host, sizeof(ForceField), 0, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	LIMAENG::genericErrorCheck("Error while moving forcefield to device\n");


	handleNLISTS(simulation, false, true);				// Fix neighborlists before running

	printf("Engine ready\n\n\n");
}





void Engine::deviceMaster() {
	LIMAENG::genericErrorCheck("Error before step!");
	step();
}

void Engine::hostMaster() {						// This is and MUST ALWAYS be called after the deviceMaster, and AFTER intStep()!
	auto t0 = std::chrono::high_resolution_clock::now();
	if ((simulation->getStep() % STEPS_PER_LOGTRANSFER) == 0) {
		offloadLoggingData();
		offloadPositionData();
		if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 0 && APPLY_THERMOSTAT) {
			applyThermostat();
		}
	}
	if ((simulation->getStep() % STEPS_PER_TRAINDATATRANSFER) == 0) {
		offloadTrainData();
	}

	
	handleNLISTS(simulation);
	

	if ((simulation->getStep() % STEPS_PER_THERMOSTAT) == 1) {	// So this runs 1 step AFTER applyThermostat
		simulation->box->thermostat_scalar = 1.f;
	}

	auto t1 = std::chrono::high_resolution_clock::now();

	int cpu_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	timings = timings + Int3(0,0,cpu_duration);
}









//--------------------------------------------------------------------------	CPU workload --------------------------------------------------------------//


void Engine::handleNLISTS(Simulation* simulation, bool async, bool force_update) {



	if ((simulation->getStep() - prev_nlist_update_step >= STEPS_PER_NLIST_UPDATE) || force_update) {
		offloadPositionDataNLIST(simulation);
		// Lots of waiting time here...
		cudaDeviceSynchronize();


		if (async) {
			std::thread nlist_worker(Engine::updateNeighborLists, simulation, nlist_data_collection, &updated_neighborlists_ready, &timings.z);
			nlist_worker.detach();
		}
		else {
			updateNeighborLists(simulation, nlist_data_collection, &updated_neighborlists_ready, &timings.z);
		}

		prev_nlist_update_step = simulation->getStep();
	}

	if (updated_neighborlists_ready) {
		pushNlistsToDevice();
	}

	if (force_update) {
		Int3 n_data(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
		//Int3 after(nlist_data_collection->solvent_neighborlists[193].n_compound_neighbors, nlist_data_collection->solvent_neighborlists[193].n_solvent_neighbors, 0);

		printf("\nEntity neighbors: %d %d\n", n_data.x, n_data.y);
	}
}

void Engine::offloadPositionDataNLIST(Simulation* simulation) {
	//cudaMemcpyAsync(compoundstatearray_host, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->box->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(nlist_data_collection->compoundstates, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(nlist_data_collection->solvents, simulation->box->solvents, sizeof(Solvent) * simulation->n_solvents, cudaMemcpyDeviceToHost);
}

void Engine::pushNlistsToDevice() {
	cudaMemcpy(simulation->box->compound_neighborlists, nlist_data_collection->compound_neighborlists, sizeof(NeighborList) * simulation->n_compounds, cudaMemcpyHostToDevice);
	cudaMemcpy(simulation->box->solvent_neighborlists, nlist_data_collection->solvent_neighborlists, sizeof(NeighborList) * simulation->n_solvents, cudaMemcpyHostToDevice);
	updated_neighborlists_ready = 0;
}

/*
void Engine::offloadPositionData(Simulation* simulation) {
	//cudaMemcpyAsync(compoundstatearray_host, simulation->box->compound_state_array, sizeof(CompoundState) * simulation->box->n_compounds, cudaMemcpyDeviceToHost);

	const int step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) * simulation->total_particles_upperbound;	// Tongue in cheek here, i think this is correct...
	cudaMemcpy(&nlist_data_collection->compoundstates[step_offset], simulation->box->compound_state_array, sizeof(CompoundState) * simulation->n_compounds, cudaMemcpyDeviceToHost);
	cudaMemcpy(&nlist_data_collection->solvents[step_offset], simulation->box->solvents, sizeof(Solvent) * simulation->n_solvents, cudaMemcpyDeviceToHost);
}
*/


/*
bool Engine::neighborWithinCutoff(Float3* pos_a, Float3* pos_b) {
	Float3 pos_b_temp = *pos_b;
	LIMAENG::applyHyperpos(pos_a, &pos_b_temp);
	double dist = (*pos_a - pos_b_temp).len();
	return (dist < CUTOFF);
}*/

bool Engine::neighborWithinCutoff(Float3* pos_a, Float3* pos_b, float cutoff_offset) {		// This is used for compounds with a confining_particle_sphere from key_particle BEFORE CUTOFF begins
	Float3 pos_b_temp = *pos_b;
	LIMAENG::applyHyperpos(pos_a, &pos_b_temp);
	double dist = (*pos_a - pos_b_temp).len();
	return (dist < (CUTOFF + cutoff_offset));
}



void Engine::cullDistantNeighbors(Simulation* simulation, NListDataCollection* nlist_data_collection) {	// Calling with nlist avoids writing function for both solvent and compound
	for (int id_self = 0; id_self < nlist_data_collection->n_compounds; id_self++) {
		NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
		float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;



		for (int j = 0; j < nlist_self->n_compound_neighbors; j++) {		// Cull compound-compound
			int id_neighbor = nlist_self->neighborcompound_ids[j];
			NeighborList* nlist_neighbor = &nlist_data_collection->compound_neighborlists[id_neighbor];
			float cutoff_add_neighbor = simulation->compounds_host[id_neighbor].confining_particle_sphere;

			if (id_self < id_neighbor) {
				if (!neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->compound_key_positions[id_neighbor], cutoff_add_self + cutoff_add_neighbor + CUTOFF)) {
					nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::COMPOUND);
					nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);				
					j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
				}
			}
		}


		for (int j = 0; j < nlist_self->n_solvent_neighbors; j++) {			// Cull compound-solvent
			int id_neighbor = nlist_self->neighborsolvent_ids[j];
			NeighborList* nlist_neighbor = &nlist_data_collection->solvent_neighborlists[id_neighbor];

			//printf("Dist: %f\n", (nlist_data_collection->compound_key_positions[id_self] - nlist_data_collection->solvent_positions[id_neighbor]).len());
			if (!neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->solvent_positions[id_neighbor], cutoff_add_self + CUTOFF)) {
				nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
				j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
			}
		}
	}


	for (int id_self = 0; id_self < nlist_data_collection->n_solvents; id_self++) {																// Cull solvent-solvent
		NeighborList* nlist_self = &nlist_data_collection->solvent_neighborlists[id_self];

		for (int j = 0; j < nlist_self->n_solvent_neighbors; j++) {			/// NOT FINISHED HERE
			int id_neighbor = nlist_self->neighborsolvent_ids[j];
			NeighborList* nlist_neighbor = &nlist_data_collection->solvent_neighborlists[id_neighbor];

			if (!neighborWithinCutoff(&nlist_data_collection->solvent_positions[id_self], &nlist_data_collection->solvent_positions[id_neighbor], CUTOFF)) {
				nlist_self->removeId(id_neighbor, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				nlist_neighbor->removeId(id_self, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				j--;	// Decrement, as the removeId puts the last element at the current and now vacant spot.
			}
			
		}
	}
}

void Engine::updateNeighborLists(Simulation* simulation, NListDataCollection* nlist_data_collection, volatile bool* finished, int* timing) {	// This is a thread worker-function, so it can't own the object, thus i pass a ref to the engine object..
	auto t0 = std::chrono::high_resolution_clock::now();
	Int3 before(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
	//Int3 before(nlist_data_collection->solvent_neighborlists[193].n_compound_neighbors, nlist_data_collection->solvent_neighborlists[193].n_solvent_neighbors, 0);
	//nlist_data_collection->preparePositionData();
	nlist_data_collection->preparePositionData(simulation->compounds_host);		// Makes key positions addressable in arrays: compound_key_positions and solvent_positions
	// First do culling of neighbors that has left CUTOFF
	cullDistantNeighbors(simulation, nlist_data_collection);


	// First add compound->solvent, compound->compound
	for (int id_self = 0; id_self < simulation->n_compounds; id_self++) {										
		NeighborList* nlist_self = &nlist_data_collection->compound_neighborlists[id_self];
		HashTable hashtable_compoundneighbors(nlist_self->neighborcompound_ids, nlist_self->n_compound_neighbors, NEIGHBORLIST_MAX_COMPOUNDS * 2);
		HashTable hashtable_solventneighbors(nlist_self->neighborsolvent_ids, nlist_self->n_solvent_neighbors, NEIGHBORLIST_MAX_SOLVENTS * 2);
		float cutoff_add_self = simulation->compounds_host[id_self].confining_particle_sphere;

		//printf("\nSize to hashtable: %d\n", nlist_self->n_solvent_neighbors);
		//printf("cutoff %f\n", cutoff_offset);

		// Go through all solvents in box!
		for (int id_candidate = 0; id_candidate < simulation->n_solvents; id_candidate++) {
			NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_candidate];
			if (neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->solvent_positions[id_candidate], cutoff_add_self + CUTOFF)) {
				if (hashtable_solventneighbors.insert(id_candidate)) {
					nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::COMPOUND);
				}
			}

		}

		// Go through all compounds in box, with higher ID than self!
		for (int id_candidate = id_self + 1; id_candidate < simulation->n_compounds; id_candidate++) {	// For finding new nearby compounds, it is faster and simpler to just check all compounds, since there are so few
			NeighborList* nlist_candidate = &nlist_data_collection->compound_neighborlists[id_candidate];
			float cutoff_add_candidate = simulation->compounds_host[id_self].confining_particle_sphere;

			//printf("Distance to neighbor compound: %f\n", (nlist_data_collection->compound_key_positions[id_self] - nlist_data_collection->compound_key_positions[id_candidate]).len());
			if (neighborWithinCutoff(&nlist_data_collection->compound_key_positions[id_self], &nlist_data_collection->compound_key_positions[id_candidate], cutoff_add_self + cutoff_add_candidate + CUTOFF)) {
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



		for (int id_candidate = id_self + 1; id_candidate < simulation->n_solvents; id_candidate++) {
			NeighborList* nlist_candidate = &nlist_data_collection->solvent_neighborlists[id_candidate];
			if (neighborWithinCutoff(&nlist_data_collection->solvent_positions[id_self], &nlist_data_collection->solvent_positions[id_candidate], CUTOFF)) {
				if (hashtable_solventneighbors.insert(id_candidate)) {
					nlist_self->addId(id_candidate, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					nlist_candidate->addId(id_self, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					if (id_self == 193) {
						//printf("Adding id %d with dist %f\n", id_candidate, (nlist_data_collection->solvent_positions[id_self] - nlist_data_collection->solvent_positions[id_candidate]).len());
					}
				}
			}
			
		}
	}

	Int3 after(nlist_data_collection->compound_neighborlists[0].n_compound_neighbors, nlist_data_collection->compound_neighborlists[0].n_solvent_neighbors, 0);
	//Int3 after(nlist_data_collection->solvent_neighborlists[193].n_compound_neighbors, nlist_data_collection->solvent_neighborlists[193].n_solvent_neighbors, 0);

	//printf("\nEntity went from %d %d neighbors to %d %d\n", before.x, before.y, after.x, after.y);
	
	auto t1 = std::chrono::high_resolution_clock::now();
	*timing = (int) std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

	*finished = 1;		// Thread terminates here!
	
}


void Engine::offloadLoggingData() {
	const int step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) * simulation->total_particles_upperbound;	// Tongue in cheek here, i think this is correct...
	cudaMemcpy(&simulation->potE_buffer[step_offset], simulation->box->potE_buffer, sizeof(double) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
}

void Engine::offloadPositionData() {
	const int step_offset = (simulation->getStep() - STEPS_PER_LOGTRANSFER) * simulation->total_particles_upperbound;	// Tongue in cheek here, i think this is correct...
	cudaMemcpy(&simulation->traj_buffer[step_offset], simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER, cudaMemcpyDeviceToHost);
}

void Engine::offloadTrainData() {
	const int step_offset = (simulation->getStep() - STEPS_PER_TRAINDATATRANSFER) * MAX_COMPOUND_PARTICLES * 6;	// fix max_compound to the actual count save LOTS of space!. Might need a file in simout that specifies cnt for loading in other programs...
	cudaMemcpy(&simulation->traindata_buffer[step_offset], simulation->box->data_GAN, sizeof(Float3) * MAX_COMPOUND_PARTICLES * 6 * STEPS_PER_TRAINDATATRANSFER, cudaMemcpyDeviceToHost);
	LIMAENG::genericErrorCheck("Cuda error during traindata offloading\n");
}



float Engine::getBoxTemperature() {
	
	const int step = simulation->getStep() - 1;

	long double temp_sum = 0;				// [k]

	long double kinE_sum = 0;

	const int step_offset_a = step * simulation->total_particles_upperbound;
	const int step_offset_b = (step - 2) * simulation->total_particles_upperbound;
	const int solvent_offset = MAX_COMPOUND_PARTICLES;



	int particles_total = simulation->n_solvents;

	for (int c = 0; c < simulation->n_compounds; c++) {
		for (int i = 0; i < simulation->compounds_host[c].n_particles; i++) {	// i gotta move this somewhere else....

			Float3 posa = simulation->traj_buffer[i + step_offset_a];
			Float3 posb = simulation->traj_buffer[i + step_offset_b];
			float kinE = LIMAENG::calcKineticEnergy(&posa, &posb, forcefield_device.particle_parameters[simulation->box->compounds[c].atom_types[i]].mass, simulation->dt);			// Doesnt work, use forcefield_host!!
			//float kinE = LIMAENG::calcKineticEnergy(&posa, &posb, forcefield_device.particle_parameters[simulation->box->compounds[c].atom_types[i]].mass, simulation->dt);			// Doesnt work, use forcefield_host!!
			kinE_sum += kinE;
			particles_total++;
		}
	}
	
	for (int i = 0; i < simulation->n_solvents; i++) {
		Float3 posa = simulation->traj_buffer[i + solvent_offset + step_offset_a];
		Float3 posb = simulation->traj_buffer[i + solvent_offset + step_offset_b];
		float kinE = LIMAENG::calcKineticEnergy(&posa, &posb, SOLVENT_MASS, simulation->dt);
		kinE_sum += kinE;
	}

	double avg_kinE = kinE_sum / (long double)particles_total;
	float temperature = avg_kinE * 2 / (3 * 8.3145);
	simulation->temperature_buffer[step / STEPS_PER_THERMOSTAT] = temperature;
	//printf("\nTemp: %f\n", temperature);
	return temperature;
}

void Engine::applyThermostat() {
	const float max_temp = 300.f;				// [k]
	float temp = getBoxTemperature();
	//printf("\n %d Temperature: %f\n", (simulation->getStep()-1) / STEPS_PER_THERMOSTAT, temp);
	if (temp > 10) {
		simulation->box->thermostat_scalar = max_temp / temp;
		//printf("Scalar: %f\n", simulation->box->thermostat_scalar);
	}
	else {
		printf("Critically low temperature encountered\n");
		simulation->box->critical_error_encountered = true;
	}
}



//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//


void Engine::step() {
	auto t0 = std::chrono::high_resolution_clock::now();

	compoundBridgeKernel <<< simulation->box->bridge_bundle->n_bridges, MAX_PARTICLES_IN_BRIDGE >>> (simulation->box);	// Must come before forceKernel()
	cudaDeviceSynchronize();
	//printf("\n\n");
	forceKernel <<< simulation->box->n_compounds, THREADS_PER_COMPOUNDBLOCK >>> (simulation->box);
	solventForceKernel <<< BLOCKS_PER_SOLVENTKERNEL, THREADS_PER_SOLVENTBLOCK >>> (simulation->box);
	cudaDeviceSynchronize();
	//printf("\n\n");

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
	LIMAENG::genericErrorCheck("Error during step or state_transfer\n");		// Temp, we want to do host stuff while waiting for async GPU operations...
	auto t2 = std::chrono::high_resolution_clock::now();


	simulation->incStep();


	int force_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
	int copy_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
	timings = timings + Int3(force_duration, copy_duration, 0);
}





// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//



__device__ void applyPBC(Float3* current_position) {	// Only changes position if position is outside of box;
	for (int dim = 0; dim < 3; dim++) {
		*current_position->placeAt(dim) += BOX_LEN * (current_position->at(dim) < 0);
		*current_position->placeAt(dim) -= BOX_LEN * (current_position->at(dim) > BOX_LEN);
	}
}



/*
constexpr double sigma = 0.3923f;										// nm, point at which d/dE or force is 0.
constexpr double epsilon = 0.5986 * 1'000.f;							// J/mol
*/

__device__ Float3 calcLJForce(Float3* pos0, Float3* pos1, float* data_ptr, float* potE, float sigma, float epsilon, bool verbose=false) {
	// Calculates LJ force on p0	//
	// input positions in cartesian coordinates [nm]
	// sigma [nm]
	// epsilon
	//	Returns force in J/mol*M		?????????????!?!?//

	double dist = (*pos0 - *pos1).len();								// [nm]
	

	double fraction = sigma / dist;										// [nm/nm], so unitless
	double f2 = fraction * fraction;
	double f6 = f2 * f2 * f2;

	double force = 24.f * (epsilon / (dist)) * f6 * (1.f - 2.f * f6) * 1e-9;	//	1e+9 * [kg*m/(mol*s^2)] * 1e-9 => [kg*nm/(mol * ns^2)]										wrong: [J/mol]/[nm] [N/mol] 

	double LJ_pot = 4.f * epsilon * (f6 * f6 - f6);						// [J/mol]
	Float3 force_unit_vector = (*pos1 - *pos0).norm();

	*potE += LJ_pot * 0.5;												// Log this value for energy monitoring purposes

#ifdef LIMA_VERBOSE
	if (verbose) {
		//printf("Dist %f       force %f      sigma: %f         epsilon %f\n", (float)dist, (float)force, sigma, epsilon);
	}
	if (threadIdx.x == 0 && blockIdx.x == 0) {
		//printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);
		//printf("Force %f %f %f\n", force_unit_vector.x * force, force_unit_vector.y * force, force_unit_vector.z * force);
	}


	{
		Float3 force_vec = force_unit_vector * force;
		if (force_vec.x != force_vec.x) {
			//printf("Force: %f\n", force);
			//force_unit_vector.print('u');
			printf("Thread %d Block %d self %f %f %f other %f %f %f\n", threadIdx.x, blockIdx.x, pos0->x, pos0->y, pos0->z, pos1->x, pos1->y, pos1->z);

		}
		if (dist < 0.1f) {
			//printf("\nThread %d Block %d step %d dist %f force %f\n", threadIdx.x, blockIdx.x, (int)data_ptr[3], (*pos0 - *pos1).len(), force);
			//printf("Self %f %f %f -> %f %f %f\n", data_ptr[0], data_ptr[1], data_ptr[2], pos0->x, pos0->y, pos0->z);
		}
	}
#endif

	return force_unit_vector * force;									//N/mol*v	(v is direction)
}

constexpr double kb = 17.5 * 1e+6;		//	J/(mol*nm^2)
__device__ void calcPairbondForces(Float3* pos_a, Float3* pos_b, float* reference_dist, Float3* results, float* potE) {
	// Calculates bond force on both particles					//
	// Calculates forces as J/mol*M								//

	Float3 difference = *pos_a - *pos_b;						//	[nm]
	double error = difference.len() - *reference_dist;			//	[nm]

	*potE += 0.5 * kb * (error * error) * 0.5;					// [J/mol]

	difference = difference.norm();								// dif_unit_vec, but shares register with dif
	double force_scalar = -kb * error;							// [N/mol] directionless, yes?

	results[0] = difference * force_scalar;
	results[1] = difference * force_scalar * -1;

#ifdef LIMA_VERBOSE
	if (error > 12.7f) {
		printf("thread %d  error %f ref %f force %f\n", threadIdx.x, error, *reference_dist, force_scalar);
		pos_a->print('a');
		pos_b->print('b');
	}
#endif
		
}


constexpr double ktheta = 65 * 1e+3;	// J/mol
__device__ void calcAnglebondForces(Float3* pos_left, Float3* pos_middle, Float3* pos_right, float* reference_angle, Float3* results, float* potE) {
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
}

__device__ void calcDihedralbondForces(Float3* pos_left, Float3* pos_lm, Float3* pos_rm, Float3* pos_right, DihedralBond* dihedral, Float3* results, float* potE) {
	Float3 normal1 = (*pos_left-*pos_lm).cross((*pos_rm-*pos_lm));
	Float3 normal2 = (*pos_lm-*pos_rm).cross((*pos_right-*pos_rm));


	float torsion = Float3::getAngle(normal1, normal2);
	float error = (torsion - dihedral->phi_0);	// *360.f / 2.f * PI;	// Check whether k_phi is in radians or degrees

	//error = (error / (2.f * PI)) * 360.f;

	*potE += 0.5 * ktheta * error * error * 0.5;
	double force_scalar = dihedral->k_phi * (error);


	results[0] = normal1 * force_scalar * -1.f;
	results[3] = normal2 * force_scalar;
	// Not sure about the final two forces, for now we'll jsut split the sum of opposite forces between them.
	results[1] = (results[0] + results[3]) * -1.f * 0.5;
	results[2] = (results[0] + results[3]) * -1.f * 0.5;

	/*
	if (blockIdx.x == 0 && threadIdx.x == 30) {
		printf("Torsion %f error %f force_scalar %f, phi0 %f kphi %f\n", torsion, error, force_scalar, dihedral->phi_0, dihedral->k_phi);
		results[0].print('0');
		results[1].print('1');
		results[3].print('3');
	}
	*/
}


__device__ Float3 computeLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self, uint8_t* atomtypes_others) {	// Assumes all positions are 
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		force += calcLJForce(self_pos, &positions[i], data_ptr, potE_sum, 
			forcefield_device.particle_parameters[atomtype_self].sigma + forcefield_device.particle_parameters[atomtypes_others[i]].sigma,
			forcefield_device.particle_parameters[atomtype_self].epsilon + forcefield_device.particle_parameters[atomtypes_others[i]].epsilon
		);
	}
	return force;
}

__device__ Float3 computeSolventToSolventLJForces(Float3* self_pos, NeighborList* nlist, Solvent* solvents, float* data_ptr, float* potE_sum) {	// Specific to solvent kernel
	Float3 force(0.f);
	for (int i = 0; i < nlist->n_solvent_neighbors; i++) {
		Solvent neighbor = solvents[nlist->neighborsolvent_ids[i]];			
		LIMAENG::applyHyperpos(&neighbor.pos, self_pos);
		force += calcLJForce(self_pos, &neighbor.pos, data_ptr, potE_sum, forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma * 2.f, forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon * 2.f);
		/*if (abs(force.len()) > 1000000) {
			(*self_pos - neighbor.pos).print('d');
			printf("F %f %f %f n %d  solvent_id %d neighbor_id %d\n", force.x, force.y, force.z, nlist->n_solvent_neighbors, threadIdx.x + blockIdx.x*blockDim.x, nlist->neighborsolvent_ids[i]);
		}*/
	}
	return force;
}
__device__ Float3 computeSolventToCompoundLJForces(Float3* self_pos, int n_particles, Float3* positions, float* data_ptr, float* potE_sum, uint8_t atomtype_self) {	// Specific to solvent kernel
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		force += calcLJForce(self_pos, &positions[i], data_ptr, potE_sum,
			forcefield_device.particle_parameters[atomtype_self].sigma + forcefield_device.particle_parameters[ATOMTYPE_SOL].sigma,
			forcefield_device.particle_parameters[atomtype_self].epsilon + forcefield_device.particle_parameters[ATOMTYPE_SOL].epsilon
		);
	}
	return force;
}


/*
__device__ Float3 computeCompoundToSolventLJForces(Float3* self_pos, int n_particles, Float3* positions, double* data_ptr, double* potE_sum) {	// Specific to solvent kernel
	Float3 force(0, 0, 0);
	for (int i = 0; i < n_particles; i++) {
		Float3 hyperpos = positions[i];			// copy, DONT ref it as all threads will cann applyHyperpos
		LIMAENG::applyHyperpos(self_pos, &hyperpos);
		force += calcLJForce(self_pos, &hyperpos, data_ptr, potE_sum);	
	}
	return force;
}
*/
template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computePairbondForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {	// only works if n threads >= n bonds
	utility_buffer[threadIdx.x] = Float3(0.f);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_singlebonds; bond_offset++) {
		PairBond* pb = nullptr;
		Float3 forces[2];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_singlebonds) {
			pb = &entity->singlebonds[bond_index];

			calcPairbondForces(
				&positions[pb->atom_indexes[0]],
				&positions[pb->atom_indexes[1]],
				//&pb->reference_dist,
				&pb->b0,
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

template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computeAnglebondForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_anglebonds; bond_offset++) {
		AngleBond* ab = nullptr;
		Float3 forces[3];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_anglebonds) {
			ab = &entity->anglebonds[bond_index];

			calcAnglebondForces(
				&positions[ab->atom_indexes[0]],
				&positions[ab->atom_indexes[1]],
				&positions[ab->atom_indexes[2]],
				//&ab->reference_angle,
				&ab->theta_0,
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

template <typename T>	// Can either be Compound or CompoundBridgeCompact
__device__ Float3 computeDihedralForces(T* entity, Float3* positions, Float3* utility_buffer, float* potE) {
	utility_buffer[threadIdx.x] = Float3(0, 0, 0);
	for (int bond_offset = 0; (bond_offset * blockDim.x) < entity->n_anglebonds; bond_offset++) {
		DihedralBond* db = nullptr;
		Float3 forces[4];
		int bond_index = threadIdx.x + bond_offset * blockDim.x;

		if (bond_index < entity->n_anglebonds) {
			db = &entity->dihedrals[bond_index];

			calcDihedralbondForces(
				&positions[db->atom_indexes[0]],
				&positions[db->atom_indexes[1]],
				&positions[db->atom_indexes[2]],
				&positions[db->atom_indexes[3]],
				db,
				forces, 
				potE
			);
		}


		for (int i = 0; i < blockDim.x; i++) {
			if (threadIdx.x == i && db != nullptr) {
				utility_buffer[db->atom_indexes[0]] += forces[0];
				utility_buffer[db->atom_indexes[1]] += forces[1];
				utility_buffer[db->atom_indexes[2]] += forces[2];
				utility_buffer[db->atom_indexes[3]] += forces[3];
			}
			__syncthreads();
		}
	}

	return utility_buffer[threadIdx.x];
}


__device__ void integratePosition(Float3* pos, Float3* pos_tsub1, Float3* force, const double mass, const double dt, float* thermostat_scalar, int p_index, bool verbose=false) {
	// Force is in ??Newton, [kg * nm /(mol*ns^2)] //


	Float3 temp = *pos;
	LIMAENG::applyHyperpos(pos, pos_tsub1);
	//*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;	// no *0.5?	// [nm] - [nm] + [kg*m*s^-2]/[kg] * [ns]^2
	*pos = *pos * 2 - *pos_tsub1 + *force * (1.f / mass) * dt * dt;		// [nm] - [nm] + [kg/mol*m*/s^2]/[kg/mol] * [s]^2 * (1e-9)^2	=> [nm]-[nm]+[]
	*pos_tsub1 = temp;

	Float3 delta_pos = *pos - *pos_tsub1;
	*pos = *pos_tsub1 + delta_pos * *thermostat_scalar;
#ifdef LIMA_VERBOSE
	if (force->len() > 200e+6) {
		printf("\nP_index %d Thread %d blockId %d\tForce %f \tFrom %f %f %f\tTo %f %f %f\n", p_index, threadIdx.x, blockIdx.x, force->len(), pos_tsub1->x, pos_tsub1->y, pos_tsub1->z, pos->x, pos->y, pos->z);		
	}
	if (verbose) {
		printf("Mass %f Force %f %f %f\n", mass, force->x, force->y, force->z);
	}
#endif
}




// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//





__global__ void forceKernel(Box* box) {
	__shared__ Compound compound;
	__shared__ CompoundState compound_state;
	__shared__ NeighborList neighborlist;
	__shared__ Float3 utility_buffer[NEIGHBORLIST_MAX_SOLVENTS];							// waaaaay too biggg
	__shared__ uint8_t utility_buffer_small[NEIGHBORLIST_MAX_SOLVENTS];


	if (threadIdx.x == 0) {
		compound.loadMeta(&box->compounds[blockIdx.x]);
		compound_state.setMeta(compound.n_particles);
		neighborlist.loadMeta(&box->compound_neighborlists[blockIdx.x]);
	}
	__syncthreads();
	compound.loadData(&box->compounds[blockIdx.x]);
	compound_state.loadData(&box->compound_state_array[blockIdx.x]);
	neighborlist.loadData(&box->compound_neighborlists[blockIdx.x]);
	__syncthreads();



	float potE_sum = 0;
	float data_ptr[4];
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999;
	data_ptr[3] = box->step + 1;

	//Float3 force(0.f);
	Float3 force = compound.forces[threadIdx.x];
	// ------------------------------------------------------------ Intramolecular Operations ------------------------------------------------------------ //
	{

		LIMAENG::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);
		__syncthreads();	// Dunno if necessary
		force += computePairbondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		force += computeDihedralForces(&compound, compound_state.positions, utility_buffer, &potE_sum);
		//box->critical_error_encountered = true;
		for (int i = 0; i < compound.n_particles; i++) {
			if (i != threadIdx.x) {
				force += calcLJForce(&compound_state.positions[threadIdx.x], &compound_state.positions[i], data_ptr, &potE_sum,
					(forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].sigma + forcefield_device.particle_parameters[compound.atom_types[i]].sigma),
					(forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].epsilon + forcefield_device.particle_parameters[compound.atom_types[i]].epsilon),
					true
					);
			}
		}
	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------- //


	// --------------------------------------------------------------- Intermolecular forces --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {

		if (i == blockIdx.x)	// Only needed untill we have proper neighbor lists
			continue;
		int n_particles = box->compound_state_array[i].n_particles;

		if (threadIdx.x < n_particles) {	
			utility_buffer[threadIdx.x] = box->compound_state_array[i].positions[threadIdx.x];
			utility_buffer_small[threadIdx.x] = box->compounds[i].atom_types[threadIdx.x];
			LIMAENG::applyHyperpos(&compound_state.positions[0], &utility_buffer[threadIdx.x]);
		}
		__syncthreads();
//		force += computeLJForces(&compound_state.positions[threadIdx.x], n_particles, utility_buffer, data_ptr, &potE_sum, compound.atom_types[threadIdx.x], utility_buffer_small);
	}
	// ------------------------------------------------------------------------------------------------------------------------------------------------------ //




	// --------------------------------------------------------------- Solvation forces --------------------------------------------------------------- //
	for (int i = threadIdx.x; i < neighborlist.n_solvent_neighbors; i += blockDim.x) {
		utility_buffer[i] = box->solvents[neighborlist.neighborsolvent_ids[i]].pos;
		LIMAENG::applyHyperpos(&compound_state.positions[0], &utility_buffer[i]);
	}
	__syncthreads();
	if (threadIdx.x < compound.n_particles) {
		force += computeSolventToCompoundLJForces(&compound_state.positions[threadIdx.x], neighborlist.n_solvent_neighbors, utility_buffer, data_ptr, &potE_sum, compound.atom_types[threadIdx.x]);
	}		
	// ------------------------------------------------------------------------------------------------------------------------------------------------ //


	
	// ------------------------------------------------------------ Integration ------------------------------------------------------------ //
	if (threadIdx.x < compound.n_particles) {
		int p_index = threadIdx.x;
		//integratePosition(&compound_state.positions[threadIdx.x], &compound.particles[threadIdx.x].pos_tsub1, &force, compound.particles[threadIdx.x].mass, box->dt, &box->thermostat_scalar, p_index );

		//integratePosition(&compound_state.positions[threadIdx.x], &compound.particles[threadIdx.x].pos_tsub1, &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, &box->thermostat_scalar, p_index);
		integratePosition(&compound_state.positions[threadIdx.x], &compound.prev_positions[threadIdx.x], &force, forcefield_device.particle_parameters[compound.atom_types[threadIdx.x]].mass, box->dt, &box->thermostat_scalar, p_index);
		
		//box->compounds[blockIdx.x].particles[threadIdx.x].pos_tsub1 = compound.particles[threadIdx.x].pos_tsub1;
		box->compounds[blockIdx.x].prev_positions[threadIdx.x] = compound.prev_positions[threadIdx.x];
	}
	__syncthreads();
	// ------------------------------------------------------------------------------------------------------------------------------------- //





	// ------------------------------------ PERIODIC BOUNDARY CONDITION ------------------------------------------------------------------------------------------------- // 
	if (threadIdx.x == 0) {
		applyPBC(&compound_state.positions[threadIdx.x]);
	}
	LIMAENG::applyHyperpos(&compound_state.positions[0], &compound_state.positions[threadIdx.x]);	// So all particles follows p0
	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------ //
	


	
	// ------------------------------------ DATA LOG ------------------------------- //
	{
		if (threadIdx.x < compound.n_particles) {
			int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;
			box->potE_buffer[threadIdx.x + blockIdx.x * MAX_COMPOUND_PARTICLES + step_offset] = potE_sum;//data_ptr[0] + data_ptr[2];
			box->traj_buffer[threadIdx.x + blockIdx.x * MAX_COMPOUND_PARTICLES + step_offset] = compound_state.positions[threadIdx.x];
		}		
		__syncthreads();
/*
		if (blockIdx.x == LOGBLOCK && threadIdx.x == LOGTHREAD && LOGTYPE == 1) {
			box->outdata[3 + box->step * 10] = data_ptr[0];	// LJ pot


			Float3 pos_temp = compound.particles[threadIdx.x].pos_tsub1;
			applyHyperpos(&compound_state.positions[threadIdx.x], &pos_temp);
			box->outdata[4 + box->step * 10] = (compound_state.positions[threadIdx.x] - pos_temp).len() / box->dt;

			box->outdata[5 + box->step * 10] = data_ptr[2];// closest particle
			box->outdata[6 + box->step * 10] = data_ptr[1];// force.len();
		}
*/
		int step_offset = (box->step % STEPS_PER_TRAINDATATRANSFER) * MAX_COMPOUND_PARTICLES * 6;
		box->data_GAN[0 + threadIdx.x * 6 + step_offset] = compound_state.positions[threadIdx.x];
//		box->data_GAN[1 + threadIdx.x * 6 + step_offset] = force_bond + force_angle;
//		box->data_GAN[2 + threadIdx.x * 6 + step_offset] = force_LJ_com;
//		box->data_GAN[3 + threadIdx.x * 6 + step_offset] = force_LJ_sol;
	}
	
	// ----------------------------------------------------------------------------- //

	if (force.len() > 200e+6) {
		printf("Critical force %f\n\n\n", force.len());
		box->critical_error_encountered = true;
	}
		
	
	box->compound_state_array_next[blockIdx.x].loadData(&compound_state);
}












__global__ void solventForceKernel(Box* box) {
#define solvent_index (blockIdx.x * blockDim.x + threadIdx.x)
#define thread_active (solvent_index < box->n_solvents)
	//__shared__ Float3 solvent_positions[THREADS_PER_SOLVENTBLOCK];
	__shared__ Float3 utility_buffer[MAX_COMPOUND_PARTICLES];
	__shared__ uint8_t utility_buffer_small[MAX_COMPOUND_PARTICLES];

	float potE_sum = 0;
	float data_ptr[4];	// Pot, force, closest particle, ?
	for (int i = 0; i < 4; i++)
		data_ptr[i] = 0;
	data_ptr[2] = 9999.f;
	Float3 force(0,0,0);

	//printf("mass %f\n", forcefield_device.particle_parameters[ATOMTYPE_SOL].mass);


	Solvent solvent;
	Float3 solvent_pos;
	if (thread_active) {
		solvent = box->solvents[solvent_index];
		solvent_pos = solvent.pos;									// Use this when applying hyperpos, NOT SOLVENT.POS as others acess this concurrently on other kernels!
	}
	
	data_ptr[0] = solvent.pos_tsub1.x;
	data_ptr[1] = solvent.pos_tsub1.y;
	data_ptr[2] = solvent.pos_tsub1.z;


	
	// --------------------------------------------------------------- Molecule Interactions --------------------------------------------------------------- //
	for (int i = 0; i < box->n_compounds; i++) {
		int n_compound_particles = box->compound_state_array[i].n_particles;

		// First all threads help loading the molecule
		if (threadIdx.x < n_compound_particles) {
			utility_buffer[threadIdx.x] = box->compound_state_array[i].positions[threadIdx.x];
			utility_buffer_small[threadIdx.x] = box->compounds[i].atom_types[threadIdx.x];
		}
		__syncthreads();

		//Only do if mol is neighbor to the particular solvent, to save time?
		//LIMAENG::applyHyperpos(&solvent_positions[threadIdx.x], &utility_buffer[0]);									// Move own particle in relation to compound-key-position
		//force += computeLJForces(&solvent_positions[threadIdx.x], n_particles, utility_buffer, data_ptr, &potE_sum);

		if (thread_active) {
			LIMAENG::applyHyperpos(&utility_buffer[0], &solvent_pos);									// Move own particle in relation to compound-key-position
			//force += computeLJForces(&solvent_pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum);
			force += computeLJForces(&solvent_pos, n_compound_particles, utility_buffer, data_ptr, &potE_sum, ATOMTYPE_SOL, utility_buffer_small);
		}
		__syncthreads();


	}
	// ----------------------------------------------------------------------------------------------------------------------------------------------------- //

	if (thread_active) {
		//force += computeSolventToSolventLJForces(&solvent.pos, threadIdx.x, box->n_solvents, solvent_positions, data_ptr, &potE_sum);
		force += computeSolventToSolventLJForces(&solvent_pos, &box->solvent_neighborlists[solvent_index], box->solvents, data_ptr, &potE_sum);
	}





	if (thread_active) {
		int p_index = MAX_COMPOUND_PARTICLES + solvent_index;
		//integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, SOLVENT_MASS, box->dt, &box->thermostat_scalar, p_index);
		integratePosition(&solvent.pos, &solvent.pos_tsub1, &force, forcefield_device.particle_parameters[ATOMTYPE_SOL].mass, box->dt, &box->thermostat_scalar, p_index);
		applyPBC(&solvent.pos);	// forcePositionToInsideBox	// This break the integration, as that doesn't accound for PBC in the vel conservation

	}


	// ------------------------------------ DATA LOG ------------------------------- //
	if (thread_active) {
		int compounds_offset = box->n_compounds * MAX_COMPOUND_PARTICLES;
		//int solvent_index = threadIdx.x;		// Temporary
		int step_offset = (box->step % STEPS_PER_LOGTRANSFER) * box->total_particles_upperbound;

		box->potE_buffer[compounds_offset + solvent_index + step_offset] = potE_sum;	//  data_ptr[0];
		//printf("pot: %f\n", box->potE_buffer[compounds_offset + solvent_index + (box->step) * box->total_particles]);
		box->traj_buffer[compounds_offset + solvent_index + step_offset] = solvent.pos;
	}

	// ----------------------------------------------------------------------------- //


	if (thread_active) {
		box->solvents_next[solvent_index] = solvent;
	}


#undef solvent_index
#undef thread_active
}






__global__ void compoundBridgeKernel(Box* box) {
//#define compound_id blockIdx.x
#define particle_id_bridge threadIdx.x
//#define particle_active particle_id_bridge < bridge.n_particles
	__shared__ CompoundBridgeCompact bridge;
	__shared__ Float3 positions[MAX_PARTICLES_IN_BRIDGE];
	__shared__ Float3 utility_buffer[MAX_PARTICLES_IN_BRIDGE];							// waaaaay too biggg
	__shared__ uint8_t utility_buffer_small[MAX_PARTICLES_IN_BRIDGE];

	if (threadIdx.x == 0) {
		bridge.loadMeta(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	}
	__syncthreads();
	bridge.loadData(&box->bridge_bundle->compound_bridges[blockIdx.x]);
	positions[particle_id_bridge] = Float3(0.f);
	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		positions[particle_id_bridge] = box->compound_state_array[p_ref->compound_id].positions[p_ref->local_id];
	}
	__syncthreads();



	float potE_sum = 0;

	Float3 force(0.f);

	// ------------------------------------------------------------ Intramolecular Operations ------------------------------------------------------------ //
	{											// So for the very first step, these ´should all be 0, but they are not??										TODO: Look into this at some point!!!!

		LIMAENG::applyHyperpos(&positions[0], &positions[particle_id_bridge]);
		__syncthreads();	
		force += computePairbondForces(&bridge, positions, utility_buffer, &potE_sum);
		force += computeAnglebondForces(&bridge, positions, utility_buffer, &potE_sum);
	}
	__syncthreads();


	if (particle_id_bridge < bridge.n_particles) {
		ParticleRefCompact* p_ref = &bridge.particle_refs[particle_id_bridge];
		box->compounds[p_ref->compound_id].forces[p_ref->local_id] = force;
		//force.print(char((int)threadIdx.x+48));
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