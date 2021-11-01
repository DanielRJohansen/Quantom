#include "Engine.cuh"




Simulation* Engine::prepSimulation(Simulation* simulation) {
	this->simulation = simulation;
	srand(290128301);
	boxbuilder.build(simulation);
	printf("Boxbuild complete!\n");

	updateNeighborLists();
	printf("Neighborlists ready\n");



	simulation->moveToDevice();
	return this->simulation;
}





void Engine::updateNeighborLists() {	// Write actual function later;
	int maxc = 1'000'000; // this is temporary!
	CompoundState* statebuffer_host = new CompoundState[maxc];
	CompoundNeighborInfo* neighborlists_host = new CompoundNeighborInfo[maxc];
	cudaMemcpy(statebuffer_host, simulation->box->compound_state_buffer, sizeof(CompoundState) * maxc, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();


	// This only needs to be done the first time... Or does it????
	for (int i = 0; i < maxc; i++) {
		neighborlists_host[i].n_neighbors = 0;
	}
		
	printf("1\n");

	// This is the temp func //
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		for (int j = 0; j < simulation->box->n_compounds; j++) {
			if (i != j) {
				CompoundNeighborInfo* nlist = neighborlists_host;
				nlist->neighborcompound_indexes[nlist->n_neighbors++] = j;
			}
		}
	}
	// --------------------- //

	cudaMemcpy(simulation->box->compound_neighborinfo_buffer, neighborlists_host, sizeof(CompoundNeighborInfo) * maxc, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
}






//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//


void Engine::step() {
	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Error before step!");
		exit(1);
	}


	auto t0 = std::chrono::high_resolution_clock::now();


	int compounds_per_sm = 1000;
	Box* box = simulation->box;	// Why the fuck do i have to do this, VisualStudio???!
	for (int i = 0; i < N_STREAMS; i++) {
		int offset = i * compounds_per_sm;
		//intramolforceKernel <<< compounds_per_sm, 3, 0, stream[i] >>> (box, offset);
	}
	forceKernel <<< box->n_compounds, 256 >>> (box);
	cudaDeviceSynchronize();


	auto t1 = std::chrono::high_resolution_clock::now();



	int blocks_handled = 0;
	while (blocks_handled < sim_blocks) {
		for (int i = 0; i < N_STREAMS; i++) {
			//stepKernel << < BLOCKS_PER_SM, MAX_FOCUS_BODIES, 0, stream[i] >> > (simulation, blocks_handled);
			blocks_handled += BLOCKS_PER_SM;
			if (blocks_handled >= sim_blocks)
				break;
		}

		cudaDeviceSynchronize();
		if (cudaGetLastError() != cudaSuccess) {
			fprintf(stderr, "Error during step :/\n");
			exit(1);
		}
	}




	auto t2 = std::chrono::high_resolution_clock::now();

	blocks_handled = 0;
	while (blocks_handled < sim_blocks) {
		for (int i = 0; i < N_STREAMS; i++) {
			//updateKernel << < BLOCKS_PER_SM, dim3(3,3,3), 0, stream[i] >> > (simulation, blocks_handled);
			blocks_handled += BLOCKS_PER_SM;
			if (blocks_handled >= sim_blocks)
				break;
		}

		cudaDeviceSynchronize();
		if (cudaGetLastError() != cudaSuccess) {
			fprintf(stderr, "Error during update :/\n");
			exit(1);
		}
	}

	auto t3 = std::chrono::high_resolution_clock::now();

	bool verbose = true;
	if (verbose) {
		int intra_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
		int inter_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		int update_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t3 - t2).count();
		timings = timings + Int3(intra_duration, inter_duration, update_duration);
		//printf("\nStep %d ys.\tUpdate: %d\n\n", step_duration, update_duration);
	}


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



																				// TODO: MAKE IS SUCH THAT A BODY CAN NEVER BE EXACTLY ON THE EDGE OF FOCUS, THUS APPEARING IN MULTIPLE GROUPS!



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


/*
constexpr float sigma = 0.3923;	//nm
constexpr float epsilon = 0.5986 * 1'000; //kJ/mol | J/mol
__device__ Float3 calcLJForce(Particle* particle0, Particle* particle1, int i, int bid) {	// Applying force to p0 only! 
	float dist = (particle0->pos - particle1->pos).len();
	float fraction = sigma / dist;

	float f2 = fraction * fraction;
	float f6 = f2 * f2 * f2;
	float f12 = f6 * f6;

	float LJ_pot = 4 * epsilon * (f12 - f6);
	Float3 force_unit_vector = (particle0->pos - particle1->pos).norm();	// + is repulsive, - is attractive


	if (LJ_pot > 20'000) {
		//body0->molecule_type = 99;
		printf("\n\n KILOFORCE! Block %d thread %d\n", bid, threadIdx.x);
		printf("other body index: %d\n", i);
		printf("Body 0 id: %d     %f %f %f\n", particle0->id, particle0->pos.x, particle0->pos.y, particle0->pos.z);
		printf("Body 1 id: %d     %f %f %f\n", particle1->id, particle1->pos.x, particle1->pos.y, particle1->pos.z);
		printf("Distance: %f\n", (particle0->pos - particle1->pos).len());

	}
	return force_unit_vector * LJ_pot;
}
*/


/*
constexpr float kb = 17.5 * 10e+6;		//	J/(mol*nm^2)
__device__ void calcPairbondForce(Compound_H2O* compound, PairBond* pairbond, float* dataptr) {
	Float3 particle1_mirrorpos = getClosestMirrorPos(compound->particles[pairbond->atom_indexes[1]].pos, compound->particles[pairbond->atom_indexes[0]].pos);
	Float3 direction = compound->particles[pairbond->atom_indexes[0]].pos - particle1_mirrorpos;
	//int focus_index = pairbond->atom_indexes[1] == threadIdx.x;	// If it is not pos 1 we get false, meaning pos 0... Not beautiful but it works.
	//direction = direction - direction * 2 * (1 * focus_index);	// Flip dir if focusatom is at index 1

	if (pairbond->atom_indexes[1] == threadIdx.x)	// Flip so we repel the particle in focus
		direction = direction * -1;
	

	float dist = direction.len();
	float dif = dist - pairbond->reference_dist;
	if (dif > 0)	// Flip to attraction
		direction = direction * -1;

	float force_scalar = 0.5 * kb * (dif * dif);
	direction = direction.norm();
	//float invert_if_attraction = -1 + (2 * (dif < 0));

	compound->particles[threadIdx.x].force = compound->particles[threadIdx.x].force + direction * force_scalar;

	if (compound->startindex_particle + threadIdx.x == LOG_P_ID) {
		*dataptr = dist;
	}
		
	Float3 p0p = compound->particles[pairbond->atom_indexes[0]].pos;
	if (force_scalar > WARN_FORCE) {
		printf("\n\n Atom id %d dist %f dif: %f FORCE %f Repulsive %d\n", compound->startindex_particle + threadIdx.x, dist, dif, force_scalar, dif < 0);
		(direction * force_scalar).print('b');
		compound->particles[threadIdx.x].force.print('B');
		//printf("p0 %f %f %f \tp1_mirror %f %f %f\n", p0p.x, p0p.y, p0p.z, particle1_mirrorpos.x, particle1_mirrorpos.y, particle1_mirrorpos.z);
	}
}


constexpr float ktheta = 65 * 10e+3;	// J/mol
__device__ void calcAngleForce(Compound_H2O* compound, AngleBond* anglebond, float* dataptr) {	// We fix the middle particle and move the other particles so they are closest as possible
	Float3 p0_mirrorpos = getClosestMirrorPos(compound->particles[anglebond->atom_indexes[0]].pos, compound->particles[anglebond->atom_indexes[1]].pos);
	Float3 p2_mirrorpos = getClosestMirrorPos(compound->particles[anglebond->atom_indexes[2]].pos, compound->particles[anglebond->atom_indexes[1]].pos);

	Float3 v1 = p0_mirrorpos - compound->particles[anglebond->atom_indexes[1]].pos;
	Float3 v2 = p2_mirrorpos - compound->particles[anglebond->atom_indexes[1]].pos;

	Float3 force_direction = p0_mirrorpos-p2_mirrorpos;
	force_direction = force_direction - force_direction * 2 * (threadIdx.x == 2);
	force_direction = force_direction.norm();

	float angle = getAngle(v1, v2);
	//printf("\nangle %f\n", angle);
	float dif = angle - anglebond->reference_theta;
	float force_scalar = 0.5 * ktheta * (dif * dif);




	float invert_if_attraction = -1 + (2 * (dif < 0));

	compound->particles[threadIdx.x].force = compound->particles[threadIdx.x].force + force_direction * force_scalar * invert_if_attraction;


	if (force_scalar > WARN_FORCE || abs(angle) < 0.1) {
		printf("\n####################\n####################\n####################");
		printf("\nParticle ID %d Angle %f Force %f\n", compound->startindex_particle + threadIdx.x, angle, force_scalar);
		(force_direction * force_scalar * invert_if_attraction).print('a');
	}
	

	if (compound->startindex_particle + threadIdx.x == LOG_P_ID) {
		*dataptr = angle;
	}		
}
*/
__device__ void integrateTimestep(CompactParticle* particle, Float3 force, float dt) {	// Kinetic formula: v = sqrt(2*K/m), m in kg
	Float3 vel_next = particle->vel_prev + (force * (1000.f/particle->mass) * dt);
	particle->pos = particle->pos + vel_next * dt;
	particle->vel_prev = vel_next;
}

// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//

__global__ void forceKernel(Box* box) {
	__shared__ Compound_H2O compound;
	//__shared__ Compound_H2O_

	if (threadIdx.x == 0) {
		compound = box->compounds[blockIdx.x];
	}
	__syncthreads();
	bool thread_compound_active = (compound.n_particles > threadIdx.x);




	if (thread_compound_active) {
		if (threadIdx.x == 0) {
			//compound.particles[threadIdx.x].pos.print();
			box->compound_state_buffer[blockIdx.x].particle_cnt = compound.n_particles;
		}
			
		box->compound_state_buffer[blockIdx.x].positions[threadIdx.x] = compound.particles[threadIdx.x].pos;
	}
		
	//box->compounds[blockIdx.x].particles[threadIdx.x].pos = compound.particles[threadIdx.x].pos;
}



































/*
__global__ void intramolforceKernel(Box* box, int offset) {	// 1 thread per particle in compound
	__shared__ Compound_H2O compound;

	uint32_t compound_index = blockIdx.x + offset;
	if (compound_index >= box->n_compounds)
		return;



	if (threadIdx.x == 0) {
		compound = box->compounds[compound_index];
	}
	__syncthreads();
	compound.particles[threadIdx.x].pos = box->particles[compound.startindex_particle + threadIdx.x].pos;
	compound.particles[threadIdx.x].force = Float3(0, 0, 0);
	__syncthreads();

	for (int i = 0; i < compound.n_pairbonds; i++) {	// Bond forces
		PairBond* bond = &compound.pairbonds[i];
		if (bond->atom_indexes[0] == threadIdx.x || bond->atom_indexes[1] == threadIdx.x) {
			calcPairbondForce(&compound, bond, &box->outdata1[box->data1_cnt]);
			if (compound.startindex_particle + threadIdx.x == LOG_P_ID)
				box->data1_cnt++;
		}
	}

	for (int i = 0; i < compound.n_anglebonds; i++) {	// Angle forces
		AngleBond* bond = &compound.anglebonds[i];
		if (bond->atom_indexes[0] == threadIdx.x || bond->atom_indexes[2] == threadIdx.x) {
			//calcAngleForce(&compound, bond, &box->outdata2[box->data2_cnt]);
			if (compound.startindex_particle + threadIdx.x == LOG_P_ID)
				box->data2_cnt++;
		}
	}
	//CompactParticle* particle = &compound.particles[threadIdx.x];
	//PairBond* pairbond = &compound.pairbonds[threadIdx.x];
	//calcPairbondForce(&compound, pairbond);	// This applies the force directly to the particles

	box->particles[compound.startindex_particle + threadIdx.x].force = compound.particles[threadIdx.x].force;
}






__global__ void stepKernel(Simulation* simulation, int offset) {
	int blockID = blockIdx.x + offset;	// Maybe save the register, and just calculate it a couple times.
	int bodyID = threadIdx.x;

	if (blockID >= simulation->box->n_blocks)
		return;
	


	// Load bodies into shared memory
	__shared__ Block block;	
	__shared__ AccessPoint accesspoint;
	
	if (threadIdx.x == 0) {
		block = simulation->box->blocks[blockID];
		accesspoint = AccessPoint();
	}
		
	__syncthreads();
	Particle particle = block.focus_particles[bodyID];

	

	// --------------------------------- ACTUAL MOLECULAR DYNAMICS HAPPENING HERE! --------------------------------- //
	// Calc all Lennard-Jones forces from focus bodies
	if (particle.active) {						// This part acounts for about 2/5 of compute time
		// I assume that all present molecules come in order!!!!!!



		particle.force = simulation->box->particles[particle.id].force;
		simulation->box->blocks[blockID].focus_particles[bodyID].color[2] = simulation->box->particles[particle.id].color[2];
		Float3 force_total;
		for (int i = 0; i < MAX_FOCUS_BODIES; i++) {
			if (block.focus_particles[i].active) {
				if (i != bodyID && particle.compoundID !=  block.focus_particles[i].compoundID) {
					force_total = force_total + calcLJForce(&particle, &block.focus_particles[i], -i, blockID);
				}
			}
			else {
				break;
			}
		}

		// Calc all forces from Near bodies
		for (int i = 0; i < MAX_NEAR_BODIES; i++) {
			if (block.near_particles[i].active && particle.compoundID != block.near_particles[i].compoundID) {
				force_total = force_total + calcLJForce(&particle, &block.near_particles[i], i + MAX_FOCUS_BODIES, blockID);
			}
			else {
				break;
			}
		}
		


		//particle.force = simulation->box->particles[particle.id].force;
		//printf("\n ID %d\tInter: %f %f %f\tIntra %f %f %f\n", particle.id, force_total.x, force_total.y, force_total.z, particle.force.x, particle.force.y, particle.force.z);
		

		//force_total = force_total + particle.force;
		particle.force = particle.force + force_total;
		
		if (particle.id == LOG_P_ID) {
			//simulation->box->outdata3[simulation->box->data3_cnt++] = force_total.len();
			simulation->box->outdata3[simulation->box->data3_cnt++] = force_total.len();
			simulation->box->outdata4[simulation->box->data4_cnt++] = particle.force.len();
		}

		if (particle.force.len() > WARN_FORCE) {
			printf("\n");
			particle.force.print('I');
			force_total.print('T');
		}
		if (particle.force.len() > END_SIM_FORCE) {
			simulation->finished = true;
			printf("Ending due to particle %d\n", particle.id);
		}


		// Integrate position  AFTER ALL BODIES IN BLOCK HAVE BEEN CALCULATED? No should not be a problem as all update their local body, 
		// before moving to shared?? Although make the local just a pointer might be faster, since a SimBody might not fit in thread registers!!!!!
		integrateTimestep(simulation, &particle);




		// Correct for bonded-body constraints? Or maybe this should be calculated as a force too??




		
		// Swap with mirror image if body has moved out of box
		Float3 hyper_pos = getHyperPosition(particle.pos);
		if ((hyper_pos - particle.pos).len() > 1) {	// If the hyperposition is different, the body is out, and we import the mirror body
			//printf("\nSwapping Body %f %f %f to hyperpos %f %f %f\n\n", body.pos.x, body.pos.y, body.pos.z, hyper_pos.x, hyper_pos.y, hyper_pos.z);
			particle.pos = hyper_pos;
		}	
		simulation->box->particles[particle.id].pos = particle.pos;

	}






	// Publish new positions for the focus bodies
	accesspoint.particles[bodyID] = particle;

	// Mark all bodies as obsolete
	simulation->box->blocks[blockID].focus_particles[bodyID].active = false;	// Need to run update kernel before rendering, or no particle will be rendered.


	__syncthreads();
	if (bodyID == 0) {
		simulation->box->accesspoint[blockID] = accesspoint;
	}
} 












__global__ void updateKernel(Simulation* simulation, int offset) {
	int blockID1 = blockIdx.x + offset;
	if (blockID1 >= simulation->box->n_blocks)
		return;

	int threadID1 = indexConversion(Int3(threadIdx.x, threadIdx.y, threadIdx.z), 3);

	

	
	__shared__ Float3 block_center;
	__shared__ Int3 blockID3;
	__shared__ int bpd;
	__shared__ char element_cnt_focus[27];
	__shared__ char element_sum_focus[27 + 1];
	__shared__ short int element_cnt_near[27];
	__shared__ short int element_sum_near[27+1];
	__shared__ char relation_array[27][MAX_FOCUS_BODIES];	// 0 = far, 1 = near, 2 = focus
	

	

	//element_sum_focus[threadID1 + 1] = -1;
	//element_sum_near[threadID1 + 1] = -1;
	for (int i = 0; i < MAX_FOCUS_BODIES; i++)
		relation_array[threadID1][i] = 0;
	if (threadID1 == 0) {
		block_center = simulation->box->blocks[blockID1].center;			// probably faster to just calculate
		bpd = simulation->blocks_per_dim;
		blockID3 = indexConversion(blockID1, bpd);
		element_sum_focus[0] = 0;
		element_sum_near[0] = 0;	
	}
	
	__syncthreads();


	AccessPoint accesspoint;		// Personal to all threads
		
		Int3 neighbor_index3 = blockID3 + (Int3(threadIdx.x, threadIdx.y, threadIdx.z) - Int3(1, 1, 1));
		neighbor_index3 = neighbor_index3 + Int3(bpd * (neighbor_index3.x == -1), bpd * (neighbor_index3.y == -1), bpd * (neighbor_index3.z == -1));
		neighbor_index3 = neighbor_index3 - Int3(bpd * (neighbor_index3.x == bpd), bpd * (neighbor_index3.y == bpd), bpd * (neighbor_index3.z == bpd));
		accesspoint = simulation->box->accesspoint[indexConversion(neighbor_index3, bpd)];
	
	

	{	// To make these to variables temporary
		int focus_cnt = 0;
		int near_cnt = 0;
		int array_index = threadID1 + 1;
		for (int i = 0; i < MAX_FOCUS_BODIES; i++) {
			Particle* particle = &accesspoint.particles[i];
			if (!particle->active)
				break;

			Float3 closest_mirror_pos = getClosestMirrorPos(particle->pos, block_center);
			int relation_type = (bodyInNear(&closest_mirror_pos, &block_center) + bodyInFocus(&particle->pos, &block_center) * 2);
			



			if (relation_type == 1)		// If the body is ONLY in near, we make a temporary copy, with the mirrors position.
				particle->pos = closest_mirror_pos;

			//if (relation_type > 1 && body->id == 0)
				//printf("\nLoading %d to block %d %d %d by thread %d %d %d from block %d %d %d\n", body->id, blockID3.x, blockID3.y, blockID3.z, threadIdx.x-1, threadIdx.y - 1, threadIdx.z - 1, neighbor_index3.x, neighbor_index3.y, neighbor_index3.z);

			relation_array[threadID1][i] = relation_type;
			near_cnt += (relation_type == 1);
			focus_cnt += (relation_type > 1);
		}

		for (int i = 0; i < 27; i++) {	// Reserve spaces 
			if (threadID1 == i) {
				element_sum_focus[array_index] = element_sum_focus[array_index-1] + focus_cnt;
				element_sum_near[array_index] = element_sum_near[array_index-1] + near_cnt;
				if (element_sum_focus[array_index] >= MAX_FOCUS_BODIES || element_sum_near[array_index] >= MAX_NEAR_BODIES) {
					printf("Block %d overflowing! near %d    focus %d\n", blockID1, element_sum_near[i + 1], element_sum_focus[i + 1]);
					simulation->finished = true;
					break;
				}
			}
			__syncthreads();
		}
		
	}

	
	{
		int focus_index = element_sum_focus[threadID1];
		int near_index = element_sum_near[threadID1];
		for (int i = 0; i < MAX_FOCUS_BODIES; i++) {
			if (relation_array[threadID1][i] > 1) 
				simulation->box->blocks[blockID1].focus_particles[focus_index++] = accesspoint.particles[i];
			else if (relation_array[threadID1][i] == 1) {
				simulation->box->blocks[blockID1].near_particles[near_index++] = accesspoint.particles[i];
				//printf("\n%d\t %d %d %d\t%d\t nearindex: %d\n", blockID1, blockID3.x, blockID3.y, blockID3.z, accesspoint.bodies[i].molecule_type, near_index);
			}	
		}



		// Handle deactivating now-obsolete bodies
		// The stepkernel will handle the marking of focus bodies, as it has the correct amount of threads. Less overhead!

		// For nearbodies we only need to terminate 1 body, this saves alot of writes to global!
		if (threadID1 == 26) {	
			if (near_index < (MAX_NEAR_BODIES))
				simulation->box->blocks[blockID1].near_particles[near_index].active = false;
			//printf("\n%f %f %f\tExporting %d\n", block_center.x, block_center.y, block_center.z, near_index);
			//printf("Block %d loaded %d bodies\n", blockID1, focus_index);
		}
	}
}

*/

