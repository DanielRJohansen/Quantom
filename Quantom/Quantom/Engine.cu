#include "Engine.cuh"






Simulation* Engine::prepSimulation(Simulation* simulation) {
	this->simulation = simulation;
	simulation->bodies = new SimBody[simulation->n_bodies];


	int n_blocks = initBlocks();
	//linkBlocks();
	//prepareEdgeBlocks();

	srand(290128301);
	int n_bodies = fillBox();

	int update_bytesize = sizeof(Block) + sizeof(Int3) + sizeof(int) + sizeof(char) * 27 * 3 + sizeof(char) + sizeof(SimBody) * 27 + sizeof(unsigned short int) * 27 * 2 + sizeof(char) * 27 *MAX_FOCUS_BODIES;

	printf("\nSimbody size: %d bytes\n", sizeof(SimBody));
	printf("Block size: %d\n", sizeof(Block));
	printf("Simulation configured with %d blocks, and %d bodies. Approximately %02.2f bodies per block. \n", n_blocks, n_bodies, ((float)n_bodies/n_blocks));
	printf("Required shared mem per block for stepKernel: %d KB\n", (int) ((sizeof(Block)) / 1000.f));
	printf("Required shared mem per block for updateKernel: %d KB\n", (int) (update_bytesize/1000.f));
	printf("Required global mem for Box: %d MB\n", (int) (sizeof(Block) * n_blocks / 1000000.f));


	prepareCudaScheduler();

	return simToDevice();
}


int Engine::countBodies() {
	int count = 0;
	for (int i = 0; i < simulation->box->n_blocks; i++) {
		for (int j = 0; j < MAX_FOCUS_BODIES; j++) {
			if (simulation->box->blocks[i].focus_bodies[j].molecule_type == UNUSED_BODY)
				break;
			count++;
		}
	}
	printf("\nBody count: %d\n\n", count);
	return count;
}

Simulation* Engine::simToDevice() {
	simulation->moveToDevice();	// Must be done before initiating raytracer!

	Simulation* temp;
	int bytesize = sizeof(Simulation);
	cudaMallocManaged(&temp, bytesize);
	cudaMemcpy(temp, simulation, bytesize, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	delete simulation;
	simulation = temp;

	return simulation;
}

int Engine::initBlocks() {
	

	block_dist = FOCUS_LEN;
	bpd = ((simulation->box_size)/ block_dist) + 1;		// Otherwise no block focus on edge particles
	int n_blocks = pow(bpd, 3);
	box_base = - (BOX_LEN/2.f);			// Is the center of first block!
	effective_box_base = box_base - FOCUS_LEN_HALF;


	printf("Blocks per dim: %d\n", bpd);

	simulation->blocks_per_dim = bpd;
	simulation->box->blocks_per_dim = bpd;
	simulation->box->n_blocks = n_blocks;
	simulation->box->blocks = new Block[n_blocks];
	//simulation->box->accesspoint = new AccessPoint[n_blocks];


	int index = 0;
	float offset = -simulation->box_size / 2 + 0.5 * BLOCK_LEN;
	for (int x = 0; x < bpd; x++) {
		for (int y = 0; y < bpd; y++) {
			for (int z = 0; z < bpd; z++) {
				Float3 center(x * block_dist + box_base, y * block_dist + box_base, z * block_dist + box_base);
				//center.print();
				simulation->box->blocks[index] = Block(center);
				index++;
			}
		}
	}
	return index;
}

int Engine::fillBox() {

	int bodies_per_dim = ceil(cbrt((float)simulation->n_bodies));
	float body_edge_mindist = 0.5;

	float dist = (BOX_LEN - body_edge_mindist*2) / (float)bodies_per_dim;	// dist_per_index
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist);

	float base = -BOX_LEN / 2.f + body_edge_mindist + dist/2.f + 0.000001;
	


	double m = 18.01528;		// g/mol
	double k_B = 8.617333262145 * 10e-5;
	double T = 293;	// Kelvin
	float mean_velocity = m / (2 * k_B * T);
	int p = 10000;

	int index = 0;
	for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
				if (index == simulation->n_bodies)
					break;

				float r1 = rand() % p / (float)p - 0.5;
				float r2 = rand() % p / (float)p - 0.5;
				float r3 = rand() % p / (float)p - 0.5;

				simulation->bodies[index].pos = Float3(base + dist * (float)x_index, base + dist * (float)y_index, base + dist * (float)z_index);

				//printf("Body pos: ");
				//simulation->bodies[index].pos.print();

				simulation->bodies[index].id = index;	// 

				simulation->bodies[index].molecule_type = 0;
				simulation->bodies[index].mass = m;

				simulation->bodies[index].vel_prev = Float3(r1, r2, r3).norm() * mean_velocity;
				simulation->bodies[index].rotation = Float3(0, 0, 0);
				simulation->bodies[index].rot_vel = Float3(r1*r2, 4*PI, 0);
				placeBody(&simulation->bodies[index++]);
			}
		}
	}
	return index;
}
/*
void Engine::placeBody(SimBody* body) {
	//const Int3 block_index = posToBlockIndex(&body->pos);
	

	int count = 0;

	SimBody temp;
	Int3 block_index_;

	for (int z_off = -1; z_off <= 1; z_off++) {
		for (int y_off = -1; y_off <= 1; y_off++) {
			for (int x_off = -1; x_off <= 1; x_off++) {


				Float3 pos_(body->pos.x + x_off * BLOCK_OVERLAP, body->pos.y + y_off * BLOCK_OVERLAP, body->pos.z + z_off * BLOCK_OVERLAP);
				block_index_ = posToBlockIndex(&pos_);
				//printf("%d %d %d\n", block_index_.x, block_index_.y, block_index_.z);
				int block_index_1d = block3dIndexTo1dIndex(block_index_);
				//printf("Block index. %d\n", block_index_1d);
				Block* block = &simulation->box->blocks[block_index_1d];
				if (block->addBody(body))
					count++;
			}
		}
	}
	//printf("Molecule placed in %d blocks\n", count);
}
*/

void Engine::placeBody(SimBody* body) {
	Int3 block_index = posToBlockIndex(&body->pos);
	//printf("Block index: %d %d %d\n", block_index.x, block_index.y, block_index.z);

	int block_index_1d = block3dIndexTo1dIndex(block_index);
	//printf("Block index: %d\n", block_index_1d);

	if (block_index_1d < 0)
		printf("Rebuild All you twat\n");

	Block* block = &simulation->box->blocks[block_index_1d];

	
	if (!block->addBody(body))
		printf("Body lost!\n");

}


void Engine::prepareCudaScheduler() {
	sim_blocks = simulation->box->n_blocks;

	for (int i = 0; i < N_STREAMS; i++)
		cudaStreamCreate(&stream[i]);

	printf("%d kernel launches necessary to step\n", (int) ceil((float)simulation->box->n_blocks / (float)BLOCKS_PER_SM));
	//gridblock_size = dim3(GRIDBLOCKS_PER_BODY, BLOCKS_PER_SM, 1);
}




//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//


void Engine::step() {

	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Error before step!");
		exit(1);
	}


	auto t0 = std::chrono::high_resolution_clock::now();

	int blocks_handled = 0;
	while (blocks_handled < sim_blocks) {
		for (int i = 0; i < N_STREAMS; i++) {
			stepKernel << < BLOCKS_PER_SM, MAX_FOCUS_BODIES, 0, stream[i] >> > (simulation, blocks_handled);
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

	auto t1 = std::chrono::high_resolution_clock::now();

	blocks_handled = 0;
	while (blocks_handled < sim_blocks) {
		for (int i = 0; i < N_STREAMS; i++) {
			updateKernel << < BLOCKS_PER_SM, dim3(3,3,3), 0, stream[i] >> > (simulation, blocks_handled);
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

	auto t2 = std::chrono::high_resolution_clock::now();

	bool verbose = true;
	if (verbose) {
		int step_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
		int update_duration = (int)std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		timings = timings + Int3(step_duration, update_duration, 0);
		//printf("\nStep %d ys.\tUpdate: %d\n\n", step_duration, update_duration);
	}


	simulation->step++;

	
}





// ------------------------------------------------------------------------------------------- DEVICE FUNCTIONS -------------------------------------------------------------------------------------------//




__device__ float cudaMax1(float a, float b) {
	if (a > b)
		return a;
	return b;
}
__device__ float cudaMin1(float a, float b) {
	if (a < b)
		return a;
	return b;
}

__device__ Float3 forceFromDist(Float3 dists) {
	return Float3(
		((edgeforce_scalar / dists.x)- edgeforce_scalar),
		((edgeforce_scalar / dists.y)- edgeforce_scalar),
		((edgeforce_scalar / dists.z)- edgeforce_scalar)
	);
}

__device__ inline Float3 calcEdgeForce(SimBody* body) {
	Float3 dimensionwise_dists_from_base = Float3(body->pos.x - BOX_BASE, body->pos.y - BOX_BASE, body->pos.z - BOX_BASE).abs();
	Float3 dimensionwise_dists_from_end = Float3(body->pos.x - BOX_LEN_HALF, body->pos.y - BOX_LEN_HALF, body->pos.z - BOX_LEN_HALF).abs();

	Float3 pos_forces = forceFromDist(dimensionwise_dists_from_base).zeroIfBelow(0);
	Float3 neg_forces = forceFromDist(dimensionwise_dists_from_end).zeroIfBelow(0);
	/*if (body->id == 1) {
		printf("\n pos: %f %f %f\tForce: %.1f %.1f %.1f\n\n", body->pos.x, body->pos.y, body->pos.z, (pos_forces - neg_forces).x, (pos_forces - neg_forces).y, (pos_forces - neg_forces).z);
	}*/
		
	return (pos_forces - neg_forces);
}
																				// TODO: MAKE IS SUCH THAT A BODY CAN NEVER BE EXACTLY ON THE EDGE OF FOCUS, THUS APPEARING IN MULTIPLE GROUPS!
__device__ bool bodyInNear(SimBody* body, Float3* block_center) {
	Float3 dist_from_center = (body->pos - *block_center).abs();
	return (dist_from_center.x < FOCUS_LEN && dist_from_center.y < FOCUS_LEN && dist_from_center.z < FOCUS_LEN);
}

__device__ bool bodyInFocus(SimBody* body, Float3* block_center) {
	Float3 dist_from_center = (body->pos - *block_center);

	return
		(dist_from_center.x < FOCUS_LEN_HALF&& dist_from_center.y < FOCUS_LEN_HALF&& dist_from_center.z < FOCUS_LEN_HALF)			// Open upper bound
		&&
		((dist_from_center.x >= -FOCUS_LEN_HALF && dist_from_center.y >= -FOCUS_LEN_HALF && dist_from_center.z >= -FOCUS_LEN_HALF));	// Closed lower bound;

	//return (dist_from_center.x < FOCUS_LEN_HALF && dist_from_center.y < FOCUS_LEN_HALF && dist_from_center.z < FOCUS_LEN_HALF);
}

__device__ int indexConversion(Int3 xyz, int elements_per_dim) {
	return int(xyz.x + xyz.y * elements_per_dim + xyz.z * elements_per_dim * elements_per_dim);
}

__device__ Int3 indexConversion(int index, int elements_per_dim) {
	return Int3(
		index % elements_per_dim,
		(index / elements_per_dim) % elements_per_dim,
		index / (elements_per_dim * elements_per_dim)
	);
}



__device__ Float3 calcLJForce(SimBody* body0, SimBody* body1, int i, int bid) {
	float dist_pow2 = (body0->pos - body1->pos).lenSquared();
	float dist_pow2_inv = 1.f / dist_pow2;
	float dist_pow6_inv = dist_pow2_inv * dist_pow2_inv * dist_pow2_inv;

	float LJ_pot = 48 * dist_pow2_inv * dist_pow6_inv * (dist_pow6_inv - 0.5f);
	Float3 force_unit_vector = (body0->pos - body1->pos).norm();
	
	if (LJ_pot > 10e+9) {
		body0->molecule_type = 99;
		printf("\n\n GIGAFORCE! Block %d thread %d\n", bid, threadIdx.x);
		printf("other body index: %d\n", i);
		printf("Body 0 id: %d     %f %f %f\n", body0->id,  body0->pos.x, body0->pos.y, body0->pos.z);
		printf("Body 1 id: %d     %f %f %f\n", body1->id, body1->pos.x, body1->pos.y, body1->pos.z);
		printf("Distance: %f\n", (body0->pos - body1->pos).len());
		//printf("Force: %f\n\n\n", LJ_pot);
		//force_unit_vector.print();
	}
	return force_unit_vector * LJ_pot;
}

__device__ void integrateTimestep(Simulation* simulation, SimBody* body, Float3 force) {
	Float3 vel_next = body->vel_prev + (force * (simulation->dt / body->mass));
	body->pos = body->pos + vel_next * simulation->dt;
	body->vel_prev = vel_next;
}




// ------------------------------------------------------------------------------------------- KERNELS -------------------------------------------------------------------------------------------//








__global__ void stepKernel(Simulation* simulation, int offset) {
	int blockID = blockIdx.x + offset;	// Maybe save the register, and just calculate it a couple times.
	int bodyID = threadIdx.x;

	if (blockID >= simulation->box->n_blocks)
		return;
	


	// Load bodies into shared memory
	__shared__ Block block;	
	__shared__ AccessPoint accesspoint;
	
	__syncthreads();
	if (threadIdx.x == 0) {
		block = simulation->box->blocks[blockID];
		accesspoint = AccessPoint();
	}
		
	__syncthreads();
	SimBody body = block.focus_bodies[bodyID];





	// ACTUAL MOLECULAR DYNAMICS HAPPENING HERE! //

	Float3 force_total = Float3(0, 0, 0);
	// Calc all Lennard-Jones forces from focus bodies
	if (body.molecule_type != UNUSED_BODY) {

		// I assume that all present molecules come in order!!!!!!
		for (int i = 0; i < MAX_FOCUS_BODIES; i++) {
			if (block.focus_bodies[i].molecule_type != UNUSED_BODY) {
				if (i != bodyID) {
					force_total = force_total + calcLJForce(&body, &block.focus_bodies[i], -i, blockID);
				}
			}
			else {
				break;
			}
		}

		// Calc all forces from Near bodies
		for (int i = 0; i < MAX_NEAR_BODIES; i++) {
			if (block.near_bodies[i].molecule_type != UNUSED_BODY) {
				force_total = force_total + calcLJForce(&body, &block.near_bodies[i], i + MAX_FOCUS_BODIES, blockID);
			}
			else {
				break;
			}
		}


	}
	
	

	//body.vel = body.vel + LJ_pot_force * simulation->dt;



	body.rotation = body.rotation + body.rot_vel * simulation->dt;				



	// Calc all forces from box edge
	force_total = force_total + calcEdgeForce(&body);
	//if (body.id == 1)
		//printf("\nTotal Force: %.2f %.2f %.2f\n", force_total.x, force_total.y, force_total.z);

	if (body.molecule_type == 99)	// TEMP
		simulation->finished = true;

	// Integrate position  AFTER ALL BODIES IN BLOCK HAVE BEEN CALCULATED? No should not be a problem as all update their local body, 
	// before moving to shared?? Although make the local just a pointer might be faster, since a SimBody might not fit in thread registers!!!!!
	//body.pos = body.pos + body.vel * simulation->dt;
	integrateTimestep(simulation, &body, force_total);




	// Correct for bonded-body constraints? Or maybe this should be calculated as a force too??




	// Ensure no position is EXACTLY on the block edge, causing it to be lost forever!




	// Publish new positions for the focus bodies
	
	//block.focus_bodies[bodyID] = body;	// Probably useless, right?
	accesspoint.bodies[bodyID] = body;






	// Mark all bodies as obsolete


	__syncthreads();
	if (bodyID == 0) {
		simulation->box->accesspoint[blockID] = accesspoint;
		//simulation->box->blocks[blockID] = block;	// Very expensive, only needed for rendering.	NOT NECESSARY IF WE RUN UPDATEKERNEL AFTER STEP BEFORE RENDER!!
	}
} 









__global__ void updateKernel(Simulation* simulation, int offset) {
	int blockID1 = blockIdx.x + offset;
	if (blockID1 >= simulation->box->n_blocks)
		return;

	int threadID1 = indexConversion(Int3(threadIdx.x, threadIdx.y, threadIdx.z), 3);

	



	
	//__shared__ Block block;
	__shared__ Float3 block_center;
	__shared__ Int3 blockID3;
	__shared__ int bpd;
	__shared__ char element_cnt_focus[27];
	__shared__ char element_sum_focus[27 + 1];
	__shared__ short int element_cnt_near[27];
	__shared__ short int element_sum_near[27+1];
	__shared__ char relation_array[27][MAX_FOCUS_BODIES];	// 0 = far, 1 = near, 2 = focus
	
	AccessPoint accesspoint;

	

	element_sum_focus[threadID1 + 1] = -1;
	element_sum_near[threadID1 + 1] = -1;
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


	{	
		Int3 neighbor_index3 = blockID3 + (Int3(threadIdx.x, threadIdx.y, threadIdx.z) - Int3(1, 1, 1));
		neighbor_index3 = neighbor_index3 + Int3(bpd * (neighbor_index3.x == 0), bpd * (neighbor_index3.y == 0), bpd * (neighbor_index3.z == 0));
		neighbor_index3 = neighbor_index3 - Int3(bpd * (neighbor_index3.x == (bpd-1)), bpd * (neighbor_index3.y == (bpd-1)), bpd * (neighbor_index3.z == (bpd-1)));
		accesspoint = simulation->box->accesspoint[indexConversion(neighbor_index3, bpd)];
	}
	
		


	
	
	{	// To make these to variables temporary
		int focus_cnt = 0;
		int near_cnt = 0;
		int array_index = threadID1 + 1;
		for (int i = 0; i < MAX_FOCUS_BODIES; i++) {
			SimBody* body = &accesspoint.bodies[i];
			if (body->molecule_type == UNUSED_BODY)
				break;

			int relation_type = (bodyInNear(body, &block_center) + bodyInFocus(body, &block_center) * 2);

			relation_array[threadID1][i] = relation_type;
			near_cnt += (relation_type == 1);
			focus_cnt += (relation_type > 1);
		}
		
		/*
		while (true) {
			if (element_sum_near[threadID1] > -1) {
				element_sum_focus[array_index] = element_sum_focus[threadID1] + focus_cnt;
				element_sum_near[array_index] = element_sum_near[threadID1] + near_cnt;
				break;
			}
		}*/
		// Safer solution below, with no while true!
		/**/
		for (int i = 0; i < 27; i++) {
			if (threadID1 == i) {
				element_sum_focus[array_index] = element_sum_focus[i] + focus_cnt;
				element_sum_near[array_index] = element_sum_near[i] + near_cnt;
				if (element_sum_focus[array_index] == MAX_FOCUS_BODIES+1 || element_sum_near[i+1] == MAX_NEAR_BODIES + 1) {
					printf("Block %d overflowing! near %d    focus %d\n", blockID1, element_sum_near[i + 1], element_sum_focus[i + 1]);
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
				simulation->box->blocks[blockID1].focus_bodies[focus_index++] = accesspoint.bodies[i];				
			else if (relation_array[threadID1][i] == 1)
				simulation->box->blocks[blockID1].near_bodies[near_index++] = accesspoint.bodies[i];
		}



		// Handle deactivating now-obsolete bodies
		if (threadID1 >= element_sum_focus[26]) {	// Each thread makes sure its corresponding focusbody is deactivated
			simulation->box->blocks[blockID1].focus_bodies[threadID1].molecule_type = UNUSED_BODY;
		}
		if (threadID1 == 26) {	// For nearbodies we only need to terminate 1 body, this saves alot of writes to global!
			if (near_index < (MAX_NEAR_BODIES))
				simulation->box->blocks[blockID1].near_bodies[near_index].molecule_type = UNUSED_BODY;
		}
	}
	

}




/*if (relation_array[threadID1][i] > 1)
	block.focus_bodies[focus_index++] = accesspoint.bodies[i];
else if (relation_array[threadID1][i] == 1)
	block.near_bodies[near_index++] = accesspoint.bodies[i];
	*/

	// Delete all previous bodies 
		/*
		if (threadID1 < MAX_FOCUS_BODIES)
			block.focus_bodies[threadID1] = SimBody();
		for (int i = threadID1; i <= MAX_NEAR_BODIES + 27; i += 27)
			if (i < MAX_NEAR_BODIES)
				block.near_bodies[i] = SimBody();
		*/

