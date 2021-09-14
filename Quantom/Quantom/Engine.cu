#include "Engine.cuh"




Simulation* Engine::prepSimulation(Simulation* simulation) {
	this->simulation = simulation;
	simulation->bodies = new SimBody[simulation->n_bodies];


	int n_blocks = initBlocks();
	//linkBlocks();
	//prepareEdgeBlocks();

	srand(290128301);
	int n_bodies = fillBox();


	printf("\nSimbody size: %d bytes\n", sizeof(SimBody));
	printf("Block size: %d\n", sizeof(Block));
	printf("Simulation configured with %d blocks, and %d bodies. Approximately %d bodies per block. \n", n_blocks, n_bodies, n_bodies/n_blocks);
	printf("Required shared mem for stepKernel: %d\n", sizeof(Block));
	printf("Required global mem for Box: %d MB\n", (int) (sizeof(Block) * n_blocks / 1000000.f));
	//exit(1);


	prepareCudaScheduler();

	return simToDevice();
}


int Engine::countBodies() {
	int count = 0;
	for (int i = 0; i < simulation->box->n_blocks; i++) {
		for (int j = 0; j < MAX_FOCUS_BODIES; j++) {
			if (simulation->box->blocks[i].focus_bodies[j].molecule_type != UNUSED_BODY)
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
				simulation->box->blocks[index++] = Block(center);
			}
		}
	}
	return index;
}

int Engine::fillBox() {

	int bodies_per_dim = ceil(cbrt((float)simulation->n_bodies));
	float body_edge_mindist = 0.5;

	printf("Bodies per dim: %d\n", bodies_per_dim);
	float dist = (BOX_LEN - body_edge_mindist*2) / (float)bodies_per_dim;	// dist_per_index
	float base = -BOX_LEN / 2.f + body_edge_mindist + dist/2.f;
	float vel_scalar = 4;

	int index = 0;
	for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
				if (index == simulation->n_bodies)
					break;

				int p = 10000;
				float r1 = rand() % p / (float)p - 0.5;
				float r2 = rand() % p / (float)p - 0.5;
				float r3 = rand() % p / (float)p - 0.5;

				simulation->bodies[index].pos = Float3(base + dist * (float)x_index, base + dist * (float)y_index, base + dist * (float)z_index);

				//printf("Body pos: ");
				//simulation->bodies[index].pos.print();

				simulation->bodies[index].molecule_type = 0;
				simulation->bodies[index].vel = Float3(r1 * vel_scalar, r2 * vel_scalar, r3 * vel_scalar);
				simulation->bodies[index].rotation = Float3(0, 0, 0);
				simulation->bodies[index].rot_vel = Float3(0, 4*PI, 0);
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
	auto start = std::chrono::high_resolution_clock::now();

	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Error before step!");
		exit(1);
	}




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

	cudaDeviceSynchronize();
	
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



	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	printf("Step time: %d ys.\tStep #%05d", duration.count(), simulation->step++);
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
		((edgeforce_scalar / dists.x)-0.25),
		((edgeforce_scalar / dists.y)-0.25),
		((edgeforce_scalar / dists.z)-0.25)
	);
}

__device__ inline Float3 calcEdgeForce(SimBody* body) {
	//Float3 dists1(body->pos.x + BOX_LEN / 2.f, body->pos.y + BOX_LEN / 2.f, body->pos.z + BOX_LEN / 2.f);
	//Float3 dists2(BOX_LEN / 2.f - body->pos.x, BOX_LEN / 2.f - body->pos.y, BOX_LEN / 2.f - body->pos.z);

	Float3 dists_from_base = Float3(body->pos.x - BOX_BASE, body->pos.y - BOX_BASE, body->pos.z - BOX_BASE).abs();
	Float3 dists_from_end = Float3(body->pos.x - BOX_LEN_HALF, body->pos.y - BOX_LEN_HALF, body->pos.z - BOX_LEN_HALF).abs();

	Float3 pos_forces = forceFromDist(dists_from_base).zeroIfBelow(0);
	Float3 neg_forces = forceFromDist(dists_from_end).zeroIfBelow(0);
	return pos_forces - neg_forces;
}
																				// TODO: MAKE IS SUCH THAT A BODY CAN NEVER BE EXACTLY ON THE EDGE OF FOCUS, THUS APPEARING IN MULTIPLE GROUPS!
__device__ bool bodyInNear(SimBody* body, Float3* block_center) {
	Float3 dist_from_center = (body->pos - *block_center).abs();
	return (dist_from_center.x < FOCUS_LEN && dist_from_center.y < FOCUS_LEN && dist_from_center.z < FOCUS_LEN);
}

__device__ bool bodyInFocus(SimBody* body, Float3* block_center) {
	Float3 dist_from_center = (body->pos - *block_center).abs();
	return (dist_from_center.x < FOCUS_LEN_HALF && dist_from_center.y < FOCUS_LEN_HALF && dist_from_center.z < FOCUS_LEN_HALF);
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
	
	if (LJ_pot > 1000000) {
		printf("\n\n MEGAFORCE! Block %d thread %d\n", bid, threadIdx.x);
		/*printf("other body id: %d\n", i);
		printf("Body 0 %f %f %f\n", body0->pos.x, body0->pos.y, body0->pos.z);
		printf("Body 1 %f %f %f\n", body1->pos.x, body1->pos.y, body1->pos.z);*/
		printf("Distance: %f\n", (body0->pos - body1->pos).len());
		//printf("Force: %f\n\n\n", LJ_pot);
		//force_unit_vector.print();
	}
	return force_unit_vector * LJ_pot;
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
	}
		
	__syncthreads();
	SimBody body = block.focus_bodies[bodyID];


	// Always need thread0 to send block global
	// Also compute for non active blocks, to avoid branching!



	Float3 LJ_pot_force = Float3(0, 0, 0);
	// Calc all forces from focus bodies
	for (int i = 0; i < MAX_FOCUS_BODIES; i++) {
		if (i != bodyID && block.focus_bodies[i].molecule_type != UNUSED_BODY) {
			LJ_pot_force = LJ_pot_force + calcLJForce(&body, &block.focus_bodies[i], i, blockID);
		}
	}
	for (int i = 0; i < MAX_NEAR_BODIES; i++) {
		if (i != bodyID && block.near_bodies[i].molecule_type != UNUSED_BODY) {
			LJ_pot_force = LJ_pot_force + calcLJForce(&body, &block.near_bodies[i], i, blockID);
		}
	}
	/*int i = 0;
	while (i < MAX_FOCUS_BODIES) {
		i = i + (i == bodyID) 
	}*/
	// Calc all forces from Near bodies
	

	body.vel = body.vel + LJ_pot_force * simulation->dt;



	//body.rotation = body.rotation + body.rot_vel * simulation->dt;				


	if (body.molecule_type != UNUSED_BODY) {
		//LJ_pot_force.print();
		//body.pos.print();
	}
		

	// Calc all forces from box edge
	body.vel = body.vel + calcEdgeForce(&body) * simulation->dt;



	// Integrate position  AFTER ALL BODIES IN BLOCK HAVE BEEN CALCULATED? No should not be a problem as all update their local body, 
	// before moving to shared?? Although make the local just a pointer might be faster, since a SimBody might not fit in thread registers!!!!!
	body.pos = body.pos + body.vel * simulation->dt;

	if (body.molecule_type != UNUSED_BODY ) {
		//printf("Block %d output body at pos: %f %f %f\n", blockID, body.pos.x, body.pos.y, body.pos.z);
	}


	// Correct for bonded-body constraints? Or maybe this should be calculated as a force too??




	// Ensure no position is EXACTLY on the block edge, causing it to be lost forever!




	// Publish new positions for the focus bodies
	
	//block.focus_bodies[bodyID] = body;	// Probably useless, right?
	
	accesspoint.bodies[bodyID] = body;
	__syncthreads();
	if (bodyID == 0) {
		simulation->box->accesspoint[blockID] = accesspoint;
		//simulation->box->blocks[blockID] = block;	// Very expensive, only needed for rendering.	NOT NECESSARY IF WE RUN UPDATEKERNEL AFTER STEP BEFORE RENDER!!
	}
} 









__global__ void updateKernel(Simulation* simulation, int offset) {
	int blockID1 = blockIdx.x + offset;
	Int3 threadID3 = Int3( threadIdx.x, threadIdx.y, threadIdx.z);
	int threadID1 = indexConversion(threadID3, 3);

	if (blockID1 >= simulation->box->n_blocks)
		return;



	
	__shared__ Block block;
	__shared__ Int3 blockID3;
	__shared__ int bpd;
	__shared__ char relationtype[27];	// 0 = far, 1 = near, 2 = focus
	__shared__ SimBody bodybuffer[27];
	
	__syncthreads();
	if (threadID1 == 0) {
		block = simulation->box->blocks[blockID1];
		block.n_bodies = 0;
		bpd = simulation->blocks_per_dim;
		blockID3 = indexConversion(blockID1, bpd);
	}
	__syncthreads();




	Int3 neighbor_index3 = blockID3 + (threadID3-Int3(1,1,1));


	AccessPoint accesspoint;
	int neighbor_index1 = -1;
	if (!(neighbor_index3.x < 0 || neighbor_index3.x >= bpd || neighbor_index3.y < 0 || neighbor_index3.y >= bpd || neighbor_index3.z < 0 || neighbor_index3.z >= bpd)) {
		 neighbor_index1 = indexConversion(neighbor_index3, bpd);
		accesspoint = simulation->box->accesspoint[neighbor_index1];
	}
		


	// Delete all previous bodies 
	if (threadID1 < MAX_FOCUS_BODIES)
		block.focus_bodies[threadID1] = SimBody();
	for (int i = threadID1; i < MAX_NEAR_BODIES + 27; i += 27) 
		if (i < MAX_NEAR_BODIES)
			block.near_bodies[threadID1] = SimBody();
	
	



	
	int focus_ptr = 0;	
	int near_ptr = 0;
	

	for (int i = 0; i < MAX_FOCUS_BODIES; i++) {
		SimBody body = accesspoint.bodies[i];

		bodybuffer[threadID1] = body;
		relationtype[threadID1] = 0;
		relationtype[threadID1] = (bodyInNear(&body, &block.center) + bodyInFocus(&body, &block.center) * 2) * (body.molecule_type != UNUSED_BODY);


		__syncthreads();
		if (threadID1 == 0) {		
			for (int j = 0; j < 27; j++) {
				if (near_ptr > 47 || focus_ptr > 15) {
					printf("\n\n\t\t\t\tOVERFLOW\n");
					printf("Block %d overflowing! near %d    focus %d\n", blockID1, near_ptr, focus_ptr);
					break;
				}
					
				
				if (relationtype[j] > 1) {
					block.focus_bodies[focus_ptr++] = bodybuffer[j];
					//printf("Block %d loading body at pos: %f %f %f\n",blockID1, block.focus_bodies[focus_ptr - 1].pos.x, block.focus_bodies[focus_ptr - 1].pos.y, block.focus_bodies[focus_ptr - 1].pos.z);
					//block.n_bodies++;
				}				
				else if (relationtype[j] == 1)
					block.near_bodies[near_ptr++] = bodybuffer[j];
					
				
			}
		}
		__syncthreads();
	}


	__syncthreads();
	if (threadID1 == 0) {
		//printf("\nBlock %d bodies: %d\n", blockID, block.n_bodies);
		simulation->box->blocks[blockID1] = block;
		//printf("Bodies: %d %d\n", focus_ptr, near_ptr);
	}


}