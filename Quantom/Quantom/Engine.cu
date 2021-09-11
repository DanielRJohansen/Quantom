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
	bpd = ((simulation->box_size + FOCUS_LEN)/ block_dist);		// Otherwise no block focus on edge particles
	int n_blocks = pow(bpd, 3);
	box_base = - (BOX_LEN/2.f);			// Is the center of first block!
	//float block_center_base = box_base + block_dist;


	printf("Blocks per dim: %d\n", bpd);

	simulation->blocks_per_dim = bpd;
	simulation->box->blocks_per_dim = bpd;
	simulation->box->n_blocks = n_blocks;
	simulation->box->blocks = new Block[n_blocks];


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
	printf("Bodies per dim: %d\n", bodies_per_dim);
	float dist = simulation->box_size / (float)bodies_per_dim;	// dist_per_index
	float base = -simulation->box_size / 2.f + dist / 2.f;

	float vel_scalar = 0.5;

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

				simulation->bodies[index].pos = Float3(base + dist * (float)x_index, base + dist * float(y_index), base + dist * float(z_index));

				//printf("Body pos: ");
				//simulation->bodies[index].pos.print();

				simulation->bodies[index].molecule_type = 0;
				simulation->bodies[index].vel = Float3(r1 * vel_scalar, r2 * vel_scalar, r3 * vel_scalar);
				simulation->bodies[index].rotation = Float3(0, 0, 0);
				simulation->bodies[index].rot_vel = Float3(0, PI, 0);
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
	int sharedmem_size = sizeof(Block);
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
	


	/*
	stepKernel <<< sim_blocks, MAX_BLOCK_BODIES, sizeof(int) >> > (simulation);
	cudaDeviceSynchronize();
*/

	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	printf("Step time: %d ys", duration.count());
}

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
		((0.5f / dists.x) + 0.5),
		((0.5f / dists.y) + 0.5),
		((0.5f / dists.z) + 0.5)
	);
}

__device__ inline Float3 calcEdgeForce(Block* block, SimBody* body) {
	Float3 dists1(body->pos.x + BOX_LEN / 2.f, body->pos.y + BOX_LEN / 2.f, body->pos.z + BOX_LEN / 2.f);
	Float3 dists2(BOX_LEN / 2.f - body->pos.x, BOX_LEN / 2.f - body->pos.y, BOX_LEN / 2.f - body->pos.z);


	Float3 pos_forces = forceFromDist(dists1).zeroIfBelow(0);
	Float3 neg_forces = forceFromDist(dists2).zeroIfBelow(0);
	return pos_forces - neg_forces;
}


enum Direction
{
	up, down, left, right, back, forward
};

__device__ void transferBody(Box* box, Block* block, SimBody* body) {
	Float3 rel_pos = body->pos - block->center;

	//switch (rel_pos.x)
}

																				// TODO: MAKE IS SUCH THAT A BODY CAN NEVER BE EXACTLY ON THE EDGE OF FOCUS, THUS APPEARING IN MULTIPLE GROUPS!
__device__ bool bodyInNear(SimBody* body, Float3* block_center) {
	Float3 dist_from_center = (body->pos - *block_center).abs();
	return (dist_from_center.x < FOCUS_LEN, dist_from_center.z < FOCUS_LEN, dist_from_center.z < FOCUS_LEN);
}

__device__ bool bodyInFocus(SimBody* body, Float3* block_center) {
	Float3 dist_from_center = (body->pos - *block_center).abs();
	return (dist_from_center.x < FOCUS_LEN_HALF, dist_from_center.z < FOCUS_LEN_HALF, dist_from_center.z < FOCUS_LEN_HALF);
}





__global__ void stepKernel(Simulation* simulation, int offset) {
	int blockID = blockIdx.x + offset;
	int bodyID = threadIdx.x;

	if (blockID >= simulation->box->n_blocks)
		return;
	




	

	


	// Load bodies into shared memory
	__shared__ Block block;	
	if (threadIdx.x == 0)
		block = simulation->box->blocks[blockID];


	int focus_ptr = 0;

	if (threadIdx.x == 0) {

		while (block.focus_bodies[focus_ptr].molecule_type != UNUSED_BODY)
			focus_ptr++;


		for (int i = 0; i < MAX_NEAR_BODIES; i++) {
			SimBody body = block.near_bodies[i];

			if (body.molecule_type == UNUSED_BODY)
				continue;				// Maybe break if always sorted??

			if (focus_ptr == 16)
				printf("FUCK! %d\n", blockID);

			if (bodyInFocus(&body, &block.center)) {
				block.focus_bodies[focus_ptr++] = body;
				block.near_bodies[i].molecule_type = UNUSED_BODY;
			}
		}
		//printf("%d focusses\n", focus_ptr);
	}
	__syncthreads();


	


	// End thread if not needed.
	if (block.focus_bodies[bodyID].molecule_type == UNUSED_BODY && bodyID != 0)		// Always need thread0 to send block global
		return;
		//printf("What the fuckj is going on. n focusses %d\n", focus_ptr);
		//return;


	// BEGIN WORK
	SimBody body = block.focus_bodies[bodyID];
	


	
	//if (block.edge_block)
		//body.vel = body.vel + calcEdgeForce(&block, &body) * 0.1;

	//if (abs(body.pos.x - block.center.x) > SOLOBLOCK_DIST) {
	//}
	


	body.rotation = body.rotation + body.rot_vel * simulation->dt;				// * dt of course!

	//body.pos = body.pos + body.vel * simulation->dt;

	block.focus_bodies[bodyID] = body;
	



	__syncthreads();
	if (bodyID == 0)
		simulation->box->blocks[blockID] = block;	// Very expensive..

	
} 