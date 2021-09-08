#include "Engine.cuh"




Simulation* Engine::prepSimulation(Simulation* simulation) {
	this->simulation = simulation;
	simulation->bodies = new SimBody[simulation->n_bodies];


	int n_blocks = initBlocks();
	linkBlocks();
	prepareEdgeBlocks();

	srand(290128301);
	int n_bodies = fillBox();


	printf("\n\nSimbody size: %d bytes\n", sizeof(SimBody));
	printf("Block size: %d calculated size: %d\n", sizeof(Block), sizeof(SimBody) * MAX_BLOCK_BODIES + sizeof(Float3) + 7 * sizeof(int));
	printf("Simulation configured with %d blocks, and %d bodies. Approximately %d bodies per block. \n", n_blocks, n_bodies, n_bodies/n_blocks);
	printf("Required shared mem for stepKernel: %d\n", sizeof(Block));
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
	

	block_dist = (BLOCK_LEN_CUDA - BLOCK_OVERLAP);
	bpd = ceil(simulation->box_size / block_dist);
	box_base = -(bpd) / 2.f * block_dist;
	float block_center_base = box_base + block_dist / 2.f;

	if (BOX_LEN_CUDA == BLOCK_LEN_CUDA)
		block_center_base = 0;
	printf("Blocks per dim: %d\n", bpd);

	simulation->blocks_per_dim = bpd;
	simulation->box->blocks_per_dim = bpd;
	int n_blocks = pow(bpd, 3);
	simulation->box->n_blocks = n_blocks;
	simulation->box->blocks = new Block[n_blocks];

	int index = 0;
	float offset = -simulation->box_size / 2 + 0.5 * BLOCK_LEN_CUDA;
	for (int x = 0; x < bpd; x++) {
		for (int y = 0; y < bpd; y++) {
			for (int z = 0; z < bpd; z++) {
				//Float3 center(x * BLOCK_LEN_CUDA + offset, y * BLOCK_LEN_CUDA + offset, z * BLOCK_LEN_CUDA + offset);
				Float3 center(x * block_dist + block_center_base, y * block_dist + block_center_base, z * block_dist + block_center_base);
				simulation->box->blocks[index++] = Block(center);
				//printf("Box center: %.1f %.1f %.1f\n", center.x, center.y, center.z);
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
				simulation->bodies[index].vel = Float3(r1 * vel_scalar, r2 * vel_scalar, r3 * vel_scalar);
				simulation->bodies[index].rotation = Float3(0, 0, 0);
				simulation->bodies[index].rot_vel = Float3(0, PI, 0);
				placeBody(&simulation->bodies[index++]);
			}
		}
	}
	return index;
}

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

/*
void Engine::placeBody(SimBody* body) {
	Int3 block_index = posToBlockIndex(&body->pos);
	

	int block_index_1d = block3dIndexTo1dIndex(block_index);
	

	Block* block = &simulation->box->blocks[block_index_1d];


	if (block->n_bodies == MAX_BLOCK_BODIES) {
		printf("Too many bodies for this block!");
		exit(1);
	}

	block->bodies[block->n_bodies] = *body;
	block->n_bodies++;
	
}
	*/
void Engine::linkBlocks() {
	for (int x_index = 0; x_index < simulation->blocks_per_dim; x_index++) {
		for (int y_index = 0; y_index < simulation->blocks_per_dim; y_index++) {
			for (int z_index = 0; z_index < simulation->blocks_per_dim; z_index++) {
				int block_index = block3dIndexTo1dIndex(Int3(x_index, y_index, z_index));
				Block* block = &simulation->box->blocks[block_index];

				block->neighbor_indexes[0] = x_index > 0 ? block3dIndexTo1dIndex(Int3(x_index - 1, y_index, z_index)) : NULL;
				block->neighbor_indexes[1] = (x_index + 1) < simulation->blocks_per_dim ? block3dIndexTo1dIndex(Int3(x_index + 1, y_index, z_index)) : NULL;

				block->neighbor_indexes[2] = y_index > 0 ? block3dIndexTo1dIndex(Int3(x_index, y_index - 1, z_index)) : NULL;
				block->neighbor_indexes[3] = (y_index + 1) < simulation->blocks_per_dim ? block3dIndexTo1dIndex(Int3(x_index, y_index + 1, z_index)) : NULL;

				block->neighbor_indexes[4] = z_index > 0 ? block3dIndexTo1dIndex(Int3(x_index, y_index, z_index - 1)) : NULL;
				block->neighbor_indexes[5] = (z_index + 1) < simulation->blocks_per_dim ? block3dIndexTo1dIndex(Int3(x_index, y_index, z_index + 1)) : NULL;
			}
		}
	}
}

void Engine::prepareEdgeBlocks() {
	for (int x_index = 0; x_index < simulation->blocks_per_dim; x_index++) {
		for (int y_index = 0; y_index < simulation->blocks_per_dim; y_index++) {
			for (int z_index = 0; z_index < simulation->blocks_per_dim; z_index++) {
				int block_index = block3dIndexTo1dIndex(Int3(x_index, y_index, z_index));
				Block* block = &simulation->box->blocks[block_index];
				if (x_index == 0 || x_index == (simulation->blocks_per_dim - 1) || y_index == 0 || y_index == (simulation->blocks_per_dim - 1) || z_index == 0 || z_index == (simulation->blocks_per_dim-1))
					block->edge_block = true;

			}
		}
	}
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
			stepKernel << < BLOCKS_PER_SM, MAX_BLOCK_BODIES, sharedmem_size, stream[i] >> > (simulation, blocks_handled);

			

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

	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "step kernel failed!");
		exit(1);
	}
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
	Float3 dists1(body->pos.x + BOX_LEN_CUDA / 2.f, body->pos.y + BOX_LEN_CUDA / 2.f, body->pos.z + BOX_LEN_CUDA / 2.f);
	Float3 dists2(BOX_LEN_CUDA / 2.f - body->pos.x, BOX_LEN_CUDA / 2.f - body->pos.y, BOX_LEN_CUDA / 2.f - body->pos.z);


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



__global__ void stepKernel(Simulation* simulation, int offset) {
	int blockID = blockIdx.x + offset;
	int bodyID = threadIdx.x;
	//return;
	

	


	// Load bodies into shared memory
	__shared__ Block block;	
	if (threadIdx.x == 0)
		block = simulation->box->blocks[blockID];
	__syncthreads();


	Block block1 = block;
	// End thread if not needed.
	if (bodyID >= block.n_bodies)
		return;


	// BEGIN WORK
	SimBody body = block.bodies[bodyID];
	



	if (block.edge_block)
		body.vel = body.vel + calcEdgeForce(&block, &body) * 0.1;

	if (abs(body.pos.x - block.center.x) > SOLOBLOCK_DIST) {

	}



	body.rotation = body.rotation + body.rot_vel * simulation->dt;				// * dt of course!
	body.pos = body.pos + body.vel * simulation->dt;


	block.bodies[bodyID] = body;




	__syncthreads();
	if (bodyID == 0)
		simulation->box->blocks[blockID] = block;	// Very expensive..

	
} 