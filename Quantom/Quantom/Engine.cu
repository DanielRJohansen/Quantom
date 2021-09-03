#include "Engine.cuh"




Simulation* Engine::prepSimulation(Simulation* simulation) {
	this->simulation = simulation;
	simulation->bodies = new SimBody[simulation->n_bodies];


	int n_blocks = initBlocks();
	linkBlocks();

	int n_bodies = fillBox();


	printf("\n\nSimbody size: %d\n", sizeof(SimBody));
	printf("Float3 size: %d", sizeof(Float3));
	printf("Block size: %d calculated size: %d\n", sizeof(Block), sizeof(SimBody) * MAX_BLOCK_BODIES + sizeof(Float3) + 7 * sizeof(int));

	printf("Simulation configured with %d blocks, and %d bodies. \n", n_blocks, n_bodies);
	printf("Required shared mem for stepKernel: %d\n", sizeof(Block));


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
	int blocks_per_dim = simulation->box_size / BLOCK_LEN_CUDA;
	simulation->blocks_per_dim = blocks_per_dim;
	int n_blocks = pow(blocks_per_dim, 3);
	simulation->box->n_blocks = n_blocks;
	simulation->box->blocks = new Block[n_blocks];
	printf("Blocks per dim: %d\n", blocks_per_dim);
	 

	int index = 0;
	float offset = -simulation->box_size / 2 + 0.5 * BLOCK_LEN_CUDA;
	for (int x = 0; x < blocks_per_dim; x++) {
		for (int y = 0; y < blocks_per_dim; y++) {
			for (int z = 0; z < blocks_per_dim; z++) {
				Float3 center(x * BLOCK_LEN_CUDA + offset, y * BLOCK_LEN_CUDA + offset, z * BLOCK_LEN_CUDA + offset);
				simulation->box->blocks[index++] = Block(center);
				//printf("%.1f %.1f %.1f\n", center.x, center.y, center.z);
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

	int index = 0;
	for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
				if (index == simulation->n_bodies)
					break;
				simulation->bodies[index].pos = Float3(base + dist * (float) x_index, base + dist * float(y_index), base + dist * float(z_index));
				simulation->bodies[index].rotation = Float3(0, 0, 0);
				simulation->bodies[index].rot_vel = Float3(0, 0, PI);
				placeBody(&simulation->bodies[index++]);
			}
		}
	}

	return index;
}
	
		
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


void Engine::prepareCudaScheduler() {
	sim_blocks = simulation->box->n_blocks;

	for (int i = 0; i < N_STREAMS; i++)
		cudaStreamCreate(&stream[i]);


	gridblock_size = dim3(GRIDBLOCKS_PER_BODY, BLOCKS_PER_SM, 1);
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
			stepKernel << < gridblock_size, THREADS_PER_GRIDBLOCK, sizeof(int), stream[i] >> > (simulation);
			blocks_handled += BLOCKS_PER_SM;

			if (blocks_handled >= sim_blocks)
				break;
		}
	}

	stepKernel <<< sim_blocks, MAX_BLOCK_BODIES, sizeof(int) >> > (simulation);
	cudaDeviceSynchronize();

	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "step kernel failed!");
		exit(1);
	}



	


	auto stop = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
	printf("Step time: %5.d ys", duration.count());
}





__global__ void stepKernel(Simulation* simulation) {
	int blockID = blockIdx.x;
	int bodyID = threadIdx.x;
	//return;
	
	if (bodyID >= simulation->box->blocks[blockID].n_bodies)
		return;

	


	// Load bodies into shared memory
	__shared__ SimBody bodies[MAX_BLOCK_BODIES];	
	//bodies[bodyID] = simulation->box->blocks[blockID].bodies[bodyID];	// transferring the entrire block is too much for 1 thread. Must find alternative.
	__syncthreads();



	// BEGIN WORK
	Block* block = &simulation->box->blocks[blockID];
	SimBody body = block->bodies[bodyID];
	//SimBody body = bodies[bodyID];

	//body.rot_vel = test;

	body.rotation = body.rotation + body.rot_vel * simulation->dt;				// * dt of course!



	block->bodies[bodyID] = body;




	__syncthreads();
	//if (bodyID == 0)
		//simulation->box->blocks[blockID] = block;

	
} 