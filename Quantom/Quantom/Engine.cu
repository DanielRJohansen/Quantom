#include "Engine.cuh"




Simulation* Engine::prepSimulation(Simulation* simulation) {
	this->simulation = simulation;
	simulation->bodies = new SimBody[simulation->n_bodies];


	int n_blocks = initBlocks();
	linkBlocks();

	int n_bodies = fillBox();


	printf("\n\nSimbody size: %d\n", sizeof(SimBody));
	printf("Double3 size: %d", sizeof(Double3));
	printf("Block size: %d calculated size: %d\n", sizeof(Block), sizeof(SimBody) * MAX_BLOCK_BODIES + sizeof(Double3) + 7 * sizeof(int));

	printf("Simulation configured with %d blocks, and %d bodies. \n", n_blocks, n_bodies);
	printf("Required shared mem for stepKernel: %d\n", sizeof(Block));


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
	double offset = -simulation->box_size / 2 + 0.5 * BLOCK_LEN_CUDA;
	for (int x = 0; x < blocks_per_dim; x++) {
		for (int y = 0; y < blocks_per_dim; y++) {
			for (int z = 0; z < blocks_per_dim; z++) {
				Double3 center(x * BLOCK_LEN_CUDA + offset, y * BLOCK_LEN_CUDA + offset, z * BLOCK_LEN_CUDA + offset);
				simulation->box->blocks[index++] = Block(center);
				//printf("%.1f %.1f %.1f\n", center.x, center.y, center.z);
			}
		}
	}
	return index;
}
		
int Engine::fillBox() {
	int bodies_per_dim = ceil(cbrt((double)simulation->n_bodies));
	printf("Bodies per dim: %d\n", bodies_per_dim);
	double dist = simulation->box_size / (double)bodies_per_dim;	// dist_per_index
	double base = -simulation->box_size / 2.f + dist / 2.f;

	int index = 0;
	for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
				if (index == simulation->n_bodies)
					break;
				simulation->bodies[index].pos = Double3(base + dist * (double) x_index, base + dist * double(y_index), base + dist * double(z_index));
				simulation->bodies[index].rotation = Double3(0, 0, 0);
				simulation->bodies[index].rot_vel = Double3(0, 0, PI);
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





		

	


//--------------------------------------------------------------------------	SIMULATION BEGINS HERE --------------------------------------------------------------//


void Engine::step() {
	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "Error before step!");
		exit(1);
	}

	auto start = std::chrono::high_resolution_clock::now();


	int n_gpublocks = simulation->box->n_blocks;
	int n_blockthreads = MAX_BLOCK_BODIES;
	//int sharedmem_bytesize = sizeof(SimBody) * MAX_BLOCK_BODIES;
	int sharedmem_bytesize = sizeof(int);
	//int sharedmem_bytesize = sizeof(Double3);

	//printf("Simulation total steps host = %d\n", simulation->n_steps);
	stepKernel <<< n_gpublocks, n_blockthreads, sharedmem_bytesize>> > (simulation);
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