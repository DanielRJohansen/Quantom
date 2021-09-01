#include "Engine.cuh"



Engine::Engine(Simulation* simulation) : simulation(simulation) {
	simulation->bodies = new SimBody[simulation->n_bodies];

	
	int n_blocks = initBlocks();
	linkBlocks();
	
	int n_bodies = fillBox();


	printf("Simbody size: %d\n", sizeof(SimBody));
	printf("Block size: %d calculated size: %d\n", sizeof(Block), sizeof(SimBody) * MAX_BLOCK_BODIES + sizeof(Double3) + 7 * sizeof(int));

	printf("Simulation configured with %d blocks, and %d bodies. \n", n_blocks, n_bodies);
	printf("Required shared mem for stepKernel: %d\n", sizeof(Block));
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

	int n_gpublocks = simulation->box->n_blocks;
	int n_blockthreads = MAX_BLOCK_BODIES;
	int sharedmem_bytesize = sizeof(Block);
	stepKernel <<< n_gpublocks, n_blockthreads, sharedmem_bytesize>> > (simulation);
	cudaDeviceSynchronize();

	cuda_status = cudaGetLastError();
	if (cuda_status != cudaSuccess) {
		fprintf(stderr, "step kernel failed!");
		exit(1);
	}


}

__global__ void stepKernel(Simulation* simulation) {
	int blockID = blockIdx.x;
	int bodyID = threadIdx.x;

	if (bodyID >= simulation->box->blocks[blockID].n_bodies)
		return;

	__shared__ Block block;
	if (bodyID == 0)
		block = simulation->box->blocks[blockID];

	__syncthreads();
} 