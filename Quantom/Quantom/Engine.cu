#include "Engine.cuh"



Engine::Engine(Simulation* simulation) : simulation(simulation) {
	simulation->bodies = new SimBody[simulation->n_bodies];

	
	int n_blocks = initBlocks();
	//printf("%d blocks\n", n_blocks);
	int n_bodies = fillBox();

	printf("Simulation configured with %d blocks, and %d bodies\n", n_blocks, n_bodies);
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
	

	int block_index_1d = block_index.z + block_index.y * simulation->blocks_per_dim + block_index.x * simulation->blocks_per_dim * simulation->blocks_per_dim;	// NOTICE THAT X IS THE "GRANDPARENT" ITERATOR
	

	Block* block = &simulation->box->blocks[block_index_1d];

	if ((body->pos - Double3(-0.9, -0.9, 0.1)).len() < 0.001) {
		printf("\n\n\n\nBody pos: %f %f %f\n", body->pos.x, body->pos.y, body->pos.z);
		printf("block index: %d %d %d\n", block_index.x, block_index.y, block_index.z);
		printf("1d index: %d\n", block_index_1d);
		printf("block center: %.2f %.2f %.2f\n", block->center.x, block->center.y, block->center.z);
		printf(" \n\n\n");
	}
	if (block->n_bodies == MAX_BLOCK_BODIES) {
		printf("Too many bodies for this block!");
		exit(1);
	}

	block->bodies[block->n_bodies] = *body;
	block->n_bodies++;
	
}
