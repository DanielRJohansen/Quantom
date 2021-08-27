#include "Engine.cuh"



Engine::Engine(Simulation* simulation) : simulation(simulation) {
	simulation->bodies = new SimBody[simulation->n_bodies];

	
	initBlocks();
	fillBox();
}

void Engine::initBlocks() {
	int blocks_per_dim = simulation->box_size / BLOCK_LEN;
	simulation->blocks_per_dim = blocks_per_dim;
	int n_blocks = pow(blocks_per_dim, 3);
	simulation->box->blocks = new Block[n_blocks];
	int index = 0;
	double offset = -simulation->box_size / 2 + 0.5 * BLOCK_LEN;
	for (int x = 0; x < blocks_per_dim; x++) {
		for (int y = 0; y < blocks_per_dim; y++) {
			for (int z = 0; z < blocks_per_dim; z++) {
				Double3 center(x * BLOCK_LEN + offset, y * BLOCK_LEN + offset, z * BLOCK_LEN + offset);
				simulation->box->blocks[index] = Block(center);
			}
		}
	}
}

void Engine::fillBox() {
	int bodies_per_dim = cbrt(simulation->n_bodies);
	double dist = simulation->box_size / double(bodies_per_dim);	// dist_per_index
	
	for (int i = 0; i < simulation->n_bodies; i++) {
		int x_index = i % bodies_per_dim;
		int y_index = (i / bodies_per_dim) % bodies_per_dim;
		int z_index = i / (bodies_per_dim * bodies_per_dim);

		simulation->bodies[i].pos = Double3(dist * double(x_index), dist * double(y_index), dist * double(z_index));
	}
}

void Engine::placeBody(SimBody* body) {
	Int3 block_index = posToBlockIndex(&body->pos);
	


}
