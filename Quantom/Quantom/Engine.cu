#include "Engine.cuh"



Engine::Engine(Simulation* simulation) : simulation(simulation) {
	simulation->bodies = new SimBody[simulation->n_bodies];
	fillBox();


}

void Engine::fillBox() {
	int bodies_per_dim = cbrt(simulation->n_bodies);
	double dist = simulation->box_size.x / double(bodies_per_dim);	// dist_per_index
	for (int i = 0; i < simulation->n_bodies; i++) {
		int x_index = i % bodies_per_dim;
		int y_index = (i / bodies_per_dim) % bodies_per_dim;
		int z_index = i / (bodies_per_dim * bodies_per_dim);

		simulation->bodies[i].pos = Double3(dist * double(x_index), dist * double(y_index), dist * double(z_index));
	}
}