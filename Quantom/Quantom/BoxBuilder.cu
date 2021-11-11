#include "BoxBuilder.cuh"



void BoxBuilder::build(Simulation* simulation) {
	simulation->box->compounds = new Compound_H2O[max_compounds];
	compoundneighborlists_host = new CompoundNeighborList[max_compounds];
	compoundstates_host = new CompoundState[max_compounds];

	cudaMalloc(&simulation->box->compound_state_array, sizeof(CompoundState) * max_compounds);
	cudaMalloc(&simulation->box->compound_neighborlist_array, sizeof(CompoundNeighborList) * max_compounds);

	simulation->box->n_compounds = solvateBox(simulation);

	compoundLinker(simulation);

	cudaMemcpy(simulation->box->compound_state_array, compoundstates_host, sizeof(CompoundState) * max_compounds, cudaMemcpyHostToDevice);
	cudaMemcpy(simulation->box->compound_neighborlist_array, compoundneighborlists_host, sizeof(CompoundNeighborList) * max_compounds, cudaMemcpyHostToDevice);




	// FUCK THIS SHIT TO ABSOLUTE HEELLLL
	int n_points = simulation->box->n_compounds * 3 * 2 * 10000;
	cudaMalloc(&simulation->box->data_buffer, sizeof(float) * n_points);	// Can only log molecules of size 3 for now...
	float* aaa = new float[n_points];
	for (int i = 0; i < n_points; i++)
		aaa[i] = 999999 * 0;
	cudaMemcpy(simulation->box->data_buffer, aaa, sizeof(float) * n_points, cudaMemcpyHostToDevice);
	cudaDeviceSynchronize();
	delete(aaa);






	Molecule water;
	for (int i = 0; i < water.n_atoms; i++) {
		simulation->box->rendermolecule.radii[i] = water.atoms[i].radius;
		for (int j = 0; j < 3; j++)
			simulation->box->rendermolecule.colors[i][j] = water.atoms[i].color[j];
	}

	simulation->box->dt = simulation->dt;

	simulation->box->moveToDevice();
}


int BoxBuilder::solvateBox(Simulation* simulation)
{
	int bodies_per_dim = ceil(cbrt((float)simulation->n_bodies));
	float dist_between_compounds = (BOX_LEN) / (float)bodies_per_dim;	// dist_per_index
	float base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);


	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (simulation->box->n_compounds == N_BODIES_START)
					break;

				Float3 compound_center = Float3(base + dist_between_compounds * (float)x_index, base + dist_between_compounds * (float)y_index, base + dist_between_compounds * (float)z_index);
				float compound_radius = 0.2;

				if (spaceAvailable(compound_center, compound_radius)) {
					simulation->box->compounds[simulation->box->n_compounds++] = createCompound(	compound_center, 
																								simulation->box->n_compounds, 
																								&simulation->box->compound_state_array[simulation->box->n_compounds],
																								&simulation->box->compound_neighborlist_array[simulation->box->n_compounds]
					);
				}
			}
		}
	}
	return simulation->box->n_compounds;
}









Compound_H2O BoxBuilder::createCompound(Float3 com, int compound_index, CompoundState* statebuffer_node, CompoundNeighborList* neighborinfo_node)	// Nodes obv. points to addresses in device global memory.
{
	Molecule water;
	for (int i = 0; i < water.n_atoms; i++) {
		//(com + water.atoms[i].pos).print('p');
		compoundstates_host[compound_index].positions[i] = com + water.atoms[i].pos;	// PLACE EACH PARTICLE IN COMPOUNDS STATE, BEFORE CREATING COMPOUNDS, LETS US IMMEDIATELY CALCULATE THE COMPOUNDS CENTER OF MASS.
		compoundstates_host[compound_index].n_particles++;
	}


	Float3 compound_united_vel = Float3(random(), random(), random()).norm() * mean_velocity;

	Compound_H2O compound(compound_index, neighborinfo_node, statebuffer_node, compoundstates_host);
	for (int i = 0; i < water.n_atoms; i++) {
		compound.particles[i] = CompactParticle(water.atoms[i].mass, compound_united_vel);
	}
	return compound;
}

bool BoxBuilder::spaceAvailable(Float3 com, float radius) {	// Too lazy to implement yet..
	return true;
}





void BoxBuilder::compoundLinker(Simulation* simulation) {
	for (int i = 0; i < max_compounds; i++) {
		compoundneighborlists_host[i].n_neighbors = 0;
	}

	for (int i = 0; i < simulation->box->n_compounds; i++) {
		for (int j = 0; j < simulation->box->n_compounds; j++) {
			if (i != j) {
				CompoundNeighborList* nlist = &compoundneighborlists_host[i];
				nlist->neighborcompound_indexes[nlist->n_neighbors++] = j;
			}
		}
	}
}
