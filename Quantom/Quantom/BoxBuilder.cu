#include "BoxBuilder.cuh"



void BoxBuilder::build(Simulation* simulation) {
	simulation->box->compounds = new Compound[max_compounds];
	//simulation->box->solvents = new Compound[max_compounds];
	compoundneighborlists_host = new CompoundNeighborList[max_compounds];
	compoundstates_host = new CompoundState[max_compounds];

	cudaMalloc(&simulation->box->compound_state_array, sizeof(CompoundState) * max_compounds);
	cudaMalloc(&simulation->box->compound_state_array_next, sizeof(CompoundState) * max_compounds);
	cudaMalloc(&simulation->box->compound_neighborlist_array, sizeof(CompoundNeighborList) * max_compounds);





	placeMainMolecule(simulation);
	placeMainMolecule(simulation);

	//simulation->box->n_compounds += solvateBox(simulation);

	compoundLinker(simulation);

	cudaMemcpy(simulation->box->compound_state_array, compoundstates_host, sizeof(CompoundState) * max_compounds, cudaMemcpyHostToDevice);
	cudaMemcpy(simulation->box->compound_neighborlist_array, compoundneighborlists_host, sizeof(CompoundNeighborList) * max_compounds, cudaMemcpyHostToDevice);




	
	int n_points = simulation->box->n_compounds * 3 * simulation->n_steps;
	cudaMalloc(&simulation->box->data_buffer, sizeof(float) * n_points);	// Can only log molecules of size 3 for now...
										




	Molecule water;
	for (int i = 0; i < water.n_atoms; i++) {
		simulation->box->rendermolecule.radii[i] = water.atoms[i].radius;
		for (int j = 0; j < 3; j++)
			simulation->box->rendermolecule.colors[i][j] = water.atoms[i].color[j];
	}

	simulation->box->dt = simulation->dt;






	printf("N points: %d", simulation->box->n_compounds * 3 * simulation->n_steps);
	cudaMalloc(&simulation->box->trajectory, sizeof(Float3) * simulation->box->n_compounds * 3 * simulation->n_steps);




	simulation->box->moveToDevice();
}


void BoxBuilder::placeMainMolecule(Simulation* simulation) {
	Float3 compound_center = Float3(BOX_LEN_HALF * 0.5 + BOX_LEN_HALF * simulation->box->n_compounds, BOX_LEN_HALF, BOX_LEN_HALF);
	float compound_radius = 0.2;

	simulation->box->compounds[simulation->box->n_compounds++] = createCompound(
		compound_center,
		simulation->box->n_compounds,
		&simulation->box->compound_state_array[simulation->box->n_compounds],
		&simulation->box->compound_neighborlist_array[simulation->box->n_compounds],
		simulation->dt
	);
}

int BoxBuilder::solvateBox(Simulation* simulation)
{
	int bodies_per_dim = ceil(cbrt((float)N_SOLVATE_MOLECULES));
	float dist_between_compounds = (BOX_LEN) / (float)bodies_per_dim;	// dist_per_index
	float base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);


	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (simulation->box->n_compounds == N_SOLVATE_MOLECULES)
					break;

				Float3 compound_center = Float3(base + dist_between_compounds * (float)x_index, base + dist_between_compounds * (float)y_index, base + dist_between_compounds * (float)z_index);
				float compound_radius = 0.2;

				if (spaceAvailable(compound_center, compound_radius)) {
					simulation->box->compounds[simulation->box->n_compounds++] = createSolvent(
						compound_center,
						simulation->box->n_compounds,
						&simulation->box->compound_state_array[simulation->box->n_compounds],
						&simulation->box->compound_neighborlist_array[simulation->box->n_compounds],
						simulation->dt
					);
				}
			}
		}
	}
	return simulation->box->n_compounds;
}







Compound BoxBuilder::createCompound(Float3 com, int compound_index, CompoundState* statebuffer_node, CompoundNeighborList* neighborinfo_node, float dt) {
	compoundstates_host[compound_index].positions[0] = com;	// PLACE EACH PARTICLE IN COMPOUNDS STATE, BEFORE CREATING COMPOUNDS, LETS US IMMEDIATELY CALCULATE THE COMPOUNDS CENTER OF MASS.
	compoundstates_host[compound_index].n_particles++;
	


	//Float3 compound_united_vel = Float3(200 - 400 * compound_index, random(), random()).norm() * v_rms;
	Float3 compound_united_vel = Float3(500 - 1000 * compound_index, 0,0) * 2;


	Compound compound(compound_index, neighborinfo_node, statebuffer_node, compoundstates_host);
	for (int i = 0; i < 1; i++) {
		Float3 atom_pos_sub1 = compoundstates_host[compound_index].positions[i] - compound_united_vel * dt;
		compound.particles[i] = CompactParticle(12.0107, atom_pos_sub1);
	}
	return compound;
}



Compound BoxBuilder::createSolvent(Float3 com, int compound_index, CompoundState* statebuffer_node, CompoundNeighborList* neighborinfo_node, float dt)	// Nodes obv. points to addresses in device global memory.
{

	Molecule water;
	for (int i = 0; i < water.n_atoms; i++) {
		//(com + water.atoms[i].pos).print('p');
		compoundstates_host[compound_index].positions[i] = com + water.atoms[i].pos;	// PLACE EACH PARTICLE IN COMPOUNDS STATE, BEFORE CREATING COMPOUNDS, LETS US IMMEDIATELY CALCULATE THE COMPOUNDS CENTER OF MASS.
		compoundstates_host[compound_index].n_particles++;
	}


	//Float3 compound_united_vel = Float3(random(), random(), random()).norm() * mean_velocity;
	Float3 compound_united_vel = Float3(random(), random(), random()).norm() * v_rms * 0;

	printf("Velocity: %f nm/ns\n", compound_united_vel.len());

	Compound compound(compound_index, neighborinfo_node, statebuffer_node, compoundstates_host);
	for (int i = 0; i < water.n_atoms; i++) {
		Float3 atom_pos_sub1 = compoundstates_host[compound_index].positions[i] - compound_united_vel * dt;
		compound.particles[i] = CompactParticle(water.atoms[i].mass, atom_pos_sub1);
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
