#include "BoxBuilder.cuh"



void BoxBuilder::build(Simulation* simulation) {
	simulation->box->compounds = new Compound[MAX_COMPOUNDS];
	//simulation->box->solvents = new Compound[MAX_COMPOUNDS];
	//compoundneighborlists_host = new CompoundNeighborList[MAX_COMPOUNDS];
	//compoundstates_host = new CompoundState[MAX_COMPOUNDS];


	simulation->box->compounds = new Compound[MAX_COMPOUNDS];
	simulation->box->solvents = new Solvent[MAX_SOLVENTS];

	simulation->box->compound_state_array = new CompoundState[MAX_COMPOUNDS];
	cudaMalloc(&simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS);
	cudaMalloc(&simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS);

	simulation->box->solvent_neighborlist_array = new SolventNeighborList[MAX_COMPOUNDS + MAX_SOLVENTS];	// These are divided since both compounds and solvents will be near many more solvents than compounds
	simulation->box->compound_neighborlist_array = new CompoundNeighborList[MAX_COMPOUNDS + MAX_SOLVENTS];



	//cudaMalloc(&simulation->box->compound_state_array, sizeof(CompoundState) * MAX_COMPOUNDS);
	//cudaMalloc(&simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS);
	//cudaMalloc(&simulation->box->compound_neighborlist_array, sizeof(CompoundNeighborList) * MAX_COMPOUNDS);

	//cudaMalloc(&simulation->box->solvents, sizeof(Solvent) * MAX_SOLVENTS);






	placeMainMolecule(simulation);
	solvateBox(simulation);	// Always do after placing compounds


	compoundLinker(simulation);
	solvateLinker(simulation);
	solvateCompoundCrosslinker(simulation);



	
	int n_points = simulation->box->n_compounds * 3 * simulation->n_steps;
	cudaMalloc(&simulation->box->data_buffer, sizeof(float) * n_points);	// Can only log molecules of size 3 for now...
										




	Molecule water;
	for (int i = 0; i < water.n_atoms; i++) {
		simulation->box->rendermolecule.radii[i] = water.atoms[i].radius;
		for (int j = 0; j < 3; j++)
			simulation->box->rendermolecule.colors[i][j] = water.atoms[i].color[j];
	}

	simulation->box->dt = simulation->dt;






	printf("N points: %d\n", simulation->box->n_compounds * 3 * simulation->n_steps);
	cudaMalloc(&simulation->box->trajectory, sizeof(Float3) * simulation->box->n_compounds * 3 * simulation->n_steps);




	simulation->box->moveToDevice();
}


void BoxBuilder::placeMainMolecule(Simulation* simulation) {
	Float3 compound_center = Float3(BOX_LEN_HALF, BOX_LEN_HALF, BOX_LEN_HALF);
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
	simulation->box->solvents = new Solvent[MAX_SOLVENTS];

	


	int bodies_per_dim = ceil(cbrt((float)N_SOLVATE_MOLECULES));
	float dist_between_compounds = (BOX_LEN) / (float)bodies_per_dim;	// dist_per_index
	float base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);


	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (simulation->box->n_solvents == N_SOLVATE_MOLECULES)
					break;

				Float3 solvent_center = Float3(base + dist_between_compounds * (float)x_index, base + dist_between_compounds * (float)y_index, base + dist_between_compounds * (float)z_index);
				float solvent_radius = 0.2;

				if (spaceAvailable(simulation->box, solvent_center, solvent_radius)) {
					simulation->box->solvents[simulation->box->n_solvents++] = createSolvent(
						solvent_center,
						simulation->dt
					);
				}
			}
		}
	}
	printf("%d solvents added to box\n", simulation->box->n_solvents);
	return simulation->box->n_solvents;
}







Compound BoxBuilder::createCompound(Float3 com, int compound_index, CompoundState* statebuffer_node, CompoundNeighborList* neighborinfo_node, float dt) {

	int n_atoms = 3;
	Float3 offsets[3] = { Float3(0,0,0), Float3(0.13, 0, 0), Float3(0, 0, -0.13) };
	for (int i = 0; i < n_atoms; i++) {
		statebuffer_node->positions[i] = com + offsets[i];	// PLACE EACH PARTICLE IN COMPOUNDS STATE, BEFORE CREATING COMPOUNDS, LETS US IMMEDIATELY CALCULATE THE COMPOUNDS CENTER OF MASS.
		statebuffer_node->n_particles++;
	}
	
	float vrms = 250;
	Float3 compound_united_vel = Float3(vrms , 0,0);
	Compound compound(compound_index, statebuffer_node);
	for (int i = 0; i < n_atoms; i++) {
		Float3 atom_pos_sub1 = statebuffer_node->positions[i] - compound_united_vel * dt;
		compound.particles[i] = CompactParticle(12.0107*1e-3, atom_pos_sub1);
		compound.n_particles++;
	}
	return compound;
}



Solvent BoxBuilder::createSolvent(Float3 com, float dt)	// Nodes obv. points to addresses in device global memory.
{
	Float3 solvent_vel = Float3(random(), random(), random()).norm() * v_rms * 0;
	return 	Solvent(com, com - solvent_vel * dt);
}

bool BoxBuilder::spaceAvailable(Box* box, Float3 com, float radius) {	// Too lazy to implement yet..
	for (int i = 0; i < box->n_compounds; i++) {
		for (int j = 0; j < box->compounds[i].n_particles; j++) {
			float dist = (box->compound_state_array[i].positions[j] - com).len();
			if (dist < radius)
				return false;
		}
	}
	return true;
}





void BoxBuilder::compoundLinker(Simulation* simulation) {
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		for (int j = 0; j < simulation->box->n_compounds; j++) {
			if (i != j) {
				simulation->box->compound_neighborlist_array[i].addIndex(j);
				//CompoundNeighborList* nlist = &compoundneighborlists_host[i];
				//CompoundNeighborList* nlist = &simulation->box->compound_neighborlist_array[i];
				//nlist->neighborcompound_indexes[nlist->n_neighbors++] = j;
			}
		}
	}
}

void BoxBuilder::solvateLinker(Simulation* simulation)
{
	for (int i = 0; i < simulation->box->n_solvents; i++) {
		Solvent* self = &simulation->box->solvents[i];
		for (int j = i; j < simulation->box->n_solvents; j++) {
			Solvent* other = &simulation->box->solvents[j];
			if (i != j) {
				if ((self->pos - other->pos).len() < CUTOFF) {
					simulation->box->solvent_neighborlist_array[i + MAX_COMPOUNDS].addIndex(j);
					simulation->box->solvent_neighborlist_array[j + MAX_COMPOUNDS].addIndex(i);
					//self->addNeighbor(j);		// WRONG INDEX
					//other->addNeighbor(i);		// WRONG INDEX
				}
			}
		}
	}
}

void BoxBuilder::solvateCompoundCrosslinker(Simulation* simulation)
{
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		Compound* compound = &simulation->box->compounds[i];
		for (int j = 0; j < simulation->box->n_solvents; j++) {
			Solvent* solvent= &simulation->box->solvents[j];
			if ((compound->center_of_mass - solvent->pos).len() < CUTOFF) {
				simulation->box->solvent_neighborlist_array[i].addIndex(j);
				simulation->box->compound_neighborlist_array[j + MAX_COMPOUNDS].addIndex(i);

				//CompoundNeighborList* nlist = &compoundneighborlists_host[i];
				//nlist->neighborcompound_indexes[nlist->n_neighbors++] = j;
			}
		}
	}
	printf("Compound 0 solvents: %d\n", simulation->box->solvent_neighborlist_array[0].n_neighbors);
}











/*
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
*/