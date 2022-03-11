#include "BoxBuilder.cuh"



void BoxBuilder::buildBox(Simulation* simulation) {

	simulation->box->compounds = new Compound[MAX_COMPOUNDS];
	simulation->box->solvents = new Solvent[MAX_SOLVENTS];

	simulation->box->compound_state_array = new CompoundState[MAX_COMPOUNDS];
	cudaMalloc(&simulation->box->compound_state_array_next, sizeof(CompoundState) * MAX_COMPOUNDS);
	cudaMalloc(&simulation->box->solvents_next, sizeof(Solvent) * MAX_SOLVENTS);




	simulation->box->solvent_neighborlists = new NeighborList[MAX_SOLVENTS];	
	simulation->box->compound_neighborlists = new NeighborList[MAX_COMPOUNDS];
	for (int i = 0; i < MAX_COMPOUNDS; i++) {
		simulation->box->compound_neighborlists->init();
	}
	for (int i = 0; i < MAX_SOLVENTS; i++) {
		simulation->box->solvent_neighborlists->init();
	}



	Molecule water;
	for (int i = 0; i < water.n_atoms; i++) {
		simulation->box->rendermolecule.radii[i] = water.atoms[i].radius;
		for (int j = 0; j < 3; j++)
			simulation->box->rendermolecule.colors[i][j] = water.atoms[i].color[j];
	}

	simulation->box->dt = simulation->dt;


}

void BoxBuilder::addSingleMolecule(Simulation* simulation, Compound* molecule)
{
	Float3 compound_center = Float3(BOX_LEN_HALF, BOX_LEN_HALF, BOX_LEN_HALF);

	Float3 offset = compound_center - molecule->getCOM();
	for (int i = 0; i < molecule->n_particles; i++) {
		molecule->particles[i].pos_tsub1 += offset;
	}


	integrateCompound(
		molecule,
		simulation
	);
}

void BoxBuilder::addScatteredMolecules(Simulation* simulation, Compound* molecule, int n_copies)
{
	placeMultipleCompoundsRandomly(simulation, molecule, N_LIPID_COPIES);
	printf("Scattered %d Compounds in box\n", simulation->box->n_compounds);
}

void BoxBuilder::finishBox(Simulation* simulation)
{
	printf("%d Compounds in box\n\n", simulation->box->n_compounds);
	solvateBox(simulation);	// Always do after placing compounds


	// Need this variable both on host and device
	simulation->total_particles_upperbound = simulation->box->n_compounds * MAX_COMPOUND_PARTICLES + simulation->box->n_solvents;											// BAD AMBIGUOUS AND WRONG CONSTANTS
	simulation->box->total_particles_upperbound = simulation->total_particles_upperbound;											// BAD AMBIGUOUS AND WRONG CONSTANTS

		/*
	compoundLinker(simulation);
	solvateLinker(simulation);
	solvateCompoundCrosslinker(simulation);
	*/






	cudaMemcpy(simulation->box->compound_state_array_next, simulation->box->compound_state_array, sizeof(CompoundState) * MAX_COMPOUNDS, cudaMemcpyHostToDevice);	// Just make sure they have the same n_particles info...
	cudaMemcpy(simulation->box->solvents_next, simulation->box->solvents, sizeof(Solvent) * MAX_SOLVENTS, cudaMemcpyHostToDevice);



	

	int n_points = simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER;
	printf("Malloc %d MB on device for data buffers\n", (sizeof(double) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER + sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER) / 1000000);
	cudaMallocManaged(&simulation->box->potE_buffer, sizeof(double) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER);	// Can only log molecules of size 3 for now...
	simulation->potE_buffer = new double[simulation->total_particles_upperbound * simulation->n_steps];

	cudaMallocManaged(&simulation->box->traj_buffer, sizeof(Float3) * simulation->total_particles_upperbound * STEPS_PER_LOGTRANSFER);
	simulation->traj_buffer = new Float3[simulation->total_particles_upperbound * simulation->n_steps];

	simulation->temperature_buffer = new float[SIMULATION_STEPS / STEPS_PER_THERMOSTAT + 1];





	long double total_bytes = sizeof(double) * 10 * simulation->n_steps
		+ sizeof(Float3) * 6 * MAX_COMPOUND_PARTICLES * STEPS_PER_TRAINDATATRANSFER;

	printf("Reserving %.2lf GB device mem for logging\n", (total_bytes) * 1e-9);
	cudaMallocManaged(&simulation->box->outdata, sizeof(double) * 10 * simulation->n_steps);	// 10 data streams for 10k steps. 1 step at a time.
	//cudaMallocManaged(&simulation->box->data_GAN, sizeof(Float3) * 6 * MAX_COMPOUND_PARTICLES * simulation->n_steps);
	cudaMallocManaged(&simulation->box->data_GAN, sizeof(Float3) * 6 * MAX_COMPOUND_PARTICLES * STEPS_PER_TRAINDATATRANSFER);
	simulation->traindata_buffer = new Float3[6 * MAX_COMPOUND_PARTICLES * simulation->n_steps];

	cudaDeviceSynchronize();
	if (cudaGetLastError() != cudaSuccess) {
		fprintf(stderr, "Error during log-data allocation\n");
		exit(1);
	}


	simulation->copyBoxVariables();
	simulation->box->moveToDevice();
	printf("Boxbuild complete!\n\n\n");
}



void BoxBuilder::placeSingleMolecule(Simulation* simulation, Compound* main_compound)
{

}

int BoxBuilder::solvateBox(Simulation* simulation)
{
	simulation->box->solvents = new Solvent[MAX_SOLVENTS];

	


	int bodies_per_dim = ceil(cbrt((double)N_SOLVATE_MOLECULES));
	double dist_between_compounds = (BOX_LEN) / (double)bodies_per_dim;	// dist_per_index
	double base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);


	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (simulation->box->n_solvents == N_SOLVATE_MOLECULES)
					break;

				Float3 solvent_center = Float3(base + dist_between_compounds * (double)x_index, base + dist_between_compounds * (double)y_index, base + dist_between_compounds * (double)z_index);
				double solvent_radius = 0.2;

				if (spaceAvailable(simulation->box, solvent_center)) {
					simulation->box->solvents[simulation->box->n_solvents++] = createSolvent(
						solvent_center,
						simulation->dt
					);
				}
				/*
				if (spaceAvailable(simulation->box, solvent_center, solvent_radius)) {
					simulation->box->solvents[simulation->box->n_solvents++] = createSolvent(
						solvent_center,
						simulation->dt
					);
				}*/
			}
		}
	}
	printf("%d solvents added to box\n", simulation->box->n_solvents);
	return simulation->box->n_solvents;
}






/*
void BoxBuilder::integrateCompound(Float3 com, int compound_index, CompoundState* statebuffer_node, double dt, Simulation* simulation) {

	int n_atoms = PARTICLES_PER_COMPOUND;
	Float3 offsets[3] = { Float3(0,0,0), Float3(0.13, 0, 0), Float3(0, 0, -0.13) };
	for (int i = 0; i < n_atoms; i++) {
		statebuffer_node->positions[i] = com + offsets[i];	// PLACE EACH PARTICLE IN COMPOUNDS STATE, BEFORE CREATING COMPOUNDS, LETS US IMMEDIATELY CALCULATE THE COMPOUNDS CENTER OF MASS.
		statebuffer_node->n_particles++;
	}
	
	//double vrms = 250;

	Float3 compound_united_vel = Float3(random(), random(), random()).norm() * v_rms * 0.1;
	Compound compound(compound_index, statebuffer_node);
	for (int i = 0; i < n_atoms; i++) {
		Float3 atom_pos_sub1 = statebuffer_node->positions[i] - compound_united_vel * dt;
		compound.particles[i] = CompactParticle(COMPOUNDPARTICLE_MASS, atom_pos_sub1);
		compound.n_particles++;
	}
	simulation->box->compounds[simulation->box->n_compounds++] = compound;

	//return compound;
}
*/
void BoxBuilder::integrateCompound(Compound* compound, Simulation* simulation)
{
	compound->init();
	CompoundState* state = &simulation->box->compound_state_array[simulation->box->n_compounds];
	Float3 compound_united_vel = Float3(random(), random(), random()).norm() * v_rms * 0;

	for (int i = 0; i < compound->n_particles; i++) {
		state->positions[i] = compound->particles[i].pos_tsub1;
		state->n_particles++;
	}


	for (int i = 0; i < compound->n_particles; i++) {
		Float3 atom_pos_sub1 = state->positions[i] - compound_united_vel * simulation->dt;
		compound->particles[i].pos_tsub1 = atom_pos_sub1;												// Overwrite prev pos here, since we have assigned the former prev pos to the state buffer.
	}

	simulation->box->compounds[simulation->box->n_compounds++] = *compound;

	//return compound;
}



Solvent BoxBuilder::createSolvent(Float3 com, double dt)	// Nodes obv. points to addresses in device global memory.
{
	Float3 solvent_vel = Float3(random(), random(), random()).norm() * v_rms;
	return Solvent(com, com - solvent_vel * dt);
}






/*
void BoxBuilder::compoundLinker(Simulation* simulation) {
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		for (int j = i+1; j < simulation->box->n_compounds; j++) {
			simulation->box->compound_neighborlists[i].addId(j, NeighborList::NEIGHBOR_TYPE::COMPOUND);
			simulation->box->compound_neighborlists[j].addId(i, NeighborList::NEIGHBOR_TYPE::COMPOUND);
		}
	}
}

void BoxBuilder::solvateLinker(Simulation* simulation)
{
	for (int i = 0; i < simulation->box->n_solvents; i++) {
		Solvent* self = &simulation->box->solvents[i];
		for (int j = i; j < simulation->box->n_solvents; j++) {														// +1 here!
			Solvent* other = &simulation->box->solvents[j];
			if (i != j) {
				if ((self->pos - other->pos).len() < (CUTOFF)) {
					simulation->box->solvent_neighborlists[i].addId(j, NeighborList::NEIGHBOR_TYPE::SOLVENT);
					simulation->box->solvent_neighborlists[j].addId(i, NeighborList::NEIGHBOR_TYPE::SOLVENT);
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
				simulation->box->compound_neighborlists[i].addId(j, NeighborList::NEIGHBOR_TYPE::SOLVENT);
				simulation->box->solvent_neighborlists[j].addId(i, NeighborList::NEIGHBOR_TYPE::COMPOUND);
			}
		}
	}
	printf("Compound n0 solvents: %d\n", simulation->box->compound_neighborlists[0].n_solvent_neighbors);
}

*/

void BoxBuilder::placeMultipleCompoundsRandomly(Simulation* simulation, Compound* template_compound, int n_copies)
{
	int copies_placed = 0;
	while (copies_placed < n_copies) {
		Compound* c = randomizeCompound(template_compound);

		if (spaceAvailable(simulation->box, c)) {
			integrateCompound(c, simulation);
			copies_placed++;
		}
		delete c;
	}


	// Temporary check that no to molecules placed are colliding.
	for (int i = 0; i < simulation->box->n_compounds; i++) {
		Compound* c = &simulation->box->compounds[i];
		for (int ii = 0; ii < simulation->box->n_compounds; ii++) {
			Compound* c2 = &simulation->box->compounds[ii];
			if (ii != i) {
				if (!verifyPairwiseParticleMindist(c, c2)) {
					printf("Illegal compound positioning %d %d\n", i, ii);
					exit(0);
				}				
			}
		}
	}
}

Compound* BoxBuilder::randomizeCompound(Compound* original_compound)
{
	Compound* compound = new Compound;
	*compound = *original_compound;

	Float3 xyz_rot = get3Random() * (2.f*PI);
	//rotateCompound(compound, xyz_rot);


	Float3 xyz_target = (get3Random() * 0.6 + Float3(0.2))* BOX_LEN;
	Float3 xyz_mov = xyz_target - calcCompoundCom(original_compound);
	moveCompound(compound, xyz_mov);

	return compound;
}

void BoxBuilder::moveCompound(Compound* compound, Float3 vector)
{
	for (int i = 0; i < compound->n_particles; i++)
		compound->particles[i].pos_tsub1 += vector;
}

Float3 BoxBuilder::calcCompoundCom(Compound* compound)
{
	Float3 com = Float3(0, 0, 0);
	for (int i = 0; i < compound->n_particles; i++)
		com += compound->particles[i].pos_tsub1;
	com *= (1.f / compound->n_particles);

	return com;
}

void BoxBuilder::rotateCompound(Compound* compound, Float3 xyz_rot)
{
	Float3 vec_to_origo = Float3(0, 0, 0) - calcCompoundCom(compound);
	moveCompound(compound, vec_to_origo);

	for (int i = 0; i < compound->n_particles; i++)
		compound->particles[i].pos_tsub1.rotateAroundOrigo(xyz_rot);

	moveCompound(compound, vec_to_origo * -1);
}

BoundingBox BoxBuilder::calcCompoundBoundingBox(Compound* compound)
{
	BoundingBox bb(Float3(9999, 9999, 9999), Float3(-9999, -9999, -9999));
	for (int i = 0; i < compound->n_particles; i++) {
		Float3 pos = compound->particles[i].pos_tsub1;
		for (int dim = 0; dim < 3; dim++) {
			*bb.min.placeAt(dim) = min(bb.min.at(dim), pos.at(dim));
			*bb.max.placeAt(dim) = max(bb.max.at(dim), pos.at(dim));
		}
	}
	return bb;
}

bool BoxBuilder::spaceAvailable(Box* box, Compound* compound)
{
	BoundingBox bb_a = calcCompoundBoundingBox(compound);
	bb_a.addPadding(MIN_NONBONDED_DIST);
	for (int c_index = 0; c_index < box->n_compounds; c_index++) {
		BoundingBox bb_b = calcCompoundBoundingBox(&box->compounds[c_index]);

		/*if (box->n_compounds == 13) {
			printf("\n c index %d\n", c_index);
			bb_a.min.print('a');
			bb_a.max.print('a');
			bb_b.min.print('b');
			bb_b.max.print('b');
		}*/

		if (bb_a.intersects(bb_b)) {
			//printf("Verifying %d %d\n", box->n_compounds, c_index);
			if (!verifyPairwiseParticleMindist(compound, &box->compounds[c_index]))
				return false;			
		}
	}
	return true;
}

bool BoxBuilder::spaceAvailable(Box* box, Float3 particle_center)
{
	for (int c_index = 0; c_index < box->n_compounds; c_index++) {
		BoundingBox bb = calcCompoundBoundingBox(&box->compounds[c_index]);
		bb.addPadding(MIN_NONBONDED_DIST);
		if (bb.pointIsInBox(particle_center)) {
			return false;
		}
	}
	return true;
}

bool BoxBuilder::verifyPairwiseParticleMindist(Compound* a, Compound* b)
{
	for (int ia = 0; ia < a->n_particles; ia++) {
		for (int ib = 0; ib < b->n_particles; ib++) {
			Float3 pos_a = a->particles[ia].pos_tsub1;
			Float3 pos_b = b->particles[ib].pos_tsub1;

			float dist = (pos_a - pos_b).len();
			if (dist < MIN_NONBONDED_DIST)
				return false;
		}
	}
	return true;
}











/*
int BoxBuilder::solvateBox(Simulation* simulation)
{
	int bodies_per_dim = ceil(cbrt((double)N_SOLVATE_MOLECULES));
	double dist_between_compounds = (BOX_LEN) / (double)bodies_per_dim;	// dist_per_index
	double base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);


	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (simulation->box->n_compounds == N_SOLVATE_MOLECULES)
					break;

				Float3 compound_center = Float3(base + dist_between_compounds * (double)x_index, base + dist_between_compounds * (double)y_index, base + dist_between_compounds * (double)z_index);
				double compound_radius = 0.2;

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