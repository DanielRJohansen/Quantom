#include "BoxBuilder.cuh"



void BoxBuilder::build(Simulation* simulation) {
	this->simulation = simulation;
	simulation->box->n_compounds = solvateBox();

	Molecule water;
	for (int i = 0; i < water.n_atoms; i++) {
		simulation->box->rendermolecule.radii[i] = water.atoms[i].radius;
		for (int j = 0; j < 3; j++)
			simulation->box->rendermolecule.colors[i][j] = water.atoms[i].color[j];
	}

	simulation->box->moveToDevice();
}


int BoxBuilder::solvateBox()
{
	int bodies_per_dim = ceil(cbrt((float)simulation->n_bodies));
	float dist_between_compounds = (BOX_LEN) / (float)bodies_per_dim;	// dist_per_index
	float base = box_base + dist_between_compounds / 2.f;
	printf("Bodies per dim: %d. Dist per dim: %.3f\n", bodies_per_dim, dist_between_compounds);



	double m = 18.01528;		// g/mol
	double k_B = 8.617333262145 * 10e-5;
	double T = 293;	// Kelvin
	float mean_velocity = m / (2 * k_B * T);

	const int max_compounds = 1'000'000;
	simulation->box->compounds = new Compound_H2O[max_compounds];
	cudaMalloc(&simulation->box->compound_state_buffer, sizeof(CompoundState) * max_compounds);
	cudaMalloc(&simulation->box->compound_neighborinfo_buffer, sizeof(CompoundNeighborInfo) * max_compounds);
	


	int solvate_cnt = 0;
	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (solvate_cnt == N_BODIES_START)
					break;

				Float3 compound_center = Float3(base + dist_between_compounds * (float)x_index, base + dist_between_compounds * (float)y_index, base + dist_between_compounds * (float)z_index);
				Float3 compound_united_vel = Float3(random(), random(), random()).norm() * mean_velocity;
				float compound_radius = 0.2;

				if (spaceAvailable(compound_center, compound_radius)) {
					Compound_H2O compound = createCompound(compound_center, solvate_cnt, &simulation->box->compound_state_buffer[solvate_cnt], &simulation->box->compound_neighborinfo_buffer[solvate_cnt]);
					simulation->box->compounds[simulation->box->n_compounds++] = compound;
					solvate_cnt++;
				}
			}
		}
	}
	return solvate_cnt;
}












Compound_H2O BoxBuilder::createCompound(Float3 com, int compound_index, CompoundState* statebuffer_node, CompoundNeighborInfo* neighborinfo_node)	// Nodes obv. points to addresses in device global memory.
{
	Molecule water;
	Compound_H2O compound(compound_index, neighborinfo_node, statebuffer_node);
	for (int i = 0; i < water.n_atoms; i++) {
		compound.particles[i] = CompactParticle(com + water.atoms[i].pos, water.atoms[i].mass);
	}
	return Compound_H2O();
}

bool BoxBuilder::spaceAvailable(Float3 com, float radius) {	// Too lazy to implement yet..
	return true;
}

