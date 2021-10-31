#include "BoxBuilder.cuh"





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


	simulation->box->compounds = new Compound_H2O[1'000'000];


	int solvate_cnt = 0;
	for (int z_index = 0; z_index < bodies_per_dim; z_index++) {
		for (int y_index = 0; y_index < bodies_per_dim; y_index++) {
			for (int x_index = 0; x_index < bodies_per_dim; x_index++) {
				if (solvate_cnt == N_BODIES_START)
					break;
				Float3 random = get3Random(10000);

				Float3 compound_base_pos = Float3(base + dist_between_compounds * (float)x_index, base + dist_between_compounds * (float)y_index, base + dist_between_compounds * (float)z_index);
				int compound_first_index = simulation->box->n_particles;
				Float3 compound_united_vel = Float3(random.x, random.y, random.z).norm() * mean_velocity;
				
				if (spaceAvailable(compound_base_pos, 0.2)) {
					Compound_H2O compound = createCompound(compound_base_pos);
					placeCompound(compound);
					solvate_cnt++;
				}
			}
		}
	}
	return solvate_cnt;
}

Compound_H2O BoxBuilder::createCompound(Float3 com)
{
	return Compound_H2O();
}
