#include "Bodies.cuh"


Molecule::Molecule() {	// Always returns a h2o molecule rn
	n_atoms = 3;
	atoms = new Atom[3];
	uint8_t white[3] = { 250, 250, 250 };
	uint8_t red[3] = { 250, 0, 0 };

	// Using Van der Waals radii.... https://en.m.wikipedia.org/wiki/Van_der_Waals_radius

	atoms[0] = Atom(Float3(0, 0, 0), 0.152, 15.999, red);
	atoms[1] = Atom(Float3(95.7 / 1000.f, 0, 0), 0.110, 1.008, white);
	atoms[2] = Atom(Float3(-23.96 / 1000.f, 92.65/1000, 0), 0.110, 1.008, white);// was 0.53
			
	for (int i = 0; i < 3; i++) {
		printf("Color %d:  %d %d %d\n", i, atoms[i].color[0], atoms[i].color[1], atoms[i].color[2]);
	}



	// Calc com
	CoM = Float3(0, 0, 0);
	//std::printf("Atom1: %f, %f, %f\n", atoms[, CoM.y, CoM.z);

	float accumulated_mass = 0;
	for (int i = 0; i < 3; i++) {
		CoM = CoM + (atoms[i].pos * atoms[i].mass);
		accumulated_mass += atoms[i].mass;
	}
	CoM = CoM * (1 / accumulated_mass);

	//printf("Mol CoM: %.4f %.4f %.4f\n", CoM.x, CoM.y, CoM.z);


	//exit(1);

	for (int i = 0; i < 3; i++) {
		atoms[i].pos = atoms[i].pos - CoM;
	}
	//std::printf("Center of mass: %f, %f, %f\n", CoM.x, CoM.y, CoM.z);

}

	