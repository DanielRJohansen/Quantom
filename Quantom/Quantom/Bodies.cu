#include "Bodies.cuh"


Molecule::Molecule() {	// Always returns a h2o molecule rn
	n_atoms = 3;
	atoms = new Atom[3];
	atoms[0] = Atom(Double3(0, 0, 0), 0.048, 15.999);
	atoms[1] = Atom(Double3(95.7 / 1000.f, 0, 0), 0.053, 1.008);
	atoms[2] = Atom(Double3(-23.96 / 1000.f, 92.65/1000, 0), 0.053, 1.008);// was 0.53

	// Calc com
	CoM = Double3(0, 0, 0);
	//std::printf("Atom1: %f, %f, %f\n", atoms[, CoM.y, CoM.z);

	double accumulated_mass = 0;
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

